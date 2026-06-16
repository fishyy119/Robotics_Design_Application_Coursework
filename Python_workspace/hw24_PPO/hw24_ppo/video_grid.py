import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, cast

VIDEO_SUFFIXES = (".mp4", ".mov", ".mkv", ".avi", ".webm")


@dataclass(frozen=True)
class VideoInfo:
    path: Path
    duration_seconds: float
    width: int
    height: int


@dataclass(frozen=True)
class VideoGridArtifacts:
    output_video: Path
    input_videos: tuple[Path, ...]


def run_video_grid(
    prefix: Path,
    grid_size: tuple[int, int],
    output: Path | None = None,
) -> VideoGridArtifacts:
    grid_rows, grid_cols = grid_size
    if grid_rows <= 0 or grid_cols <= 0:
        raise ValueError("Grid size must be positive.")

    output_path = resolve_output_path(prefix, grid_rows, grid_cols, output)
    input_videos = collect_matching_videos(prefix, output_path)
    if len(input_videos) > grid_rows * grid_cols:
        print(
            f"Grid {grid_rows}x{grid_cols} can place at most {grid_rows * grid_cols} videos, "
            f"but found {len(input_videos)}."
        )
        print(f"Using the first {len(input_videos)} videos for the grid.")
        input_videos = input_videos[: grid_rows * grid_cols]

    ffmpeg_path = require_binary("ffmpeg")
    ffprobe_path = require_binary("ffprobe")
    video_infos = tuple(probe_video(video_path, ffprobe_path) for video_path in input_videos)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    command = build_ffmpeg_command(
        ffmpeg_path=ffmpeg_path,
        video_infos=video_infos,
        grid_rows=grid_rows,
        grid_cols=grid_cols,
        output_path=output_path,
    )
    subprocess.run(command, check=True)

    return VideoGridArtifacts(
        output_video=output_path,
        input_videos=input_videos,
    )


def resolve_output_path(
    prefix: Path,
    grid_rows: int,
    grid_cols: int,
    output: Path | None,
) -> Path:
    if output is not None:
        return output

    search_dir, file_prefix = split_prefix(prefix)
    prefix_stem = file_prefix if file_prefix else search_dir.name
    return search_dir / f"{prefix_stem}_grid_{grid_rows}x{grid_cols}.mp4"


def collect_matching_videos(prefix: Path, output_path: Path) -> tuple[Path, ...]:
    search_dir, file_prefix = split_prefix(prefix)
    if not search_dir.exists():
        raise FileNotFoundError(f"Search directory not found: {search_dir}")

    resolved_output = output_path.resolve()
    matched = [
        candidate
        for candidate in search_dir.iterdir()
        if candidate.is_file()
        and candidate.suffix.lower() in VIDEO_SUFFIXES
        and candidate.name.startswith(file_prefix)
        and candidate.resolve() != resolved_output
    ]
    matched.sort(key=lambda path: natural_sort_key(path.name))

    if not matched:
        raise FileNotFoundError(f"No video files found for prefix: {prefix}")

    return tuple(matched)


def split_prefix(prefix: Path) -> tuple[Path, str]:
    resolved_prefix = prefix.expanduser()
    if resolved_prefix.exists() and resolved_prefix.is_dir():
        return resolved_prefix, ""

    return resolved_prefix.parent if resolved_prefix.parent != Path("") else Path("."), resolved_prefix.name


def natural_sort_key(text: str) -> tuple[tuple[int, str | int], ...]:
    parts = re.split(r"(\d+)", text.lower())
    key_parts: list[tuple[int, str | int]] = []
    for part in parts:
        if part.isdigit():
            key_parts.append((1, int(part)))
        else:
            key_parts.append((0, part))
    return tuple(key_parts)


def require_binary(name: str) -> Path:
    binary = shutil.which(name)
    if binary is None:
        raise FileNotFoundError(f"Required binary not found in PATH: {name}")
    return Path(binary)


def probe_video(video_path: Path, ffprobe_path: Path) -> VideoInfo:
    command = [
        str(ffprobe_path),
        "-v",
        "error",
        "-select_streams",
        "v:0",
        "-show_entries",
        "stream=width,height:format=duration",
        "-of",
        "json",
        str(video_path),
    ]
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = parse_json_object(result.stdout, context=f"ffprobe output for {video_path}")

    streams = require_list_field(payload, "streams")
    if not streams:
        raise ValueError(f"No video stream found in: {video_path}")

    stream_info = require_mapping(streams[0], context=f"first stream in {video_path}")
    format_info = require_mapping(payload["format"], context=f"format section in {video_path}")

    width = require_int_field(stream_info, "width")
    height = require_int_field(stream_info, "height")
    duration_seconds = require_float_field(format_info, "duration")

    if width <= 0 or height <= 0:
        raise ValueError(f"Invalid video size for {video_path}: {width}x{height}")
    if duration_seconds <= 0.0:
        raise ValueError(f"Invalid video duration for {video_path}: {duration_seconds}")

    return VideoInfo(
        path=video_path,
        duration_seconds=duration_seconds,
        width=width,
        height=height,
    )


def build_ffmpeg_command(
    ffmpeg_path: Path,
    video_infos: tuple[VideoInfo, ...],
    grid_rows: int,
    grid_cols: int,
    output_path: Path,
) -> list[str]:
    if not video_infos:
        raise ValueError("At least one input video is required.")

    tile_width = max(info.width for info in video_infos)
    tile_height = max(info.height for info in video_infos)
    max_duration = max(info.duration_seconds for info in video_infos)

    command = [str(ffmpeg_path), "-y"]
    for info in video_infos:
        command.extend(["-i", str(info.path)])

    filter_complex = build_filter_complex(
        video_infos=video_infos,
        tile_width=tile_width,
        tile_height=tile_height,
        max_duration=max_duration,
        grid_rows=grid_rows,
        grid_cols=grid_cols,
    )
    command.extend(
        [
            "-filter_complex",
            filter_complex,
            "-map",
            "[vout]",
            "-an",
            "-c:v",
            "libx264",
            "-crf",
            "18",
            "-preset",
            "medium",
            "-movflags",
            "+faststart",
            str(output_path),
        ]
    )
    return command


def build_filter_complex(
    video_infos: tuple[VideoInfo, ...],
    tile_width: int,
    tile_height: int,
    max_duration: float,
    grid_rows: int,
    grid_cols: int,
) -> str:
    filter_parts: list[str] = []

    for index, info in enumerate(video_infos):
        stop_duration = max(0.0, max_duration - info.duration_seconds)
        filters = ["setpts=PTS-STARTPTS"]
        if stop_duration > 1e-6:
            filters.append(f"tpad=stop_mode=clone:stop_duration={stop_duration:.6f}")
        filters.extend(
            [
                f"scale={tile_width}:{tile_height}:force_original_aspect_ratio=decrease",
                f"pad={tile_width}:{tile_height}:(ow-iw)/2:(oh-ih)/2:black",
                "setsar=1",
                "fps=30",
            ]
        )
        filter_parts.append(f"[{index}:v]{','.join(filters)}[v{index}]")

    canvas_width = grid_cols * tile_width
    canvas_height = grid_rows * tile_height

    if len(video_infos) == 1:
        filter_parts.append(f"[v0]pad={canvas_width}:{canvas_height}:0:0:black,format=yuv420p[vout]")
        return ";".join(filter_parts)

    stacked_inputs = "".join(f"[v{index}]" for index in range(len(video_infos)))
    layout = "|".join(
        f"{(index % grid_cols) * tile_width}_{(index // grid_cols) * tile_height}" for index in range(len(video_infos))
    )
    filter_parts.append(f"{stacked_inputs}xstack=inputs={len(video_infos)}:layout={layout}[stacked]")
    filter_parts.append(f"[stacked]pad={canvas_width}:{canvas_height}:0:0:black,format=yuv420p[vout]")
    return ";".join(filter_parts)


def parse_json_object(raw_text: str, context: str) -> Mapping[str, object]:
    parsed: object = json.loads(raw_text)
    return require_mapping(parsed, context=context)


def require_mapping(value: object, context: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Expected JSON object for {context}.")
    return cast(Mapping[str, object], value)


def require_list_field(data: Mapping[str, object], key: str) -> list[object]:
    value = data.get(key)
    if not isinstance(value, list):
        raise ValueError(f"Expected list field {key}.")
    return value


def require_int_field(data: Mapping[str, object], key: str) -> int:
    value = data.get(key)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Expected integer field {key}.")
    return value


def require_float_field(data: Mapping[str, object], key: str) -> float:
    value = data.get(key)
    if isinstance(value, bool):
        raise ValueError(f"Expected numeric field {key}.")
    if isinstance(value, int | float):
        return float(value)
    if isinstance(value, str):
        return float(value)
    raise ValueError(f"Expected numeric field {key}.")
