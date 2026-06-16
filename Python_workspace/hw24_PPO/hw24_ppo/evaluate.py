# pyright: standard

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, cast

import gymnasium as gym
import imageio.v2 as imageio
import numpy as np
import torch
from numpy.typing import ArrayLike

from .config import EvalConfig, deserialize_model_config, deserialize_train_config
from .env import make_environment, to_observation_array, validate_environment
from .networks import MLPActorCritic
from .storage import append_csv_row, initialize_csv, load_checkpoint
from .types import ObsArray, ObsBatchTensor

EVAL_FIELDNAMES = ["episode", "return", "length"]


@dataclass(frozen=True)
class EvaluationArtifacts:
    eval_csv: Path
    video_paths: tuple[Path, ...]


def run_evaluation(
    checkpoint_path: Path,
    eval_config: EvalConfig,
    device_override: str | None = None,
) -> EvaluationArtifacts:
    checkpoint_dir = checkpoint_path.resolve().parent
    load_device = torch.device(device_override) if device_override is not None else torch.device("cpu")
    checkpoint = load_checkpoint(checkpoint_path, load_device)
    validate_checkpoint(checkpoint)

    model_config = deserialize_model_config(cast(Mapping[str, object], checkpoint["model_config"]))
    train_config = deserialize_train_config(cast(Mapping[str, object], checkpoint["train_config"]))
    runtime_device = torch.device(device_override or "cpu")

    record_env = make_environment(train_config.env_id, render_mode="rgb_array")
    display_env = make_environment(train_config.env_id, render_mode="human") if eval_config.render else None
    env_spec = validate_environment(record_env)
    model = MLPActorCritic(
        obs_dim=env_spec.obs_dim,
        act_dim=env_spec.act_dim,
        config=model_config,
    ).to(runtime_device)
    model.load_state_dict(cast(dict[str, object], checkpoint["model_state_dict"]))
    model.eval()

    eval_csv = checkpoint_dir / "eval.csv"
    initialize_csv(eval_csv, EVAL_FIELDNAMES)

    video_dir = checkpoint_dir / "eval_videos"
    video_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_name = checkpoint_path.stem
    policy_name = "deterministic" if eval_config.deterministic else "stochastic"
    fps = resolve_render_fps(record_env)

    returns: list[float] = []
    lengths: list[int] = []
    video_paths: list[Path] = []

    for episode_index in range(1, eval_config.episodes + 1):
        episode_seed = train_config.seed + episode_index - 1
        observation, _ = record_env.reset(seed=episode_seed)
        if display_env is not None:
            display_env.reset(seed=episode_seed)

        obs = to_observation_array(observation, env_spec.obs_dim)
        frames: list[np.ndarray] = [capture_frame(record_env)]

        done = False
        episode_return = 0.0
        episode_length = 0

        while not done:
            obs_tensor = observation_to_tensor(obs, runtime_device)
            if eval_config.deterministic:
                action_tensor = model.act_deterministic(obs_tensor)
            else:
                action_tensor, _, _ = model.act(obs_tensor)

            action = int(action_tensor[0].item())
            next_observation, reward, terminated, truncated, _ = record_env.step(action)
            if display_env is not None:
                display_env.step(action)

            frames.append(capture_frame(record_env))
            obs = to_observation_array(next_observation, env_spec.obs_dim)
            episode_return += float(reward)
            episode_length += 1
            done = bool(terminated or truncated or episode_length >= train_config.max_ep_len)

        returns.append(episode_return)
        lengths.append(episode_length)
        append_csv_row(
            eval_csv,
            EVAL_FIELDNAMES,
            {
                "episode": episode_index,
                "return": episode_return,
                "length": episode_length,
            },
        )

        video_path = build_episode_video_path(
            video_dir=video_dir,
            checkpoint_name=checkpoint_name,
            policy_name=policy_name,
            episode_index=episode_index,
            episode_total=eval_config.episodes,
        )
        write_video(video_path=video_path, frames=frames, fps=fps)
        video_paths.append(video_path)
        print(f"Saved evaluation video: {video_path}")

    record_env.close()
    if display_env is not None:
        display_env.close()

    mean_return = float(np.mean(returns))
    std_return = float(np.std(returns))
    min_return = float(np.min(returns))
    max_return = float(np.max(returns))
    mean_length = float(np.mean(lengths))

    print(
        "Evaluation | "
        f"policy={policy_name} | "
        f"episodes={eval_config.episodes} | "
        f"mean_return={mean_return:.2f} | "
        f"std_return={std_return:.2f} | "
        f"min_return={min_return:.2f} | "
        f"max_return={max_return:.2f} | "
        f"mean_length={mean_length:.2f}"
    )
    print(f"Evaluation videos saved to: {video_dir}")

    return EvaluationArtifacts(
        eval_csv=eval_csv,
        video_paths=tuple(video_paths),
    )


def validate_checkpoint(checkpoint: Mapping[str, object]) -> None:
    required_fields = [
        "model_state_dict",
        "model_config",
        "train_config",
        "epoch",
        "global_step",
        "best_mean_return",
    ]
    for field in required_fields:
        if field not in checkpoint:
            raise ValueError(f"Checkpoint is missing required field: {field}")


def observation_to_tensor(obs: ObsArray, device: torch.device) -> ObsBatchTensor:
    return torch.as_tensor(
        np.expand_dims(obs, axis=0),
        dtype=torch.float32,
        device=device,
    )


def capture_frame(env: gym.Env[object, object]) -> np.ndarray:
    frame = env.render()
    if frame is None:
        raise RuntimeError("RGB frame capture failed during evaluation.")
    array = np.asarray(frame, dtype=np.uint8)
    if array.ndim != 3:
        raise RuntimeError(f"Expected RGB frame with 3 dimensions, got {array.shape}.")
    return array


def resolve_render_fps(env: gym.Env[object, object]) -> int:
    metadata = getattr(env, "metadata", {})
    if isinstance(metadata, Mapping):
        fps_value = metadata.get("render_fps")
        if isinstance(fps_value, int) and fps_value > 0:
            return fps_value
    return 50


def build_episode_video_path(
    video_dir: Path,
    checkpoint_name: str,
    policy_name: str,
    episode_index: int,
    episode_total: int,
) -> Path:
    return video_dir / (f"{checkpoint_name}_{policy_name}_" f"episode_{episode_index:02d}_of_{episode_total:02d}.mp4")


def write_video(
    video_path: Path,
    frames: list[np.ndarray],
    fps: int,
) -> None:
    compatible_frames = [pad_frame_for_video(frame) for frame in frames]
    imageio.mimsave(video_path, cast(list[ArrayLike], compatible_frames), fps=fps)


def pad_frame_for_video(frame: np.ndarray) -> np.ndarray:
    height, width, _ = frame.shape
    padded_height = round_up_to_multiple(height, 16)
    padded_width = round_up_to_multiple(width, 16)
    if padded_height == height and padded_width == width:
        return frame

    pad_height = padded_height - height
    pad_width = padded_width - width
    return np.pad(
        frame,
        ((0, pad_height), (0, pad_width), (0, 0)),
        mode="edge",
    )


def round_up_to_multiple(value: int, multiple: int) -> int:
    remainder = value % multiple
    if remainder == 0:
        return value
    return value + multiple - remainder
