from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ArtifactPaths:
    run_dir: Path
    config_json: Path
    metrics_csv: Path
    checkpoint_latest: Path
    checkpoint_best: Path
    eval_csv: Path
    training_curves_png: Path


def get_project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def get_artifacts_root() -> Path:
    return get_project_root() / "artifacts"


def get_run_artifact_paths(run_name: str) -> ArtifactPaths:
    run_dir = get_artifacts_root() / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    return ArtifactPaths(
        run_dir=run_dir,
        config_json=run_dir / "config.json",
        metrics_csv=run_dir / "metrics.csv",
        checkpoint_latest=run_dir / "checkpoint_latest.pt",
        checkpoint_best=run_dir / "checkpoint_best.pt",
        eval_csv=run_dir / "eval.csv",
        training_curves_png=run_dir / "training_curves.png",
    )


def resolve_plot_output(metrics_csv: Path, output: Path | None) -> Path:
    if output is not None:
        return output
    return metrics_csv.with_name("training_curves.png")
