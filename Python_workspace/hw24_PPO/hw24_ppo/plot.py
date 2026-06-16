# pyright: standard

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

from .config import PlotConfig
from .paths import resolve_plot_output

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402


@dataclass(frozen=True)
class PlotArtifacts:
    overview_png: Path
    return_only_png: Path


def run_plot(
    metrics_csv: Path,
    output: Path | None = None,
    config: PlotConfig | None = None,
) -> PlotArtifacts:
    del config
    if not metrics_csv.exists():
        raise FileNotFoundError(f"Metrics CSV not found: {metrics_csv}")

    rows = load_metrics(metrics_csv)
    if not rows:
        raise ValueError("Metrics CSV is empty.")

    epochs = np.array([row["epoch"] for row in rows], dtype=float)
    mean_returns = np.array([row["mean_episode_return"] for row in rows], dtype=float)
    mean_lengths = np.array([row["mean_episode_length"] for row in rows], dtype=float)
    policy_loss = np.array([row["policy_loss"] for row in rows], dtype=float)
    value_loss = np.array([row["value_loss"] for row in rows], dtype=float)
    approx_kl = np.array([row["approx_kl"] for row in rows], dtype=float)
    entropy = np.array([row["entropy"] for row in rows], dtype=float)
    clip_fraction = np.array([row["clip_fraction"] for row in rows], dtype=float)

    configure_plot_style()
    output_path, return_only_path = resolve_plot_outputs(metrics_csv, output)

    figure, axes_grid = plt.subplots(2, 2, figsize=(12, 8))
    axes = np.asarray(axes_grid, dtype=object).reshape(-1).tolist()

    axes[0].plot(epochs, mean_returns)
    axes[0].set_title("Mean Episode Return")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Return")
    axes[0].grid(True, linestyle="--", linewidth=0.6, alpha=0.7)

    axes[1].plot(epochs, mean_lengths)
    axes[1].set_title("Mean Episode Length")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Length")
    axes[1].grid(True, linestyle="--", linewidth=0.6, alpha=0.7)

    axes[2].plot(epochs, policy_loss, label="Policy Loss")
    axes[2].plot(epochs, value_loss, label="Value Loss")
    axes[2].set_title("Losses")
    axes[2].set_xlabel("Epoch")
    axes[2].legend()
    axes[2].grid(True, linestyle="--", linewidth=0.6, alpha=0.7)

    axes[3].plot(epochs, approx_kl, label="Approx KL")
    axes[3].plot(epochs, entropy, label="Entropy")
    axes[3].plot(epochs, clip_fraction, label="Clip Fraction")
    axes[3].set_title("Policy Diagnostics")
    axes[3].set_xlabel("Epoch")
    axes[3].legend()
    axes[3].grid(True, linestyle="--", linewidth=0.6, alpha=0.7)

    figure.tight_layout()
    figure.savefig(output_path, dpi=200)
    plt.close(figure)

    return_only_figure, return_only_axis = plt.subplots(1, 1, figsize=(8, 5))
    return_only_axis.plot(epochs, mean_returns)
    return_only_axis.set_title("Mean Episode Return")
    return_only_axis.set_xlabel("Epoch")
    return_only_axis.set_ylabel("Return")
    return_only_axis.set_ylim(top=300)
    return_only_axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.7)

    return_only_figure.tight_layout()
    return_only_figure.savefig(return_only_path, dpi=200)
    plt.close(return_only_figure)

    return PlotArtifacts(
        overview_png=output_path,
        return_only_png=return_only_path,
    )


def load_metrics(metrics_csv: Path) -> list[dict[str, float]]:
    with metrics_csv.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)
        rows: list[dict[str, float]] = []
        for row in reader:
            rows.append({key: float(value) for key, value in row.items() if key is not None and value is not None})
        return rows


def configure_plot_style() -> None:
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.serif"] = ["Times New Roman"]


def resolve_plot_outputs(metrics_csv: Path, output: Path | None) -> tuple[Path, Path]:
    overview_path = resolve_plot_output(metrics_csv, output)
    suffix = overview_path.suffix if overview_path.suffix else ".png"
    return_only_path = overview_path.with_name(f"{overview_path.stem}_return_only{suffix}")
    return overview_path, return_only_path
