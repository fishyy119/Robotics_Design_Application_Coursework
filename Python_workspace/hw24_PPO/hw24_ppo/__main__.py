import argparse
from pathlib import Path
from typing import Mapping, cast

import torch

from .config import (
    EvalConfig,
    ModelConfig,
    PlotConfig,
    TrainConfig,
    deserialize_model_config,
    deserialize_train_config,
)
from .evaluate import run_evaluation
from .plot import run_plot
from .storage import load_checkpoint
from .train import run_training, validate_resume_checkpoint
from .video_grid import run_video_grid


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m hw24_ppo")
    subparsers = parser.add_subparsers(dest="command", required=True)

    default_eval_config = EvalConfig()

    train_parser = subparsers.add_parser("train", help="Train a PPO agent.")
    train_parser.add_argument("--run-name")
    train_parser.add_argument("--resume-from")
    train_parser.add_argument("--seed", type=int)
    train_parser.add_argument("--device")
    train_parser.add_argument("--epochs", type=int)
    train_parser.add_argument(
        "--steps-per-epoch",
        type=int,
    )
    train_parser.add_argument(
        "--hidden-sizes",
        type=int,
        nargs="+",
    )

    eval_parser = subparsers.add_parser("eval", help="Evaluate a saved PPO agent.")
    eval_parser.add_argument("--checkpoint", required=True)
    eval_parser.add_argument("--episodes", type=int, default=default_eval_config.episodes)
    eval_parser.add_argument("--device")
    eval_parser.add_argument("--render", action="store_true")
    eval_parser.add_argument("--stochastic", action="store_true")

    plot_parser = subparsers.add_parser("plot", help="Plot training curves from CSV.")
    plot_parser.add_argument("--metrics-csv", required=True)
    plot_parser.add_argument("--output")

    stitch_parser = subparsers.add_parser("stitch-videos", help="Stitch matched videos into a grid with ffmpeg.")
    stitch_parser.add_argument("--prefix", required=True)
    stitch_parser.add_argument("--grid-size", type=int, nargs=2, metavar=("ROWS", "COLS"), required=True)
    stitch_parser.add_argument("--output")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    default_model_config = ModelConfig()
    default_train_config = TrainConfig()

    if args.command == "train":
        resume_path = Path(str(args.resume_from)) if args.resume_from is not None else None
        model_config, train_config, run_name = resolve_training_inputs(
            run_name_raw=args.run_name,
            resume_path=resume_path,
            seed_raw=args.seed,
            device_raw=args.device,
            epochs_raw=args.epochs,
            steps_per_epoch_raw=args.steps_per_epoch,
            hidden_sizes_raw=args.hidden_sizes,
            default_model_config=default_model_config,
            default_train_config=default_train_config,
        )
        paths = run_training(
            run_name=run_name,
            train_config=train_config,
            model_config=model_config,
            resume_from=resume_path,
        )
        print(f"Training artifacts saved to: {paths.run_dir}")
        return

    if args.command == "eval":
        eval_config = EvalConfig(
            episodes=int(args.episodes),
            deterministic=not bool(args.stochastic),
            render=bool(args.render),
        )
        eval_artifacts = run_evaluation(
            checkpoint_path=Path(str(args.checkpoint)),
            eval_config=eval_config,
            device_override=str(args.device) if args.device is not None else None,
        )
        print(f"Evaluation CSV saved to: {eval_artifacts.eval_csv}")
        return

    if args.command == "plot":
        output_path = Path(str(args.output)) if args.output is not None else None
        plot_artifacts = run_plot(
            metrics_csv=Path(str(args.metrics_csv)),
            output=output_path,
            config=PlotConfig(),
        )
        print(f"Training curves saved to: {plot_artifacts.overview_png}")
        print(f"Return-only curve saved to: {plot_artifacts.return_only_png}")
        return

    if args.command == "stitch-videos":
        output_path = Path(str(args.output)) if args.output is not None else None
        stitch_artifacts = run_video_grid(
            prefix=Path(str(args.prefix)),
            grid_size=(int(args.grid_size[0]), int(args.grid_size[1])),
            output=output_path,
        )
        print(f"Stitched video saved to: {stitch_artifacts.output_video}")
        print(f"Collected {len(stitch_artifacts.input_videos)} videos.")
        return

    raise ValueError(f"Unsupported command: {args.command}")


def resolve_training_inputs(
    run_name_raw: str | None,
    resume_path: Path | None,
    seed_raw: int | None,
    device_raw: str | None,
    epochs_raw: int | None,
    steps_per_epoch_raw: int | None,
    hidden_sizes_raw: list[int] | None,
    default_model_config: ModelConfig,
    default_train_config: TrainConfig,
) -> tuple[ModelConfig, TrainConfig, str]:
    if resume_path is None:
        if run_name_raw is None:
            raise ValueError("--run-name is required when starting a new training run.")

        hidden_sizes = (
            tuple(int(size) for size in hidden_sizes_raw)
            if hidden_sizes_raw is not None
            else default_model_config.hidden_sizes
        )
        model_config = ModelConfig(hidden_sizes=hidden_sizes)
        train_config = TrainConfig(
            seed=seed_raw if seed_raw is not None else default_train_config.seed,
            device=str(device_raw) if device_raw is not None else default_train_config.device,
            epochs=epochs_raw if epochs_raw is not None else default_train_config.epochs,
            steps_per_epoch=(
                steps_per_epoch_raw if steps_per_epoch_raw is not None else default_train_config.steps_per_epoch
            ),
        )
        return model_config, train_config, str(run_name_raw)

    checkpoint = load_checkpoint(resume_path, torch.device("cpu"))
    validate_resume_checkpoint(checkpoint)
    checkpoint_model_config = deserialize_model_config(require_mapping_field(checkpoint, "model_config"))
    checkpoint_train_config = deserialize_train_config(require_mapping_field(checkpoint, "train_config"))

    hidden_sizes = (
        tuple(int(size) for size in hidden_sizes_raw)
        if hidden_sizes_raw is not None
        else checkpoint_model_config.hidden_sizes
    )
    model_config = ModelConfig(
        hidden_sizes=hidden_sizes,
        activation=checkpoint_model_config.activation,
    )
    train_config = TrainConfig(
        env_id=checkpoint_train_config.env_id,
        seed=seed_raw if seed_raw is not None else checkpoint_train_config.seed,
        device=str(device_raw) if device_raw is not None else checkpoint_train_config.device,
        epochs=epochs_raw if epochs_raw is not None else default_train_config.epochs,
        steps_per_epoch=(
            steps_per_epoch_raw if steps_per_epoch_raw is not None else checkpoint_train_config.steps_per_epoch
        ),
        gamma=checkpoint_train_config.gamma,
        lam=checkpoint_train_config.lam,
        clip_ratio=checkpoint_train_config.clip_ratio,
        pi_lr=checkpoint_train_config.pi_lr,
        vf_lr=checkpoint_train_config.vf_lr,
        train_pi_iters=checkpoint_train_config.train_pi_iters,
        train_v_iters=checkpoint_train_config.train_v_iters,
        target_kl=checkpoint_train_config.target_kl,
        max_ep_len=checkpoint_train_config.max_ep_len,
        timeout_penalty=checkpoint_train_config.timeout_penalty,
        use_curriculum=checkpoint_train_config.use_curriculum,
        curriculum_heights=checkpoint_train_config.curriculum_heights,
        curriculum_progress_rewards=checkpoint_train_config.curriculum_progress_rewards,
        curriculum_x_tolerance_scale=checkpoint_train_config.curriculum_x_tolerance_scale,
        curriculum_angle_tolerance_scale=checkpoint_train_config.curriculum_angle_tolerance_scale,
        curriculum_horizontal_speed_tolerance=(checkpoint_train_config.curriculum_horizontal_speed_tolerance),
        save_every=checkpoint_train_config.save_every,
    )
    run_name = str(run_name_raw) if run_name_raw is not None else resume_path.resolve().parent.name
    return model_config, train_config, run_name


def require_mapping_field(
    data: Mapping[str, object],
    key: str,
) -> Mapping[str, object]:
    value = data[key]
    if not isinstance(value, Mapping):
        raise ValueError(f"Checkpoint field {key} must be a mapping.")
    return cast(Mapping[str, object], value)


if __name__ == "__main__":
    main()
