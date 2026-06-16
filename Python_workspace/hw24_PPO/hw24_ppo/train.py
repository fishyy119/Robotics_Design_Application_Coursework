# pyright: standard

import math
import random
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np
import torch
from torch.optim import Adam, Optimizer

from .buffer import RolloutBatch, RolloutBuffer
from .config import (
    ModelConfig,
    TrainConfig,
    build_config_payload,
    deserialize_model_config,
    serialize_model_config,
    serialize_train_config,
)
from .curriculum import CurriculumStepResult, build_curriculum
from .env import (
    RewardBreakdown,
    add_reward_breakdowns,
    make_environment,
    make_reward_tracker,
    mean_reward_breakdown,
    to_observation_array,
    validate_environment,
)
from .networks import MLPActorCritic
from .paths import ArtifactPaths, get_run_artifact_paths
from .storage import (
    append_csv_row,
    initialize_csv,
    load_checkpoint,
    save_checkpoint,
    write_json,
)
from .types import ObsArray, ObsBatchTensor

METRICS_FIELDNAMES = [
    "epoch",
    "env_steps",
    "episodes_finished",
    "mean_episode_return",
    "mean_episode_length",
    "policy_loss",
    "value_loss",
    "approx_kl",
    "entropy",
    "clip_fraction",
    "elapsed_seconds",
]


@dataclass(frozen=True)
class UpdateMetrics:
    policy_loss: float
    value_loss: float
    approx_kl: float
    entropy: float
    clip_fraction: float


def run_training(
    run_name: str,
    train_config: TrainConfig,
    model_config: ModelConfig,
    resume_from: Path | None = None,
) -> ArtifactPaths:
    paths = get_run_artifact_paths(run_name)
    seed_everything(train_config.seed)
    device = torch.device(train_config.device)

    env = make_environment(train_config.env_id)
    env_spec = validate_environment(env)
    env.action_space.seed(train_config.seed)
    reward_tracker = make_reward_tracker(train_config.env_id)
    curriculum = build_curriculum(train_config)

    model = MLPActorCritic(
        obs_dim=env_spec.obs_dim,
        act_dim=env_spec.act_dim,
        config=model_config,
    ).to(device)
    policy_optimizer = Adam(model.actor.parameters(), lr=train_config.pi_lr)
    value_optimizer = Adam(model.critic.parameters(), lr=train_config.vf_lr)
    buffer = RolloutBuffer(
        obs_dim=env_spec.obs_dim,
        size=train_config.steps_per_epoch,
        gamma=train_config.gamma,
        lam=train_config.lam,
    )

    if resume_from is not None:
        start_epoch, global_step, best_mean_return = restore_training_state(
            checkpoint_path=resume_from,
            model=model,
            policy_optimizer=policy_optimizer,
            value_optimizer=value_optimizer,
            expected_model_config=model_config,
            runtime_device=device,
        )
    else:
        start_epoch = 1
        global_step = 0
        best_mean_return = -math.inf

    write_training_metadata(
        paths=paths,
        run_name=run_name,
        model_config=model_config,
        train_config=train_config,
        resume_from=resume_from,
    )
    if not paths.metrics_csv.exists() or resume_from is None:
        initialize_csv(paths.metrics_csv, METRICS_FIELDNAMES)

    last_checkpoint: dict[str, object] | None = None
    start_time = time.time()

    observation, _ = env.reset(seed=train_config.seed)
    obs = to_observation_array(observation, env_spec.obs_dim)
    if reward_tracker is not None:
        reward_tracker.reset()
    if curriculum is not None:
        curriculum.reset_episode()
    episode_return = 0.0
    episode_length = 0
    episode_reward_breakdown = RewardBreakdown()

    end_epoch = start_epoch + train_config.epochs - 1
    for epoch in range(start_epoch, end_epoch + 1):
        epoch_returns: list[float] = []
        epoch_lengths: list[int] = []
        epoch_reward_breakdowns: list[RewardBreakdown] | None = [] if reward_tracker is not None else None

        for step in range(train_config.steps_per_epoch):
            obs_tensor = observation_to_tensor(obs, device)
            action_tensor, value_tensor, log_prob_tensor = model.act(obs_tensor)

            action = int(action_tensor[0].item())
            value = float(value_tensor[0].item())
            log_prob = float(log_prob_tensor[0].item())

            next_observation, reward, terminated, truncated, _ = env.step(action)
            next_obs = to_observation_array(next_observation, env_spec.obs_dim)
            episode_length += 1
            timeout = episode_length >= train_config.max_ep_len
            timeout_ended = bool((truncated or timeout) and not terminated)
            step_reward = float(reward)
            curriculum_step: CurriculumStepResult | None = None

            if curriculum is not None:
                curriculum_step = curriculum.evaluate_step(
                    observation=next_obs,
                    terminated=bool(terminated),
                )
                step_reward += curriculum_step.progress_reward

            if timeout_ended:
                step_reward += train_config.timeout_penalty

            if reward_tracker is not None:
                step_breakdown = reward_tracker.step(
                    action=action,
                    observation=next_obs,
                    reward=float(reward),
                    terminated=bool(terminated),
                )
                if curriculum_step is not None and curriculum_step.progress_reward != 0.0:
                    step_breakdown = add_reward_breakdowns(
                        step_breakdown,
                        RewardBreakdown(curriculum_progress=curriculum_step.progress_reward),
                    )
                if timeout_ended:
                    step_breakdown = add_reward_breakdowns(
                        step_breakdown,
                        RewardBreakdown(timeout=train_config.timeout_penalty),
                    )
                episode_reward_breakdown = add_reward_breakdowns(
                    episode_reward_breakdown,
                    step_breakdown,
                )

            episode_return += step_reward
            global_step += 1

            buffer.store(
                obs=obs,
                action=action,
                reward=step_reward,
                value=value,
                log_prob=log_prob,
            )

            obs = next_obs
            episode_done = bool(terminated or truncated or timeout)
            epoch_ended = step == train_config.steps_per_epoch - 1

            if episode_done or epoch_ended:
                if terminated:
                    last_value = 0.0
                else:
                    last_value_tensor = model.critic(observation_to_tensor(obs, device))
                    last_value = float(last_value_tensor[0].item())

                buffer.finish_path(last_value=last_value)

                if episode_done:
                    epoch_returns.append(episode_return)
                    epoch_lengths.append(episode_length)
                    if epoch_reward_breakdowns is not None:
                        epoch_reward_breakdowns.append(episode_reward_breakdown)
                    reset_observation, _ = env.reset()
                    obs = to_observation_array(reset_observation, env_spec.obs_dim)
                    if reward_tracker is not None:
                        reward_tracker.reset()
                    if curriculum is not None:
                        curriculum.reset_episode()
                    episode_return = 0.0
                    episode_length = 0
                    episode_reward_breakdown = RewardBreakdown()

        batch = buffer.get_tensors(device)
        update_metrics = update_model(
            model=model,
            batch=batch,
            policy_optimizer=policy_optimizer,
            value_optimizer=value_optimizer,
            train_config=train_config,
        )

        mean_episode_return = safe_mean(epoch_returns)
        mean_episode_length = safe_mean(epoch_lengths)
        elapsed_seconds = float(time.time() - start_time)

        row = {
            "epoch": epoch,
            "env_steps": global_step,
            "episodes_finished": len(epoch_returns),
            "mean_episode_return": mean_episode_return,
            "mean_episode_length": mean_episode_length,
            "policy_loss": update_metrics.policy_loss,
            "value_loss": update_metrics.value_loss,
            "approx_kl": update_metrics.approx_kl,
            "entropy": update_metrics.entropy,
            "clip_fraction": update_metrics.clip_fraction,
            "elapsed_seconds": elapsed_seconds,
        }
        append_csv_row(paths.metrics_csv, METRICS_FIELDNAMES, row)

        checkpoint = build_checkpoint(
            model=model,
            policy_optimizer=policy_optimizer,
            value_optimizer=value_optimizer,
            model_config=model_config,
            train_config=train_config,
            epoch=epoch,
            global_step=global_step,
            best_mean_return=best_mean_return,
        )
        last_checkpoint = checkpoint

        if not math.isnan(mean_episode_return) and mean_episode_return > best_mean_return:
            best_mean_return = mean_episode_return
            checkpoint["best_mean_return"] = best_mean_return
            save_checkpoint(paths.checkpoint_best, checkpoint)

        checkpoint["best_mean_return"] = best_mean_return
        if epoch % train_config.save_every == 0 or epoch == end_epoch:
            save_checkpoint(paths.checkpoint_latest, checkpoint)

        print(format_epoch_summary(epoch, row))
        print(format_reward_breakdown_summary(epoch_reward_breakdowns))

    if last_checkpoint is None:
        raise RuntimeError("Training finished without producing a checkpoint.")

    if not paths.checkpoint_latest.exists():
        save_checkpoint(paths.checkpoint_latest, last_checkpoint)
    if not paths.checkpoint_best.exists():
        save_checkpoint(paths.checkpoint_best, last_checkpoint)

    env.close()
    return paths


def restore_training_state(
    checkpoint_path: Path,
    model: MLPActorCritic,
    policy_optimizer: Optimizer,
    value_optimizer: Optimizer,
    expected_model_config: ModelConfig,
    runtime_device: torch.device,
) -> tuple[int, int, float]:
    checkpoint = load_checkpoint(checkpoint_path, runtime_device)
    validate_resume_checkpoint(checkpoint)

    checkpoint_model_config = deserialize_model_config(cast(dict[str, object], checkpoint["model_config"]))
    if checkpoint_model_config != expected_model_config:
        raise ValueError("Resume checkpoint model_config does not match the requested model configuration.")

    model.load_state_dict(cast(dict[str, object], checkpoint["model_state_dict"]))
    policy_optimizer.load_state_dict(cast(dict[str, object], checkpoint["optimizer_state_dict_pi"]))
    value_optimizer.load_state_dict(cast(dict[str, object], checkpoint["optimizer_state_dict_v"]))

    last_epoch = require_checkpoint_int(checkpoint, "epoch")
    global_step = require_checkpoint_int(checkpoint, "global_step")
    best_mean_return = require_checkpoint_float(checkpoint, "best_mean_return")
    return last_epoch + 1, global_step, best_mean_return


def write_training_metadata(
    paths: ArtifactPaths,
    run_name: str,
    model_config: ModelConfig,
    train_config: TrainConfig,
    resume_from: Path | None,
) -> None:
    payload = build_config_payload(
        run_name=run_name,
        model_config=model_config,
        train_config=train_config,
    )
    if resume_from is not None:
        payload["resume_from"] = str(resume_from.resolve())
    write_json(paths.config_json, payload)


def validate_resume_checkpoint(checkpoint: Mapping[str, object]) -> None:
    required_fields = [
        "model_state_dict",
        "optimizer_state_dict_pi",
        "optimizer_state_dict_v",
        "model_config",
        "train_config",
        "epoch",
        "global_step",
        "best_mean_return",
    ]
    for field in required_fields:
        if field not in checkpoint:
            raise ValueError(f"Resume checkpoint is missing required field: {field}.")


def require_checkpoint_int(checkpoint: Mapping[str, object], key: str) -> int:
    value = checkpoint[key]
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"Checkpoint field {key} must be an integer.")
    return value


def require_checkpoint_float(checkpoint: Mapping[str, object], key: str) -> float:
    value = checkpoint[key]
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"Checkpoint field {key} must be numeric.")
    return float(value)


def update_model(
    model: MLPActorCritic,
    batch: RolloutBatch,
    policy_optimizer: Optimizer,
    value_optimizer: Optimizer,
    train_config: TrainConfig,
) -> UpdateMetrics:
    last_policy_loss = 0.0
    last_value_loss = 0.0
    last_approx_kl = 0.0
    last_entropy = 0.0
    last_clip_fraction = 0.0

    for _ in range(train_config.train_pi_iters):
        policy_optimizer.zero_grad()
        policy_loss, approx_kl, entropy, clip_fraction = compute_policy_loss(
            model=model,
            batch=batch,
            clip_ratio=train_config.clip_ratio,
        )

        last_policy_loss = float(policy_loss.item())
        last_approx_kl = approx_kl
        last_entropy = entropy
        last_clip_fraction = clip_fraction

        if approx_kl > 1.5 * train_config.target_kl:
            break

        policy_loss.backward()
        policy_optimizer.step()

    for _ in range(train_config.train_v_iters):
        value_optimizer.zero_grad()
        value_loss = compute_value_loss(model=model, batch=batch)
        last_value_loss = float(value_loss.item())
        value_loss.backward()
        value_optimizer.step()

    return UpdateMetrics(
        policy_loss=last_policy_loss,
        value_loss=last_value_loss,
        approx_kl=last_approx_kl,
        entropy=last_entropy,
        clip_fraction=last_clip_fraction,
    )


def compute_policy_loss(
    model: MLPActorCritic,
    batch: RolloutBatch,
    clip_ratio: float,
) -> tuple[torch.Tensor, float, float, float]:
    log_prob, entropy, _ = model.evaluate_actions(batch.obs, batch.act)
    ratio = torch.exp(log_prob - batch.logp)
    unclipped = ratio * batch.adv
    clipped = torch.clamp(ratio, 1.0 - clip_ratio, 1.0 + clip_ratio) * batch.adv
    policy_loss = -torch.min(unclipped, clipped).mean()

    approx_kl = float((batch.logp - log_prob).mean().item())
    entropy_mean = float(entropy.mean().item())
    clipped_mask = (ratio > (1.0 + clip_ratio)) | (ratio < (1.0 - clip_ratio))
    clip_fraction = float(clipped_mask.float().mean().item())

    return policy_loss, approx_kl, entropy_mean, clip_fraction


def compute_value_loss(model: MLPActorCritic, batch: RolloutBatch) -> torch.Tensor:
    values = model.critic(batch.obs)
    return ((values - batch.ret) ** 2).mean()


def build_checkpoint(
    model: MLPActorCritic,
    policy_optimizer: Optimizer,
    value_optimizer: Optimizer,
    model_config: ModelConfig,
    train_config: TrainConfig,
    epoch: int,
    global_step: int,
    best_mean_return: float,
) -> dict[str, object]:
    return {
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict_pi": policy_optimizer.state_dict(),
        "optimizer_state_dict_v": value_optimizer.state_dict(),
        "model_config": serialize_model_config(model_config),
        "train_config": serialize_train_config(train_config),
        "epoch": epoch,
        "global_step": global_step,
        "best_mean_return": best_mean_return,
    }


def observation_to_tensor(obs: ObsArray, device: torch.device) -> ObsBatchTensor:
    obs_batch = np.expand_dims(obs, axis=0)
    return torch.as_tensor(obs_batch, dtype=torch.float32, device=device)


def safe_mean(values: list[float] | list[int]) -> float:
    if not values:
        return float("nan")
    return float(np.mean(values))


def format_epoch_summary(epoch: int, row: Mapping[str, int | float]) -> str:
    return (
        f"Epoch {epoch:03d} | "
        f"steps={row['env_steps']} | "
        f"episodes={row['episodes_finished']} | "
        f"return={format_metric(float(row['mean_episode_return']))} | "
        f"length={format_metric(float(row['mean_episode_length']))} | "
        f"pi={float(row['policy_loss']):.4f} | "
        f"v={float(row['value_loss']):.4f} | "
        f"kl={float(row['approx_kl']):.4f} | "
        f"ent={float(row['entropy']):.4f} | "
        f"clip={float(row['clip_fraction']):.4f}"
    )


def format_metric(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return f"{value:.2f}"


def format_reward_breakdown_summary(
    epoch_reward_breakdowns: list[RewardBreakdown] | None,
) -> str:
    if epoch_reward_breakdowns is None:
        return "Reward | breakdown=unavailable"

    mean_breakdown = mean_reward_breakdown(epoch_reward_breakdowns)
    if mean_breakdown is None:
        return (
            "Reward | "
            "pos=nan | vel=nan | angle=nan | legs=nan | "
            "main=nan | side=nan | curr_prog=nan | terminal=nan | timeout=nan"
        )

    return (
        "Reward | "
        f"pos={mean_breakdown.position:.2f} | "
        f"vel={mean_breakdown.velocity:.2f} | "
        f"angle={mean_breakdown.angle:.2f} | "
        f"legs={mean_breakdown.leg_contact:.2f} | "
        f"main={mean_breakdown.main_engine:.2f} | "
        f"side={mean_breakdown.side_engine:.2f} | "
        f"curr_prog={mean_breakdown.curriculum_progress:.2f} | "
        f"terminal={mean_breakdown.terminal:.2f} | "
        f"timeout={mean_breakdown.timeout:.2f}"
    )


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
