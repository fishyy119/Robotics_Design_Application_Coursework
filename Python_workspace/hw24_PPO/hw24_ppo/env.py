# pyright: standard

from dataclasses import dataclass
from typing import cast

import gymnasium as gym
import numpy as np
from gymnasium.spaces import Box, Discrete

from .types import ObsArray


@dataclass(frozen=True)
class EnvSpec:
    obs_dim: int
    act_dim: int


@dataclass(frozen=True)
class RewardBreakdown:
    position: float = 0.0
    velocity: float = 0.0
    angle: float = 0.0
    leg_contact: float = 0.0
    main_engine: float = 0.0
    side_engine: float = 0.0
    curriculum_progress: float = 0.0
    terminal: float = 0.0
    timeout: float = 0.0


class LunarLanderRewardTracker:
    def __init__(self) -> None:
        self._previous_terms: RewardBreakdown | None = None

    def reset(self) -> None:
        self._previous_terms = None

    def step(
        self,
        action: int,
        observation: ObsArray,
        reward: float,
        terminated: bool,
    ) -> RewardBreakdown:
        if terminated:
            return RewardBreakdown(terminal=reward)

        current_terms = extract_lunar_lander_shaping_terms(observation)
        main_engine = -0.30 if action == 2 else 0.0
        side_engine = -0.03 if action in (1, 3) else 0.0

        if self._previous_terms is None:
            breakdown = RewardBreakdown(
                main_engine=main_engine,
                side_engine=side_engine,
            )
        else:
            breakdown = RewardBreakdown(
                position=current_terms.position - self._previous_terms.position,
                velocity=current_terms.velocity - self._previous_terms.velocity,
                angle=current_terms.angle - self._previous_terms.angle,
                leg_contact=current_terms.leg_contact - self._previous_terms.leg_contact,
                main_engine=main_engine,
                side_engine=side_engine,
            )

        self._previous_terms = current_terms
        return breakdown


def make_environment(
    env_id: str,
    render_mode: str | None = None,
) -> gym.Env[object, object]:
    env = cast(gym.Env[object, object], gym.make(env_id, render_mode=render_mode))
    validate_environment(env)
    return env


def validate_environment(env: gym.Env[object, object]) -> EnvSpec:
    observation_space = env.observation_space
    action_space = env.action_space

    if not isinstance(action_space, Discrete):
        raise TypeError("This PPO implementation only supports discrete action spaces.")
    if not isinstance(observation_space, Box):
        raise TypeError("This PPO implementation only supports Box observation spaces.")
    if len(observation_space.shape) != 1:
        raise TypeError("This PPO implementation only supports one-dimensional vector observations.")

    obs_dim = int(observation_space.shape[0])
    act_dim = int(action_space.n)
    return EnvSpec(obs_dim=obs_dim, act_dim=act_dim)


def to_observation_array(observation: object, obs_dim: int) -> ObsArray:
    array = np.asarray(observation, dtype=np.float32)
    if array.shape != (obs_dim,):
        raise ValueError(f"Expected observation shape ({obs_dim},), got {array.shape}.")
    return cast(ObsArray, array)


def add_reward_breakdowns(left: RewardBreakdown, right: RewardBreakdown) -> RewardBreakdown:
    return RewardBreakdown(
        position=left.position + right.position,
        velocity=left.velocity + right.velocity,
        angle=left.angle + right.angle,
        leg_contact=left.leg_contact + right.leg_contact,
        main_engine=left.main_engine + right.main_engine,
        side_engine=left.side_engine + right.side_engine,
        curriculum_progress=left.curriculum_progress + right.curriculum_progress,
        terminal=left.terminal + right.terminal,
        timeout=left.timeout + right.timeout,
    )


def mean_reward_breakdown(
    breakdowns: list[RewardBreakdown],
) -> RewardBreakdown | None:
    if not breakdowns:
        return None

    count = float(len(breakdowns))
    totals = RewardBreakdown()
    for breakdown in breakdowns:
        totals = add_reward_breakdowns(totals, breakdown)

    return RewardBreakdown(
        position=totals.position / count,
        velocity=totals.velocity / count,
        angle=totals.angle / count,
        leg_contact=totals.leg_contact / count,
        main_engine=totals.main_engine / count,
        side_engine=totals.side_engine / count,
        curriculum_progress=totals.curriculum_progress / count,
        terminal=totals.terminal / count,
        timeout=totals.timeout / count,
    )


def make_reward_tracker(env_id: str) -> LunarLanderRewardTracker | None:
    if env_id == "LunarLander-v3":
        return LunarLanderRewardTracker()
    return None


def extract_lunar_lander_shaping_terms(observation: ObsArray) -> RewardBreakdown:
    x_position = float(observation[0])
    y_position = float(observation[1])
    x_velocity = float(observation[2])
    y_velocity = float(observation[3])
    angle = float(observation[4])
    left_leg_contact = float(observation[6])
    right_leg_contact = float(observation[7])

    return RewardBreakdown(
        position=-100.0 * float(np.sqrt(x_position**2 + y_position**2)),
        velocity=-100.0 * float(np.sqrt(x_velocity**2 + y_velocity**2)),
        angle=-100.0 * abs(angle),
        leg_contact=10.0 * left_leg_contact + 10.0 * right_leg_contact,
    )
