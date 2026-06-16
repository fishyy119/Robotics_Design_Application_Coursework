from dataclasses import dataclass

from .config import TrainConfig
from .types import ObsArray


@dataclass(frozen=True)
class CurriculumStepResult:
    progress_reward: float = 0.0


class LunarLanderCurriculum:
    def __init__(self, config: TrainConfig) -> None:
        if not config.curriculum_heights:
            raise ValueError("curriculum_heights must not be empty.")
        if len(config.curriculum_progress_rewards) != len(config.curriculum_heights):
            raise ValueError("curriculum_progress_rewards must have the same length as curriculum_heights.")

        self.heights = config.curriculum_heights
        self.progress_rewards = config.curriculum_progress_rewards
        self.x_tolerance_scale = config.curriculum_x_tolerance_scale
        self.angle_tolerance_scale = config.curriculum_angle_tolerance_scale
        self.horizontal_speed_tolerance = config.curriculum_horizontal_speed_tolerance
        if self.x_tolerance_scale <= 0.0:
            raise ValueError("curriculum_x_tolerance_scale must be positive.")
        if self.angle_tolerance_scale <= 0.0:
            raise ValueError("curriculum_angle_tolerance_scale must be positive.")
        self.reached_milestones: list[bool] = [False] * len(self.heights)

    def reset_episode(self) -> None:
        self.reached_milestones = [False] * len(self.heights)

    def evaluate_step(
        self,
        observation: ObsArray,
        terminated: bool,
    ) -> CurriculumStepResult:
        if terminated:
            return CurriculumStepResult()

        reward_total = 0.0
        for index, target_height in enumerate(self.heights):
            if self.reached_milestones[index]:
                continue
            if not milestone_reached(
                observation=observation,
                target_height=target_height,
                x_tolerance_scale=self.x_tolerance_scale,
                angle_tolerance_scale=self.angle_tolerance_scale,
                horizontal_speed_tolerance=self.horizontal_speed_tolerance,
            ):
                continue

            self.reached_milestones[index] = True
            reward_total += self.progress_rewards[index]

        return CurriculumStepResult(progress_reward=reward_total)


def build_curriculum(config: TrainConfig) -> LunarLanderCurriculum | None:
    if not config.use_curriculum:
        return None
    if config.env_id != "LunarLander-v3":
        return None
    return LunarLanderCurriculum(config)


def milestone_reached(
    observation: ObsArray,
    target_height: float,
    x_tolerance_scale: float,
    angle_tolerance_scale: float,
    horizontal_speed_tolerance: float,
) -> bool:
    x_position = abs(float(observation[0]))
    y_position = float(observation[1])
    x_velocity = abs(float(observation[2]))
    angle = abs(float(observation[4]))
    x_tolerance = x_tolerance_scale * target_height
    angle_tolerance = angle_tolerance_scale * target_height

    return bool(
        y_position < target_height
        and x_position < x_tolerance
        and angle < angle_tolerance
        and x_velocity < horizontal_speed_tolerance
    )
