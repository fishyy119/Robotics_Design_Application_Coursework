from dataclasses import asdict, dataclass
from typing import Mapping, cast


@dataclass(frozen=True)
class ModelConfig:
    hidden_sizes: tuple[int, ...] = (256, 256)
    activation: str = "relu"


@dataclass(frozen=True)
class TrainConfig:
    env_id: str = "LunarLander-v3"
    seed: int = 42
    device: str = "cpu"
    steps_per_epoch: int = 8192
    epochs: int = 200
    gamma: float = 0.99
    lam: float = 0.95
    clip_ratio: float = 0.05
    pi_lr: float = 2e-6
    vf_lr: float = 2e-6
    train_pi_iters: int = 10
    train_v_iters: int = 10
    target_kl: float = 0.015
    max_ep_len: int = 512
    timeout_penalty: float = -30.0
    use_curriculum: bool = True
    curriculum_heights: tuple[float, ...] = (1.0, 0.5)
    curriculum_progress_rewards: tuple[float, ...] = (30.0, 60.0)
    curriculum_x_tolerance_scale: float = 0.5
    curriculum_angle_tolerance_scale: float = 0.5
    curriculum_horizontal_speed_tolerance: float = 0.05
    save_every: int = 10


@dataclass(frozen=True)
class EvalConfig:
    episodes: int = 10
    deterministic: bool = True
    render: bool = False


@dataclass(frozen=True)
class PlotConfig:
    pass


def serialize_model_config(config: ModelConfig) -> dict[str, object]:
    return {
        "hidden_sizes": list(config.hidden_sizes),
        "activation": config.activation,
    }


def serialize_train_config(config: TrainConfig) -> dict[str, object]:
    return asdict(config)


def build_config_payload(
    run_name: str,
    model_config: ModelConfig,
    train_config: TrainConfig,
) -> dict[str, object]:
    return {
        "run_name": run_name,
        "model_config": serialize_model_config(model_config),
        "train_config": serialize_train_config(train_config),
    }


def deserialize_model_config(data: Mapping[str, object]) -> ModelConfig:
    hidden_sizes_raw = data.get("hidden_sizes")
    if not isinstance(hidden_sizes_raw, list):
        raise ValueError("Checkpoint model_config.hidden_sizes is missing or invalid.")

    hidden_sizes_list: list[int] = []
    for value in cast(list[object], hidden_sizes_raw):
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError("Checkpoint model_config.hidden_sizes must contain integers.")
        hidden_sizes_list.append(value)

    activation = _require_str(data, "activation")
    hidden_sizes = tuple(hidden_sizes_list)
    return ModelConfig(hidden_sizes=hidden_sizes, activation=activation)


def deserialize_train_config(data: Mapping[str, object]) -> TrainConfig:
    return TrainConfig(
        env_id=_require_str(data, "env_id"),
        seed=_require_int(data, "seed"),
        device=_require_str(data, "device"),
        steps_per_epoch=_require_int(data, "steps_per_epoch"),
        epochs=_require_int(data, "epochs"),
        gamma=_require_float(data, "gamma"),
        lam=_require_float(data, "lam"),
        clip_ratio=_require_float(data, "clip_ratio"),
        pi_lr=_require_float(data, "pi_lr"),
        vf_lr=_require_float(data, "vf_lr"),
        train_pi_iters=_require_int(data, "train_pi_iters"),
        train_v_iters=_require_int(data, "train_v_iters"),
        target_kl=_require_float(data, "target_kl"),
        max_ep_len=_require_int(data, "max_ep_len"),
        timeout_penalty=_require_float(data, "timeout_penalty"),
        use_curriculum=_require_bool(data, "use_curriculum"),
        curriculum_heights=_require_float_tuple(data, "curriculum_heights"),
        curriculum_progress_rewards=_require_float_tuple(data, "curriculum_progress_rewards"),
        curriculum_x_tolerance_scale=_require_float(data, "curriculum_x_tolerance_scale"),
        curriculum_angle_tolerance_scale=_require_float(data, "curriculum_angle_tolerance_scale"),
        curriculum_horizontal_speed_tolerance=_require_float(data, "curriculum_horizontal_speed_tolerance"),
        save_every=_require_int(data, "save_every"),
    )


def _require_value(data: Mapping[str, object], key: str) -> object:
    if key not in data:
        raise ValueError(f"Checkpoint train_config is missing required field: {key}.")
    return data[key]


def _require_str(data: Mapping[str, object], key: str) -> str:
    value = _require_value(data, key)
    if not isinstance(value, str):
        raise ValueError(f"Checkpoint field {key} must be a string.")
    return value


def _require_int(data: Mapping[str, object], key: str) -> int:
    value = _require_value(data, key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"Checkpoint field {key} must be an integer.")
    return value


def _require_float(data: Mapping[str, object], key: str) -> float:
    value = _require_value(data, key)
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError(f"Checkpoint field {key} must be numeric.")
    return float(value)


def _require_bool(data: Mapping[str, object], key: str) -> bool:
    value = _require_value(data, key)
    if not isinstance(value, bool):
        raise ValueError(f"Checkpoint field {key} must be a boolean.")
    return value


def _require_float_tuple(data: Mapping[str, object], key: str) -> tuple[float, ...]:
    value = _require_value(data, key)
    if not isinstance(value, list | tuple):
        raise ValueError(f"Checkpoint field {key} must be a float sequence.")

    sequence = cast(list[object] | tuple[object, ...], value)
    numbers: list[float] = []
    for item in sequence:
        if isinstance(item, bool) or not isinstance(item, int | float):
            raise ValueError(f"Checkpoint field {key} must contain numeric values.")
        numbers.append(float(item))
    return tuple(numbers)
