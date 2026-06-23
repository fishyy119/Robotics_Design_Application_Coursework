# pyright:standard
"""
hw25 正式代码：基于 LIPM 的开环 / 闭环 MPC CoM 轨迹对比.

该脚本保留参考代码中的线性倒立摆建模与 CoP 约束形式，
并补充滚动时域的闭环 MPC，用于比较状态检测误差下的鲁棒性差异。
"""

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import matrix_power
from numpy.typing import NDArray
from quadprog import solve_qp
from scipy.linalg import toeplitz

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


@dataclass(frozen=True)
class MPCConfig:
    alpha: float = 1.0e-1
    beta: float = 1.0e-2
    gamma: float = 1.0e-3
    h: float = 0.80
    g: float = 9.81
    foot_length: float = 0.20
    foot_width: float = 0.10
    delta_t: float = 0.1
    step_time: float = 0.8
    step_length: float = 0.21
    no_desired_steps: int = 6
    x0: tuple[float, float] = (0.0, 0.0)
    y0: tuple[float, float] = (-0.09, 0.0)
    desired_forward_velocity: float = 0.2625
    seed: int = 42
    initial_state_offset: tuple[float, float, float, float] = (1.0e-4, 2.0e-4, -1.0e-4, -2.0e-4)
    process_noise_std: tuple[float, float, float, float] = (1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8)
    measurement_noise_std: tuple[float, float, float, float] = (1.0e-8, 1.0e-8, 1.0e-8, 1.0e-8)

    @property
    def no_steps_per_T(self) -> int:
        return int(round(self.step_time / self.delta_t))

    @property
    def walking_time(self) -> int:
        return self.no_desired_steps * self.no_steps_per_T

    @property
    def nominal_x0(self) -> NDArray:
        return np.array(self.x0, dtype=float)

    @property
    def nominal_y0(self) -> NDArray:
        return np.array(self.y0, dtype=float)

    @property
    def disturbed_x0(self) -> NDArray:
        return self.nominal_x0 + np.array(self.initial_state_offset[:2], dtype=float)

    @property
    def disturbed_y0(self) -> NDArray:
        return self.nominal_y0 + np.array(self.initial_state_offset[2:], dtype=float)


@dataclass
class SimulationResult:
    cop: NDArray
    x_state: NDArray
    y_state: NDArray
    measured_x_state: NDArray | None = None
    measured_y_state: NDArray | None = None


def manual_foot_placement(foot_step_0: NDArray, fixed_step_x: float, no_steps: int) -> NDArray:
    foot_steps = np.zeros((no_steps, 2))
    for i in range(no_steps):
        if i == 0:
            foot_steps[i, :] = foot_step_0
        else:
            foot_steps[i, 0] = foot_steps[i - 1, 0] + fixed_step_x
            foot_steps[i, 1] = -foot_steps[i - 1, 1]
    return foot_steps


def create_cop_trajectory(foot_steps: NDArray, no_steps_per_T: int) -> NDArray:
    walking_time = foot_steps.shape[0] * no_steps_per_T
    z_ref = np.zeros((walking_time, 2))
    index = 0
    for foot in foot_steps:
        z_ref[index : index + no_steps_per_T, :] = foot
        index += no_steps_per_T
    return z_ref


def discrete_lip_dynamics(delta_t: float, g: float, h: float) -> tuple[NDArray, NDArray]:
    omega = math.sqrt(g / h)
    a_d = np.array(
        [
            [math.cosh(omega * delta_t), math.sinh(omega * delta_t) / omega],
            [omega * math.sinh(omega * delta_t), math.cosh(omega * delta_t)],
        ]
    )
    b_d = np.array([1.0 - math.cosh(omega * delta_t), -omega * math.sinh(omega * delta_t)])
    return a_d, b_d


def compute_recursive_matrices(
    delta_t: float, g: float, h: float, horizon: int
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    a_d, b_d = discrete_lip_dynamics(delta_t, g, h)
    p_ps = np.zeros((horizon, 2))
    p_vs = np.zeros((horizon, 2))
    temp_pu = np.zeros(horizon)
    temp_vu = np.zeros(horizon)

    for i in range(horizon):
        a_d_pow = matrix_power(a_d, i + 1)
        p_ps[i, :] = a_d_pow[0, :]
        p_vs[i, :] = a_d_pow[1, :]

        temp_u = matrix_power(a_d, i) @ b_d
        temp_pu[i] = temp_u[0]
        temp_vu[i] = temp_u[1]

    p_pu = toeplitz(temp_pu) * np.tri(horizon, horizon)
    p_vu = toeplitz(temp_vu) * np.tri(horizon, horizon)
    return p_ps, p_vs, p_pu, p_vu


def add_zmp_constraints(horizon: int, foot_length: float, foot_width: float, z_ref: NDArray) -> tuple[NDArray, NDArray]:
    a = np.zeros((4 * horizon, 2 * horizon))
    b = np.zeros(4 * horizon)

    a[0:horizon, 0:horizon] = np.eye(horizon)
    a[horizon : 2 * horizon, 0:horizon] = -np.eye(horizon)
    a[2 * horizon : 3 * horizon, horizon : 2 * horizon] = np.eye(horizon)
    a[3 * horizon : 4 * horizon, horizon : 2 * horizon] = -np.eye(horizon)

    half_length = 0.5 * foot_length
    half_width = 0.5 * foot_width
    b[0:horizon] = z_ref[:, 0] - half_length
    b[horizon : 2 * horizon] = -z_ref[:, 0] - half_length
    b[2 * horizon : 3 * horizon] = z_ref[:, 1] - half_width
    b[3 * horizon : 4 * horizon] = -z_ref[:, 1] - half_width
    return a, b


def compute_objective_terms(
    config: MPCConfig,
    p_ps: NDArray,
    p_pu: NDArray,
    p_vs: NDArray,
    p_vu: NDArray,
    x_hat: NDArray,
    y_hat: NDArray,
    z_ref: NDArray,
) -> tuple[NDArray, NDArray]:
    horizon = z_ref.shape[0]
    q = np.zeros((2 * horizon, 2 * horizon))
    p = np.zeros(2 * horizon)

    q_block = config.alpha * np.eye(horizon) + config.beta * (p_pu.T @ p_pu) + config.gamma * (p_vu.T @ p_vu)
    q[0:horizon, 0:horizon] = q_block
    q[horizon : 2 * horizon, horizon : 2 * horizon] = q_block

    x_ref = z_ref[:, 0]
    y_ref = z_ref[:, 1]
    x_dot_ref = np.full(horizon, config.desired_forward_velocity)
    y_dot_ref = np.zeros(horizon)

    p[0:horizon] = (
        config.gamma * (p_vu.T @ (p_vs @ x_hat) - p_vu.T @ x_dot_ref)
        + config.beta * (p_pu.T @ (p_ps @ x_hat) - p_pu.T @ x_ref)
        - config.alpha * z_ref[:, 0]
    )
    p[horizon : 2 * horizon] = (
        config.gamma * (p_vu.T @ (p_vs @ y_hat) - p_vu.T @ y_dot_ref)
        + config.beta * (p_pu.T @ (p_ps @ y_hat) - p_pu.T @ y_ref)
        - config.alpha * z_ref[:, 1]
    )
    return q, p


def solve_preview_mpc(
    config: MPCConfig,
    x_hat: NDArray,
    y_hat: NDArray,
    z_ref: NDArray,
    matrix_cache: dict[int, tuple[NDArray, NDArray, NDArray, NDArray]],
) -> NDArray:
    horizon = z_ref.shape[0]
    if horizon not in matrix_cache:
        matrix_cache[horizon] = compute_recursive_matrices(config.delta_t, config.g, config.h, horizon)

    p_ps, p_vs, p_pu, p_vu = matrix_cache[horizon]
    q, p = compute_objective_terms(config, p_ps, p_pu, p_vs, p_vu, x_hat, y_hat, z_ref)
    a_zmp, b_zmp = add_zmp_constraints(horizon, config.foot_length, config.foot_width, z_ref)
    return solve_qp(q, -p, a_zmp.T, b_zmp)[0]


def rollout_lipm(
    a_d: NDArray,
    b_d: NDArray,
    cop: NDArray,
    x0: NDArray,
    y0: NDArray,
    disturbances: NDArray | None = None,
) -> tuple[NDArray, NDArray]:
    horizon = cop.shape[0]
    x_state = np.zeros((horizon + 1, 2))
    y_state = np.zeros((horizon + 1, 2))
    x_state[0, :] = x0
    y_state[0, :] = y0

    if disturbances is None:
        disturbances = np.zeros((horizon, 4))

    for k in range(horizon):
        x_state[k + 1, :] = a_d @ x_state[k, :] + b_d * cop[k, 0] + disturbances[k, 0:2]
        y_state[k + 1, :] = a_d @ y_state[k, :] + b_d * cop[k, 1] + disturbances[k, 2:4]
    return x_state, y_state


def simulate_open_loop(
    config: MPCConfig,
    z_ref: NDArray,
    disturbances: NDArray,
    initial_measurement_error: NDArray,
    matrix_cache: dict[int, tuple[NDArray, NDArray, NDArray, NDArray]],
) -> tuple[SimulationResult, SimulationResult]:
    a_d, b_d = discrete_lip_dynamics(config.delta_t, config.g, config.h)
    nominal_u = solve_preview_mpc(config, config.nominal_x0, config.nominal_y0, z_ref, matrix_cache)
    actual_x0 = config.disturbed_x0
    actual_y0 = config.disturbed_y0
    measured_x0 = actual_x0 + initial_measurement_error[0:2]
    measured_y0 = actual_y0 + initial_measurement_error[2:4]
    disturbed_u = solve_preview_mpc(config, measured_x0, measured_y0, z_ref, matrix_cache)
    horizon = z_ref.shape[0]
    nominal_cop = np.column_stack((nominal_u[0:horizon], nominal_u[horizon : 2 * horizon]))
    disturbed_cop = np.column_stack((disturbed_u[0:horizon], disturbed_u[horizon : 2 * horizon]))

    nominal_x_state, nominal_y_state = rollout_lipm(a_d, b_d, nominal_cop, config.nominal_x0, config.nominal_y0)
    disturbed_x_state, disturbed_y_state = rollout_lipm(
        a_d,
        b_d,
        disturbed_cop,
        actual_x0,
        actual_y0,
        disturbances=disturbances,
    )
    nominal_result = SimulationResult(cop=nominal_cop, x_state=nominal_x_state, y_state=nominal_y_state)
    disturbed_result = SimulationResult(cop=disturbed_cop, x_state=disturbed_x_state, y_state=disturbed_y_state)
    return nominal_result, disturbed_result


def simulate_closed_loop(
    config: MPCConfig,
    z_ref: NDArray,
    disturbances: NDArray,
    measurement_noise: NDArray,
    matrix_cache: dict[int, tuple[NDArray, NDArray, NDArray, NDArray]],
) -> SimulationResult:
    horizon_total = z_ref.shape[0]
    a_d, b_d = discrete_lip_dynamics(config.delta_t, config.g, config.h)
    cop = np.zeros((horizon_total, 2))
    x_state = np.zeros((horizon_total + 1, 2))
    y_state = np.zeros((horizon_total + 1, 2))
    measured_x_state = np.zeros((horizon_total + 1, 2))
    measured_y_state = np.zeros((horizon_total + 1, 2))

    x_state[0, :] = config.disturbed_x0
    y_state[0, :] = config.disturbed_y0

    for k in range(horizon_total):
        measured_x_state[k, :] = x_state[k, :] + measurement_noise[k, 0:2]
        measured_y_state[k, :] = y_state[k, :] + measurement_noise[k, 2:4]

        u_k = solve_preview_mpc(config, measured_x_state[k, :], measured_y_state[k, :], z_ref[k:, :], matrix_cache)
        horizon_k = horizon_total - k
        cop[k, 0] = u_k[0]
        cop[k, 1] = u_k[horizon_k]

        x_state[k + 1, :] = a_d @ x_state[k, :] + b_d * cop[k, 0] + disturbances[k, 0:2]
        y_state[k + 1, :] = a_d @ y_state[k, :] + b_d * cop[k, 1] + disturbances[k, 2:4]

    measured_x_state[horizon_total, :] = x_state[horizon_total, :] + measurement_noise[horizon_total, 0:2]
    measured_y_state[horizon_total, :] = y_state[horizon_total, :] + measurement_noise[horizon_total, 2:4]
    return SimulationResult(
        cop=cop,
        x_state=x_state,
        y_state=y_state,
        measured_x_state=measured_x_state,
        measured_y_state=measured_y_state,
    )


def compute_position_error(reference: SimulationResult, candidate: SimulationResult) -> NDArray:
    reference_pos = np.column_stack((reference.x_state[:, 0], reference.y_state[:, 0]))
    candidate_pos = np.column_stack((candidate.x_state[:, 0], candidate.y_state[:, 0]))
    return np.linalg.norm(candidate_pos - reference_pos, axis=1)


def first_threshold_crossing(error: NDArray, delta_t: float, threshold: float) -> float:
    above = np.nonzero(error > threshold)[0]
    if above.size == 0:
        return -1
    return float(above[0] * delta_t)


def summarize_results(
    reference: SimulationResult, open_loop: SimulationResult, closed_loop: SimulationResult, delta_t: float
) -> dict[str, float]:
    open_error = compute_position_error(reference, open_loop)
    closed_error = compute_position_error(reference, closed_loop)
    return {
        "open_loop_rms_error_m": float(np.sqrt(np.mean(open_error**2))),
        "closed_loop_rms_error_m": float(np.sqrt(np.mean(closed_error**2))),
        "open_loop_max_error_m": float(np.max(open_error)),
        "closed_loop_max_error_m": float(np.max(closed_error)),
        "open_loop_final_error_m": float(open_error[-1]),
        "closed_loop_final_error_m": float(closed_error[-1]),
        "open_loop_first_cross_0p10m_s": first_threshold_crossing(open_error, delta_t, 0.10),
        "closed_loop_first_cross_0p10m_s": first_threshold_crossing(closed_error, delta_t, 0.10),
    }


def plot_results(
    config: MPCConfig,
    foot_steps: NDArray,
    z_ref: NDArray,
    nominal_result: SimulationResult,
    open_loop_result: SimulationResult,
    closed_loop_result: SimulationResult,
    output_dir: Path,
) -> list[Path]:
    time_state = np.arange(config.walking_time + 1) * config.delta_t
    time_input = np.arange(config.walking_time) * config.delta_t

    nominal_error = np.full(config.walking_time + 1, 1.0e-12)
    open_error = compute_position_error(nominal_result, open_loop_result)
    closed_error = compute_position_error(nominal_result, closed_loop_result)
    combined_x = np.concatenate(
        (nominal_result.x_state[:, 0], open_loop_result.x_state[:, 0], closed_loop_result.x_state[:, 0], z_ref[:, 0])
    )
    combined_y = np.concatenate(
        (nominal_result.y_state[:, 0], open_loop_result.y_state[:, 0], closed_loop_result.y_state[:, 0], z_ref[:, 1])
    )
    x_margin = 0.08
    y_margin = 0.04
    x_min = float(np.min(combined_x) - x_margin)
    x_max = float(np.max(combined_x) + x_margin)
    y_min = float(np.min(combined_y) - y_margin)
    y_max = float(np.max(combined_y) + y_margin)
    nominal_com_alpha = 0.5
    disturbed_com_alpha = 0.7
    closed_com_alpha = 1.0

    saved_paths: list[Path] = []

    fig_x, ax_x = plt.subplots(figsize=(7.2, 4.8))
    ax_x.plot(
        time_state,
        nominal_result.x_state[:, 0],
        label="Open-loop CoM (no disturbance)",
        linewidth=2.0,
        alpha=nominal_com_alpha,
    )
    ax_x.plot(
        time_state,
        open_loop_result.x_state[:, 0],
        label="Open-loop CoM (with disturbance)",
        linestyle="--",
        alpha=disturbed_com_alpha,
    )
    ax_x.plot(
        time_state,
        closed_loop_result.x_state[:, 0],
        label="Closed-loop CoM (with disturbance)",
        linestyle="-.",
        alpha=closed_com_alpha,
    )
    ax_x.step(time_input, z_ref[:, 0], where="post", label="CoP reference")
    ax_x.set_title("CoM motion in x")
    ax_x.set_xlabel("Time (s)")
    ax_x.set_ylabel("x (m)")
    ax_x.grid(True, alpha=0.3)
    ax_x.legend()
    path_x = output_dir / "hw25_com_x.png"
    fig_x.tight_layout()
    fig_x.savefig(path_x, dpi=200, bbox_inches="tight")
    plt.close(fig_x)
    saved_paths.append(path_x)

    fig_y, ax_y = plt.subplots(figsize=(7.2, 4.8))
    ax_y.plot(
        time_state,
        nominal_result.y_state[:, 0],
        label="Open-loop CoM (no disturbance)",
        linewidth=2.0,
        alpha=nominal_com_alpha,
    )
    ax_y.plot(
        time_state,
        open_loop_result.y_state[:, 0],
        label="Open-loop CoM (with disturbance)",
        linestyle="--",
        alpha=disturbed_com_alpha,
    )
    ax_y.plot(
        time_state,
        closed_loop_result.y_state[:, 0],
        label="Closed-loop CoM (with disturbance)",
        linestyle="-.",
        alpha=closed_com_alpha,
    )
    ax_y.step(time_input, z_ref[:, 1], where="post", label="CoP reference")
    ax_y.set_title("CoM motion in y")
    ax_y.set_xlabel("Time (s)")
    ax_y.set_ylabel("y (m)")
    ax_y.grid(True, alpha=0.3)
    ax_y.legend()
    path_y = output_dir / "hw25_com_y.png"
    fig_y.tight_layout()
    fig_y.savefig(path_y, dpi=200, bbox_inches="tight")
    plt.close(fig_y)
    saved_paths.append(path_y)

    fig_xy, ax_xy = plt.subplots(figsize=(7.2, 5.6))
    ax_xy.plot(
        nominal_result.x_state[:, 0],
        nominal_result.y_state[:, 0],
        label="Open-loop CoM (no disturbance)",
        linewidth=2.0,
        alpha=nominal_com_alpha,
    )
    ax_xy.plot(
        open_loop_result.x_state[:, 0],
        open_loop_result.y_state[:, 0],
        label="Open-loop CoM (with disturbance)",
        linestyle="--",
        alpha=disturbed_com_alpha,
    )
    ax_xy.plot(
        closed_loop_result.x_state[:, 0],
        closed_loop_result.y_state[:, 0],
        label="Closed-loop CoM (with disturbance)",
        linestyle="-.",
        alpha=closed_com_alpha,
    )
    ax_xy.scatter(
        open_loop_result.cop[:, 0],
        open_loop_result.cop[:, 1],
        label="Open-loop CoP (with disturbance)",
        s=22,
        marker="o",
    )
    ax_xy.scatter(
        closed_loop_result.cop[:, 0],
        closed_loop_result.cop[:, 1],
        label="Closed-loop CoP (with disturbance)",
        s=26,
        marker="x",
    )
    for foot in foot_steps:
        foot_patch = patches.Rectangle(
            (foot[0] - config.foot_length / 2.0, foot[1] - config.foot_width / 2.0),
            config.foot_length,
            config.foot_width,
            linewidth=1.0,
            edgecolor="black",
            facecolor="none",
            linestyle=":",
        )
        ax_xy.add_patch(foot_patch)
    ax_xy.set_title("Trajectories in x-y plane")
    ax_xy.set_xlabel("x (m)")
    ax_xy.set_ylabel("y (m)")
    ax_xy.set_xlim(x_min, x_max)
    ax_xy.set_ylim(y_min, y_max)
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(loc="upper left")
    ax_xy.set_aspect("equal", adjustable="box")
    path_xy = output_dir / "hw25_xy_plane.png"
    fig_xy.tight_layout()
    fig_xy.savefig(path_xy, dpi=200, bbox_inches="tight")
    plt.close(fig_xy)
    saved_paths.append(path_xy)

    fig_err, ax_err = plt.subplots(figsize=(7.2, 4.8))
    ax_err.plot(time_state, nominal_error, label="Open-loop (no disturbance)", linewidth=2.0)
    ax_err.plot(time_state, open_error, label="Open-loop (with disturbance)", linestyle="--")
    ax_err.plot(time_state, closed_error, label="Closed-loop (with disturbance)", linestyle="-.")
    ax_err.set_title("Deviation from nominal CoM trajectory")
    ax_err.set_xlabel("Time (s)")
    ax_err.set_ylabel("Position error norm (m)")
    ax_err.set_yscale("log")
    ax_err.grid(True, alpha=0.3)
    ax_err.legend()
    path_err = output_dir / "hw25_error_log.png"
    fig_err.tight_layout()
    fig_err.savefig(path_err, dpi=200, bbox_inches="tight")
    plt.close(fig_err)
    saved_paths.append(path_err)

    return saved_paths


def save_summary(config: MPCConfig, summary: dict[str, float], output_path: Path) -> None:
    payload = {
        "config": asdict(config),
        "summary": summary,
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    output_dir = Path(__file__).resolve().parent / "artifacts"
    output_dir.mkdir(parents=True, exist_ok=True)

    config = MPCConfig()
    rng = np.random.default_rng(config.seed)

    foot_step_0 = np.array([0.0, config.y0[0]], dtype=float)
    foot_steps = manual_foot_placement(foot_step_0, config.step_length, config.no_desired_steps)
    z_ref = create_cop_trajectory(foot_steps, config.no_steps_per_T)

    disturbances = rng.normal(
        loc=0.0,
        scale=np.array(config.process_noise_std, dtype=float),
        size=(config.walking_time, 4),
    )
    measurement_noise = rng.normal(
        loc=0.0,
        scale=np.array(config.measurement_noise_std, dtype=float),
        size=(config.walking_time + 1, 4),
    )

    matrix_cache: dict[int, tuple[NDArray, NDArray, NDArray, NDArray]] = {}
    nominal_result, open_loop_result = simulate_open_loop(
        config, z_ref, disturbances, measurement_noise[0, :], matrix_cache
    )
    closed_loop_result = simulate_closed_loop(config, z_ref, disturbances, measurement_noise, matrix_cache)
    summary = summarize_results(nominal_result, open_loop_result, closed_loop_result, config.delta_t)

    summary_path = output_dir / "hw25_open_vs_closed_loop_mpc_summary.json"
    figure_paths = plot_results(
        config, foot_steps, z_ref, nominal_result, open_loop_result, closed_loop_result, output_dir
    )
    save_summary(config, summary, summary_path)

    for figure_path in figure_paths:
        print("Saved figure:", figure_path)
    print("Saved summary:", summary_path)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
