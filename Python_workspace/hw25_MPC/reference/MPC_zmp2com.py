# pyright:standard
# #   download from:  https://github.com/machines-in-motion/lmpc_walking/tree/master/second_order
#
#    LMPC_walking is a python software implementation of some of the linear MPC
#    algorithms based presented in:
#    https://groups.csail.mit.edu/robotics-center/public_papers/Wieber15.pdf
#    Copyright (C) 2019 @ahmad gazar

#    LMPC_walking is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    LMPC_walking is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


# time vs CoP and CoM in x: 'A.K.A run rabbit run !'
# -------------------------------------------------
def plot_x(time, walking_time, min_admissible_CoP, max_admissible_cop, Z_x, X, Z_ref):
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")
    ZMP_x_fig = plt.figure()
    plt.plot(time, Z_x, label=r"\textbf{computed CoP}")
    plt.plot(time, Z_ref[:, 0])
    plt.plot(time, X[:, 0], "g", label=r"\textbf{CoM}")
    plt.plot(time, min_admissible_CoP[0:walking_time, 0], "r", linestyle="--", linewidth=0.5, label=r"\textbf{min CoP}")
    plt.plot(time, max_admissible_cop[0:walking_time, 0], "m", linestyle="--", linewidth=0.5, label=r"\textbf{max CoP}")

    plt.xlabel(r"\textbf{time} (s)")
    plt.ylabel(r"\textbf{x, z} (m)")
    plt.legend()
    # plt.title("time from"  ,str(time[0]), str(time[15]))


# time VS CoP and CoM in y: 'A.K.A what goes up must go down'
# ----------------------------------------------------------
def plot_y(time, walking_time, min_admissible_CoP, max_admissible_cop, Z_y, Y, Z_ref):
    ZMP_y_fig = plt.figure()
    plt.plot(time, Z_y, label=r"\textbf{computed CoP}")
    plt.plot(time, Z_ref[:, 1])
    plt.plot(time, Y[0:walking_time, 0], "g", label=r"\textbf{CoM}")
    plt.plot(time, min_admissible_CoP[0:walking_time, 1], "r", linestyle="--", linewidth=0.5, label=r"\textbf{min CoP}")
    plt.plot(time, max_admissible_cop[0:walking_time, 1], "m", linestyle="--", linewidth=0.5, label=r"\textbf{max CoP}")

    plt.xlabel(r"\textbf{time} (s)")
    plt.ylabel(r"\textbf{y, z} (m)")
    plt.legend()


# plot CoP, CoM in x Vs Cop, CoM in y:
# ------------------------------------
def plot_xy(time, walking_time, foot_length, foot_width, Z_ref, Z_x, Z_y, X, Y):
    ZMP_CoP_xy_fig = plt.figure()
    plt.plot(Z_x, Z_y, "r", label=r"\textbf{computed CoP}")
    plt.plot(X[0:walking_time, 0], Y[0:walking_time, 0], "lime", label=r"\textbf{CoM}")
    currentAxis = plt.gca()
    for i in range(walking_time):
        current_foot = patches.Rectangle(
            (Z_ref[i, 0] - foot_length / 2, Z_ref[i, 1] - foot_width / 2),
            foot_length,
            foot_width,
            linewidth=0.8,
            linestyle="-.",
            edgecolor="b",
            facecolor="none",
        )
        currentAxis.add_patch(current_foot)
    currentAxis.set_xlim((-0.5, 5.0))
    currentAxis.set_ylim((-0.5, 0.8))
    plt.xlabel(r"\textbf{x} (m)")
    plt.ylabel(r"\textbf{y} (m)")
    plt.legend()
    plt.show()


# plot every single horizon for debugging:
# ---------------------------------------
def plot_horizons(desired_walking_time, N, desired_Z_ref, horizon_data, foot_length, foot_width):
    for i in range(desired_walking_time):
        time_k = horizon_data[i]["time_k"]
        Z_ref_k = horizon_data[i]["zmp_reference"]
        X_k = horizon_data[i]["X_k"]
        Y_k = horizon_data[i]["Y_k"]
        Z_x_k = horizon_data[i]["Z_x_k"]
        Z_y_k = horizon_data[i]["Z_y_k"]

        min_admissible_CoP = Z_ref_k - np.tile([foot_length / 2, foot_width / 2], (N, 1))
        max_admissible_cop = Z_ref_k + np.tile([foot_length / 2, foot_width / 2], (N, 1))

        # plot_x(time_k, N, min_admissible_CoP, max_admissible_cop, \
        #                  Z_x_k, X_k, Z_ref_k)
        plot_y(time_k, N, min_admissible_CoP, max_admissible_cop, Z_y_k, Y_k, Z_ref_k)
        # plot_xy(time_k, N, foot_length, foot_width, desired_Z_ref, \
        #                   Z_x_k, Z_y_k, X_k, Y_k)


import math

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import matrix_power
from quadprog import solve_qp
from scipy.linalg import toeplitz

# cost weights in the objective function:
# ---------------------------------------
alpha = 10 ** (-1)  # CoP error squared cost weight
beta = 0 * 10 ** (-3)  # CoM position error squared cost weight
gamma = 10 ** (-3)  # CoM velocity error squared cost weight

# Inverted pendulum parameters:
# ----------------------------
h = 0.80  # fixed CoM height (assuming walking on a flat terrain)
g = 9.81  # norm of the gravity vector
foot_length = 0.20  # foot size in the x-direction
foot_width = 0.10  # foot size in the y-direciton

# MPC Parameters:
# --------------
delta_t = 0.1  # sampling time interval
step_time = 0.8  # time needed for every step
no_steps_per_T = int(round(step_time / delta_t))

# walking parameters:
# ------------------
step_length = 0.21  # fixed step length in the xz-plane
no_desired_steps = 6  # number of desired walking steps
desired_walking_time = no_desired_steps * no_steps_per_T  # number of desired walking intervals
N = desired_walking_time  # preceding horizon

# CoM initial state: [x_0, xdot_0].T
#                    [y_0, ydot_0].T
# ----------------------------------
x_0 = np.array([0.0, 0.0])
y_0 = np.array([-0.09, 0.0])

step_width = 2 * np.absolute(y_0[0])


# Desctiption:
# -----------
# this function assembles A and b encapsulating the ZMP inequality constraints
# in the Quadratic Program of this form (quadprog solver format):
#        minimize_u
#            (1/2) * u.T * Q * u - p.T * u
#        subject to
#            A.T * u >= b
# The foot shape is assumed to be rectangular

# Parameters:
# ----------
# N           : preceding horizon length                (scalar)
# foot_length : length of the foot along the x-axis     (scalar)
# foot_width  : length of the foot along the y-axis     (scalar)
# Z_ref_k     : [z_ref_x_k  , z_ref_y_k  ]   CoP reference trajectory
#                   .       ,    .                      (Nx2 numpy.array)
#                   .       ,    .
#               [z_ref_x_k+N, z_ref_y_k+N]
# x_hat_k     : [x_k, x_dot_k].T current CoM state in x (2, numpy.array)
# y_hat_k     : [y_k, y_dot_k].T current CoM state in y (2, numpy.array)

# Returns:
# -------
# A  : (4Nx2N numpy.array)
#      matrix defining defining the linear terms in the CoP control inputs in
#      x and y directions taking this form:
#      Z_x_k <= Zx_ref_k + foot_length/2
#      Z_x_k >= Zx_ref_k - foot_length/2
#      Z_y_k <= Zy_ref_k + foot_width/2
#      Z_y_k >= Zy_ref_k - foot_width/2
# b  : (4N, numpy.array)
#      vector defining the remaining terms


def add_ZMP_constraints(N, foot_length, foot_width, Z_ref_k, x_hat_k, y_hat_k):

    # pre-allocate memory
    A = np.zeros((4 * N, 2 * N))
    b = np.zeros((4 * N))
    foot_length_N = np.zeros((N))
    foot_width_N = np.zeros((N))

    # x-direction
    A[0:N, 0:N] = np.eye(N)
    A[N : 2 * N, 0:N] = -np.eye(N)

    # y-directions
    A[2 * N : 3 * N, N : 2 * N] = np.eye(N)
    A[3 * N : 4 * N, N : 2 * N] = -np.eye(N)

    foot_length_N = np.tile(foot_length, (N))
    foot_width_N = np.tile(foot_width, (N))

    b[0:N] = Z_ref_k[:, 0] - (0.5 * foot_length_N)
    b[N : 2 * N] = -Z_ref_k[:, 0] - (0.5 * foot_length_N)
    b[2 * N : 3 * N] = Z_ref_k[:, 1] - (0.5 * foot_width_N)
    b[3 * N : 4 * N] = -Z_ref_k[:, 1] - (0.5 * foot_width_N)

    return A, b


# Desctiption:
# -----------
# this function assembles A_eq matrix and b_eq vector encapsulating the CoM
# equality terminal constraints at the end of the preceding horizon

# Parameters:
# ----------
# N              : preceding horizon length                 (scalar)
# terminal_index : array index of the terminal constraint   (scalar)
# x_hat_k        : [x^_k, x^dot_k].T current CoM state in x (2, numpy.array)
# y_hat_k        : [y^_k, y^dot_k].T current CoM state in y (2, numpy.array)
# x_terminal     : [x_t, x_dot_t].T terminal CoM state in x (2, numpy.array)
# y_terminal     : [y_t, y_dot_t].T terminal CoM state in y (2, numpy.array)
# P_ps, P_vs     : CoM position recursive dynamics matrices (Nx2 numpy.array,
#                  x^_k+1 = P_ps x^_k + P_pu z_k            (NXN numpy.array)
# P_pu , P_vu    : CoM velocity recursive dynamics matrix   (Nx2 numpy.array,
#                  x^dot_k+1 = P_vs x^dot_k + P_vu z_k      (NXN numpy.array)

# Returns:
# -------
# A_eq  : (4x2N numpy.array)
#         matrix defining the linear terms in the CoM state. However, since the
#         MPC problem decision variables are only the CoP control inputs
#         then the CoM equality terminal constraints can be formulated in terms
#         of the CoP control inputs in x and y as follows:
#         P_ps x^_k + P_pu[terminal_index,:] z_x_k = x_terminal[0]
#         P_vs x^_k + P_vu[terminal_index,:] z_x_k = x_terminal[0]

# b_eq  : (4, numpy.array)
#         vector defining the remaining terms


def add_terminal_constraints(N, terminal_index, x_hat_k, y_hat_k, x_terminal, y_terminal, P_ps, P_vs, P_pu, P_vu):

    # pre-allocate memory
    A_eq = np.zeros((4, 2 * N))
    b_eq = np.zeros((4))

    P_ps_current = P_ps[terminal_index, :]
    P_vs_current = P_vs[terminal_index, :]

    P_pu_current = P_pu[terminal_index, :]
    P_vu_current = P_vu[terminal_index, :]

    # x-direction
    A_eq[0, 0:N] = P_pu_current
    A_eq[1, 0:N] = P_vu_current

    # y-direction
    A_eq[2, N : 2 * N] = P_pu_current
    A_eq[3, N : 2 * N] = P_vu_current

    b_eq[0] = x_terminal[0] - np.dot(P_ps_current, x_hat_k)
    b_eq[1] = x_terminal[1] - np.dot(P_vs_current, x_hat_k)
    b_eq[2] = y_terminal[0] - np.dot(P_ps_current, y_hat_k)
    b_eq[3] = y_terminal[1] - np.dot(P_vs_current, y_hat_k)

    return A_eq, b_eq


# Description:
# -----------
# this function compute the canonical quadratic objective term Q (hessian)
# and the linear objective term p.T (gradient) of this cost function:
# min zx_k, zy_k
#    beta/2||x-x_r||^2 + gamma/2||x_dot-x_dot_r||^2 + alpha ||z_x-zx_r||^2
#  + beta/2||y-y_r||^2 + gamma/2||y_dot-y_dot_r||^2 + alpha ||z_y-zy_r||^2

# Parameters:
# ----------
# alpha          : CoP error squared cost weight            (scalar)
# beta           : CoM position error squared cost weight   (scalar)
# gamma          : CoM velocity error squared cost weight   (scalar)
# N              : preceding horizon                        (scalar)
# P_ps, P_vs     : CoM position recursive dynamics matrices (Nx2 numpy.array,
#                  x^_k+1 = P_ps x^_k + P_pu z_k             NXN numpy.array)
# P_pu , P_vu    : CoM velocity recursive dynamics matrix   (Nx2 numpy.array,
#                  x^dot_k+1 = P_vs x^dot_k + P_vu z_k       NXN numpy.array)
# x_hat_k        : [x^_k, x^dot_k].T current CoM state in x (2, numpy.array)
# y_hat_k        : [y^_k, y^dot_k].T current CoM state in y (2, numpy.array)
# Z_ref_k        : [z_ref_x_k  , z_ref_y_k  ]   CoP reference trajectory
#                   .       ,    .                          (Nx2 numpy.array)
#                   .       ,    .
#                  [z_ref_x_k+N, z_ref_y_k+N]

# Returns:
# -------
# Q     : Hessian  (2Nx2N numpy.array)
# p_k   : Gradient (2N,  numpy.array)


def compute_objective_terms(
    alpha,
    beta,
    gamma,
    step_duration,
    no_steps_per_T,
    N,
    stride_length,
    stride_width,
    P_ps,
    P_pu,
    P_vs,
    P_vu,
    x_hat_k,
    y_hat_k,
    Z_ref_k,
):

    # pre-allocate memory
    Q = np.zeros((2 * N, 2 * N))
    p_k = np.zeros((2 * N))
    Q_prime = np.zeros((N, N))

    Q_prime = alpha * np.eye(N) + (beta * np.dot(P_pu.T, P_pu)) + (gamma * np.dot(P_vu.T, P_vu))
    Q[0:N, 0:N] = Q_prime  # x-direction
    Q[N : 2 * N, N : 2 * N] = Q_prime  # y-direction

    x_r_N = np.zeros((N))
    y_r_N = np.zeros((N))
    x_dotr_N = np.zeros((N))
    y_dotr_N = np.zeros((N))

    x_r_N = np.tile(stride_length / no_steps_per_T, N)
    y_r_N = np.tile(stride_width / no_steps_per_T, N)
    x_dotr_N = np.tile(stride_length / step_duration, N)
    y_dotr_N = np.tile(stride_width / step_duration, N)
    # print('HHHHHHHHHHHH',Z_ref_k.shape)
    # print("P_ps shape:", P_ps.shape)
    # print("P_vs shape:", P_vs.shape)
    # print("P_pu shape:", P_pu.shape)
    # print("P_vu shape:", P_vu.shape)
    # print(f"p_k[0:N] shape: {p_k[0:N].shape}")
    # print(f"gamma term shape: {(gamma * (np.dot(P_vu.T, (np.dot(P_vs, x_hat_k)))- np.dot(P_vu.T, x_dotr_N))).shape}")
    # print(f"beta term shape: {(beta  * (np.dot(P_pu.T, (np.dot(P_ps, x_hat_k)))- np.dot(P_pu.T, x_r_N))).shape}")
    # print(f"alpha term shape: {(alpha*Z_ref_k[:,0]).shape}")

    p_k[0:N] = (
        gamma * (np.dot(P_vu.T, (np.dot(P_vs, x_hat_k))) - np.dot(P_vu.T, x_dotr_N))
        + beta * (np.dot(P_pu.T, (np.dot(P_ps, x_hat_k))) - np.dot(P_pu.T, x_r_N))
        - alpha * Z_ref_k[:, 0]
    )  # x-direction

    p_k[N : 2 * N] = (
        gamma * (np.dot(P_vu.T, (np.dot(P_vs, y_hat_k))) - np.dot(P_vu.T, y_dotr_N))
        + beta * (np.dot(P_pu.T, (np.dot(P_ps, y_hat_k))) - np.dot(P_pu.T, y_r_N))
        - alpha * Z_ref_k[:, 1]
    )  # y-direction
    return Q, p_k


# Generate trajectory using 3rd order polynomial with following constraints:
# x(0)=x0, x(T)=x1, dx(0)=dx(T)=0
# x(t) = a + b t + c t^2 + d t^3
# x(0) = a = x0
# dx(0) = b = 0
# dx(T) = 2 c T + 3 d T^2 = 0 => c = -3 d T^2 / (2 T) = -(3/2) d T
# x(T) = x0 + c T^2 + d T^3 = x1
#        x0 -(3/2) d T^3 + d T^3 = x1
#        -0.5 d T^3 = x1 - x0
#        d = 2 (x0-x1) / T^3
# c = -(3/2) T 2 (x0-x1) / (T^3) = 3 (x1-x0) / T^2
def compute_3rd_order_poly_traj(x0, x1, T, dt):
    a = x0
    b = np.zeros_like(x0)
    c = 3 * (x1 - x0) / (T**2)
    d = 2 * (x0 - x1) / (T**3)
    N = int(T / dt)
    n = x0.shape[0]
    x = np.zeros((n, N))
    dx = np.zeros((n, N))
    ddx = np.zeros((n, N))
    for i in range(N):
        t = i * dt
        x[:, i] = a + b * t + c * t**2 + d * t**3
        dx[:, i] = b + 2 * c * t + 3 * d * t**2
        ddx[:, i] = 2 * c + 6 * d * t
    return x, dx, ddx


def compute_foot_traj(foot_steps, N, dt, step_time, step_height, first_phase):
    x = np.zeros((3, N + 1))
    dx = np.zeros((3, N + 1))
    ddx = np.zeros((3, N + 1))
    N_step = int(step_time / dt)
    offset = 0
    if first_phase == "swing":
        offset = N_step
        x[0, :N_step] = foot_steps[0, 0]
        x[1, :N_step] = foot_steps[0, 1]

    for s in range(foot_steps.shape[0]):
        i = offset + s * 2 * N_step
        x[0, i : i + N_step] = foot_steps[s, 0]
        x[1, i : i + N_step] = foot_steps[s, 1]
        if s < foot_steps.shape[0] - 1:
            next_step = foot_steps[s + 1, :]
        elif first_phase == "swing":
            break
        else:
            next_step = foot_steps[s, :]
            step_height = 0.0
        (
            x[:2, i + N_step : i + 2 * N_step],
            dx[:2, i + N_step : i + 2 * N_step],
            ddx[:2, i + N_step : i + 2 * N_step],
        ) = compute_3rd_order_poly_traj(foot_steps[s, :], next_step, step_time, dt)

        (
            x[2, i + N_step : i + int(1.5 * N_step)],
            dx[2, i + N_step : i + int(1.5 * N_step)],
            ddx[2, i + N_step : i + int(1.5 * N_step)],
        ) = compute_3rd_order_poly_traj(np.array([0.0]), np.array([step_height]), 0.5 * step_time, dt)

        (
            x[2, i + int(1.5 * N_step) : i + 2 * N_step],
            dx[2, i + int(1.5 * N_step) : i + 2 * N_step],
            ddx[2, i + int(1.5 * N_step) : i + 2 * N_step],
        ) = compute_3rd_order_poly_traj(np.array([step_height]), np.array([0.0]), 0.5 * step_time, dt)

    return x, dx, ddx


def interpolate_lipm_traj(T_step, nb_steps, dt_mpc, dt_ctrl, com_z, g, com_state_x, com_state_y, cop_ref, cop_x, cop_y):
    # INTERPOLATE WITH TIME STEP OF CONTROLLER (TSID)
    N = nb_steps * int(round(T_step / dt_mpc))  # number of time steps for traj-opt
    N_ctrl = int((N * dt_mpc) / dt_ctrl)  # number of time steps for TSID
    com = np.empty((3, N_ctrl + 1)) * np.nan
    dcom = np.zeros((3, N_ctrl + 1))
    ddcom = np.zeros((3, N_ctrl + 1))
    cop = np.empty((2, N_ctrl + 1)) * np.nan
    foot_steps = np.empty((2, N_ctrl + 1)) * np.nan
    contact_phase = (N_ctrl + 1) * ["right"]
    com[2, :] = com_z

    N_inner = int(N_ctrl / N)
    for i in range(N):
        com[0, i * N_inner] = com_state_x[i, 0]
        com[1, i * N_inner] = com_state_y[i, 0]
        dcom[0, i * N_inner] = com_state_x[i, 1]
        dcom[1, i * N_inner] = com_state_y[i, 1]
        if i > 0:
            if np.linalg.norm(cop_ref[i, :] - cop_ref[i - 1, :]) < 1e-10:
                contact_phase[i * N_inner] = contact_phase[i * N_inner - 1]
            else:
                if contact_phase[(i - 1) * N_inner] == "right":
                    contact_phase[i * N_inner] = "left"
                elif contact_phase[(i - 1) * N_inner] == "left":
                    contact_phase[i * N_inner] = "right"

        for j in range(N_inner):
            ii = i * N_inner + j
            A, B = discrete_LIP_dynamics((j + 1) * dt_ctrl, g, com_z)
            foot_steps[:, ii] = cop_ref[i, :].T
            cop[0, ii] = cop_x[i]
            cop[1, ii] = cop_y[i]
            x_next = A.dot(com_state_x[i, :]) + B.dot(cop[0, ii])
            y_next = A.dot(com_state_y[i, :]) + B.dot(cop[1, ii])
            com[0, ii + 1] = x_next[0]
            com[1, ii + 1] = y_next[0]
            dcom[0, ii + 1] = x_next[1]
            dcom[1, ii + 1] = y_next[1]
            ddcom[:2, ii] = g / com_z * (com[:2, ii] - cop[:, ii])

            if j > 0:
                contact_phase[ii] = contact_phase[ii - 1]
    return com, dcom, ddcom, cop, contact_phase, foot_steps


# Description:
# -----------
# this function returns the discrete dynamics matrices A_d and B_d
# of the linear inverted pendulum x+ = A_d x + B_d x

# Parameters:
# ----------
# delta_t: sampling time
# g      : norm of the gravity acceleration vector
# h      : fixed height of the CoM assuming walking on a flat terrain

# Returns:
# -------
# A_d (2x2 numpy.array)
# B_d (2,  numpy.array)


def discrete_LIP_dynamics(delta_t, g, h):
    w = math.sqrt(g / h)
    A_d = np.array(
        [
            [math.cosh(w * delta_t), (1 / w) * math.sinh(w * delta_t)],
            [w * math.sinh(w * delta_t), math.cosh(w * delta_t)],
        ]
    )

    B_d = np.array([1 - math.cosh(w * delta_t), -w * math.sinh(w * delta_t)])

    return A_d, B_d


# Description:
# -----------
# this function computes the integration of the discrete dynamic matrices of the
# linear inverted pendulum. the matrices are constructed once offline since
# all parameters used in the computation are fixed.
#  x^_k+1    = P_ps x^_k    + P_pu z_k
#  x^dot_k+1 = P_vs x^dot_k + P_vu z_k

# Parameters:
# ----------
# delta_t: sampling time
# g      : norm of the gravity acceleration vector
# h      : fixed height of the CoM assuming walking on a flat terrain

# Returns:
# -------
# P_ps, P_vs  : position and velocity partitions of the states recursive
#                dynamics matrices (Nx2 numpy.array, NXN numpy.array)
# P_pu , P_vu : position and velocity partitions of the control inputs recursive
#                dynamics matrix   (Nx2 numpy.array, NXN numpy.array)


def compute_recursive_matrices(delta_t, g, h, N):
    [A_d, B_d] = discrete_LIP_dynamics(delta_t, g, h)

    # pre-allocate memmory
    P_ps = np.zeros((N, 2))
    P_vs = np.zeros((N, 2))
    temp_pu = np.zeros((N))
    temp_vu = np.zeros((N))

    for i in range(N):
        A_d_pow = matrix_power(A_d, i + 1)  # numpy.linalg.matrix_power 是一个用于计算矩阵的幂的函数
        P_ps[i, 0:2] = A_d_pow[0, :]  # A_d_pow[0,:]取第 0 行的所有列
        # P_ps:Nx2 numpy.array
        P_vs[i, 0:2] = A_d_pow[1, :]  # A_d_pow[1,:]取第 1 行的所有列
        temp_u = np.dot(matrix_power(A_d, i), B_d)
        temp_pu[i] = temp_u[0]
        temp_vu[i] = temp_u[1]
    # with open("P_ps_pv.txt", "w", encoding="utf-8") as f:
    #    print('loop number = ', i, '\n', file=f)
    #    print('P_ps = ', P_ps, '\n', file=f)
    #    print('P_vs = ', P_vs, '\n', file=f)
    """
        print('loop number = ', i, '\n')
        print('A_d_pow  = ', A_d_pow, '\n')
        print('P_ps = ', P_ps, '\n')
        print('temp_pu[i] = ', temp_pu[i], '\n')
        print('temp_vu[i] = ', temp_vu[i], '\n')
    """
    P_pu = toeplitz(temp_pu) * np.tri(N, N)
    P_vu = toeplitz(temp_vu) * np.tri(N, N)
    """
    with open("P_pu_vu.txt", "w", encoding="utf-8") as f:
        print('P_pu = ', P_pu, '\n',file=f)
        print('P_vu = ', P_vu, '\n',file=f)
    """
    print("P_pu = ", P_pu, "\n")
    print("P_vu = ", P_vu, "\n")
    return P_ps, P_vs, P_pu, P_vu


# Description:
# -----------
# this function computes the integration of the recursive discrete dynamics
# of the linear inverted pendulum based on the current initial CoM state.

# Parameters:
# ----------
# P_ps, P_vs  : position and velocity partitions of the states recursive
#                dynamics matrices          (Nx2 numpy.array, NXN numpy.array)
# P_pu , P_vu : position and velocity partitions of the control inputs recursive
#                dynamics matrix            (Nx2 numpy.array, NXN numpy.array)
# N           : preceding horizon           (scalar)
# x_hat_k     : [x^_k, x^dot_k].T current CoM state in x (2, numpy.array)
# y_hat_k     : [y^_k, y^dot_k].T current CoM state in y (2, numpy.array)
# U_k         : [zx_k, ... , zx_k+N, ... ,zy_k, ..., zy_k+N].T
#               current CoP control inputs in x, y directions
#               along the horizon (2N, numpy.array)

# Returns:
# -------
# X: [x_k+1   , x_dot_k+1  ]   recursive CoM dynamics in x direction
#       .     ,   .            (Nx2 numpy.array)
#       .     ,   .
#    [x_k+1+N , x_dot_k+1+N]

# Y: [y_k+1   , y_dot_k+1  ]   recursive CoM dynamics in y direction
#       .     ,   .            (Nx2 numpy.array)
#       .     ,   .
#    [y_k+1+N , y_dot_k+1+N]


def compute_recursive_dynamics(P_ps, P_vs, P_pu, P_vu, N, x_hat_k, y_hat_k, U_k):
    # pre-allocate memory
    X = np.zeros((N, 2))
    Y = np.zeros((N, 2))

    # evaluate your CoM states in the x-direction along the horizon
    X[0:N, 0] = np.dot(P_ps, x_hat_k) + np.dot(P_pu, U_k[0:N])  # x
    X[0:N, 1] = np.dot(P_vs, x_hat_k) + np.dot(P_vu, U_k[0:N])  # x_dot

    # evaluate your CoM states in the y-direction along the horizon
    Y[0:N, 0] = np.dot(P_ps, y_hat_k) + np.dot(P_pu, U_k[N : 2 * N])  # y
    Y[0:N, 1] = np.dot(P_vs, y_hat_k) + np.dot(P_vu, U_k[N : 2 * N])  # y_dot

    return X, Y


# Description:
# -----------
# this function implements a desired zig-zag fixed foot step plan located
# in the middle of the robot's foot starting with the right foot

# Parameters:
# ----------
#  foot_step_0 : initial foot step location
#            [foot_step_x0, foot_step_y0].T (2x1 numpy.array)
#  no_steps    : number of desired walking foot steps  (scalar)

# Returns:
# -------
#  Foot_steps = [Foot_steps_x, Foot_steps_y].T   foot steps locations
#                                                (no_steps x 2 numpy.array)


def manual_foot_placement(foot_step_0, fixed_step_x, no_steps):
    Foot_steps = np.zeros((no_steps, 2))
    for i in range(Foot_steps.shape[0]):
        if i == 0:
            Foot_steps[i, :] = foot_step_0
        else:
            Foot_steps[i, 0] = Foot_steps[i - 1, 0] + fixed_step_x
            Foot_steps[i, 1] = -Foot_steps[i - 1, 1]
    return Foot_steps


# Description:
# -----------
# this function computes a CoP reference trajectory based on a desired
# fixed foot step plan, a desired foot step duration and a sampling time

# Parameters:
# ----------
#  no_steps      : number of desired walking foot steps  (scalar)
#  Foot_steps    := [Foot_steps_x, Foot_steps_y]   foot steps locations
#                                                  (no_steps x 2 numpy.array)
# walking_time   : desired walking time duration   (scalar)
# no_steps_per_T : step_duration/T  (scalar)

# Returns:
# -------
# Z_ref  := [z_ref_x_k             , z_ref_y_k              ]   CoP reference trajectory
#                   .              ,    .                       (walking_timex2 numpy.array)
#                   .              ,    .
#           [z_ref_x_k+walking_time, z_ref_y_k+walking_time]


def create_CoP_trajectory(no_steps, Foot_steps, walking_time, no_steps_per_T):
    Z_ref = np.zeros((walking_time, 2))
    j = 0
    for i in range(Foot_steps.shape[0]):
        Z_ref[j : j + no_steps_per_T, :] = Foot_steps[i, :]
        j = j + no_steps_per_T
    return Z_ref


# compute CoP reference trajectory:
# --------------------------------
foot_step_0 = np.array([0.0, -0.09])  # initial foot step position in x-y

desiredFoot_steps = manual_foot_placement(foot_step_0, step_length, no_desired_steps)
desired_Z_ref = create_CoP_trajectory(no_desired_steps, desiredFoot_steps, desired_walking_time, no_steps_per_T)

# used in case you want to have terminal constraints
# -------------------------------------------------
x_terminal = np.array([desired_Z_ref[N - 1, 0], 0.0])  # CoM terminal constraint in x : [x, xdot].T
y_terminal = np.array([desired_Z_ref[N - 1, 1], 0.0])  # CoM terminal constraint in y : [y, ydot].T
no_terminal_constraints = 4
terminal_index = N - 1

# construct your preview system: 'Go pokemon !'
# --------------------------------------------
[P_ps, P_vs, P_pu, P_vu] = compute_recursive_matrices(delta_t, g, h, N)
[Q, p_k] = compute_objective_terms(
    alpha,
    beta,
    gamma,
    step_time,
    no_steps_per_T,
    N,
    step_length,
    step_width,
    P_ps,
    P_pu,
    P_vs,
    P_vu,
    x_0,
    y_0,
    desired_Z_ref,
)
[A_zmp, b_zmp] = add_ZMP_constraints(N, foot_length, foot_width, desired_Z_ref, x_0, y_0)

# used in case you want to add both terminal add_ZMP_constraints
# --------------------------------------------------------------
[A_terminal, b_terminal] = add_terminal_constraints(
    N, terminal_index, x_0, y_0, x_terminal, y_terminal, P_ps, P_vs, P_pu, P_vu
)
A = np.concatenate((A_terminal, A_zmp), axis=0)
b = np.concatenate((b_terminal, b_zmp), axis=0)

# call quadprog solver:
# --------------------
# U = solve_qp(Q, -p_k, A.T, b, no_terminal_constraints)[0]  # uncomment to solve with 4 equality terminal constraints
U = solve_qp(Q, -p_k, A_zmp.T, b_zmp)[0]  # solve only with only CoP inequality constraints
Z_x_total = U[0:N]
Z_y_total = U[N : 2 * N]

# Trajectory optimization: (based on the initial state x_hat_0, y_hat_0)
# -------------------------------------------------------------------------
[X_total, Y_total] = compute_recursive_dynamics(P_ps, P_vs, P_pu, P_vu, N, x_0, y_0, U)
# ------------------------------------------------------------------------------
# visualize your open-loop trajectory:
# ------------------------------------------------------------------------------
time = np.arange(0, round(desired_walking_time * delta_t, 2), delta_t)
min_admissible_CoP = desired_Z_ref - np.tile([foot_length / 2, foot_width / 2], (desired_walking_time, 1))
max_admissible_cop = desired_Z_ref + np.tile([foot_length / 2, foot_width / 2], (desired_walking_time, 1))

# time vs CoP and CoM in x: 'A.K.A run rabbit run !'
# -------------------------------------------------
plot_x(time, desired_walking_time, min_admissible_CoP, max_admissible_cop, Z_x_total, X_total, desired_Z_ref)


# time VS CoP and CoM in y: 'A.K.A what goes up must go down'
# ----------------------------------------------------------
plot_y(time, desired_walking_time, min_admissible_CoP, max_admissible_cop, Z_y_total, Y_total, desired_Z_ref)

# plot CoP, CoM in x Vs Cop, CoM in y:
# -----------------------------------
plot_xy(time, desired_walking_time, foot_length, foot_width, desired_Z_ref, Z_x_total, Z_y_total, X_total, Y_total)
