# pyright:standard
import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def cartpole_optimize():
    """
    倒立摆-小车最优控制问题
    使用scipy.optimize.minimize优化器求解
    """

    # 系统参数
    T = 2.0  # 总时间
    N = 30  # 分段数
    m1 = 1.0  # 小车质量
    m2 = 0.3  # 摆锤质量
    g = 9.81  # 重力加速度
    l = 0.5  # 摆杆长度

    # 约束参数
    d = 1.0  # 目标位置
    d_max = 2.0  # 位置约束
    u_max = 20.0  # 控制力约束

    # 时间网格
    t = np.linspace(0, T, N + 1)
    h = np.diff(t)  # 时间步长

    # 优化变量维度
    # x = [q1_0, q1_1, ..., q1_N, q2_0, q2_1, ..., q2_N,
    #      dq1_0, dq1_1, ..., dq1_N, dq2_0, dq2_1, ..., dq2_N,
    #      u_0, u_1, ..., u_N]
    n_states = 4  # 状态变量数
    n_controls = 1  # 控制变量数
    n_vars = (N + 1) * n_states + (N + 1) * n_controls

    # 初始猜测值
    x0 = np.zeros(n_vars)

    # 为状态变量设置合理的初始猜测
    for k in range(N + 1):
        # 位置：从0到d的线性插值
        x0[k] = d * k / N

        # 角度：从0到π的平滑过渡
        x0[N + 1 + k] = np.pi * (1 - np.cos(np.pi * k / N)) / 2

        # 速度：初始为0
        x0[2 * (N + 1) + k] = 0
        x0[3 * (N + 1) + k] = 0

        # 控制力：简单的正弦波
        x0[4 * (N + 1) + k] = 5 * np.sin(2 * np.pi * k / N)

    # 边界约束
    bounds = []

    # 位置约束
    for i in range(N + 1):
        bounds.append((-d_max, d_max))

    # 角度约束（无限制）
    for i in range(N + 1):
        bounds.append((None, None))

    # 速度约束（无限制）
    for i in range(2 * (N + 1)):
        bounds.append((None, None))

    # 控制力约束
    for i in range(N + 1):
        bounds.append((-u_max, u_max))

    # 约束函数
    constraints = []

    # 动力学约束
    def dynamics_constraint(x):
        """动力学约束：x_{k+1} = x_k + (h_k/2) * (f_{k+1} + f_k)"""
        q1 = x[: (N + 1)]
        q2 = x[(N + 1) : (2 * (N + 1))]
        dq1 = x[(2 * (N + 1)) : (3 * (N + 1))]
        dq2 = x[(3 * (N + 1)) : (4 * (N + 1))]
        u = x[(4 * (N + 1)) :]

        constraints = []
        for k in range(N):
            # 当前时刻状态
            x_k = np.array([q1[k], q2[k], dq1[k], dq2[k]])
            u_k = u[k]

            # 下一时刻状态
            x_kp1 = np.array([q1[k + 1], q2[k + 1], dq1[k + 1], dq2[k + 1]])
            u_kp1 = u[k + 1]

            # 计算动力学
            f_k = cartpole_dynamics(x_k, u_k, m1, m2, g, l)
            f_kp1 = cartpole_dynamics(x_kp1, u_kp1, m1, m2, g, l)

            # 梯形积分约束
            constraint = x_kp1 - x_k - h[k] / 2 * (f_k + f_kp1)
            constraints.extend(constraint)

        return np.array(constraints)

    # 边界条件约束
    def boundary_constraint(x):
        """边界条件约束"""
        q1 = x[: (N + 1)]
        q2 = x[(N + 1) : (2 * (N + 1))]
        dq1 = x[(2 * (N + 1)) : (3 * (N + 1))]
        dq2 = x[(3 * (N + 1)) : (4 * (N + 1))]

        # 初始条件：x(0) = [0, 0, 0, 0]
        # 终端条件：x(T) = [d, π, 0, 0]
        constraints = [
            q1[0],  # q1(0) = 0
            q2[0],  # q2(0) = 0
            dq1[0],  # dq1(0) = 0
            dq2[0],  # dq2(0) = 0
            q1[-1] - d,  # q1(T) = d
            q2[-1] - np.pi,  # q2(T) = π
            dq1[-1],  # dq1(T) = 0
            dq2[-1],  # dq2(T) = 0
        ]

        return np.array(constraints)

    # 添加约束
    constraints.append({"type": "eq", "fun": dynamics_constraint})

    constraints.append({"type": "eq", "fun": boundary_constraint})

    # 目标函数
    def objective_function(x):
        """目标函数：控制力平方和的积分"""
        u_start = 4 * (N + 1)
        u = x[u_start:]

        J = 0
        for k in range(N):
            J += h[k] / 2 * (u[k] ** 2 + u[k + 1] ** 2)

        return J

    # 优化选项
    options = {"maxiter": 1000, "disp": True, "ftol": 1e-6, "xtol": 1e-6}

    # 运行优化
    result = minimize(objective_function, x0, method="SLSQP", bounds=bounds, constraints=constraints, options=options)

    # 检查优化结果
    if result.success:
        print("优化成功完成！")
        print(f"目标函数值: {result.fun:.6f}")
        print(f"迭代次数: {result.nit}")
    else:
        print("优化未成功完成")
        print(f"退出原因: {result.message}")

    # 提取优化结果
    x_opt = result.x
    q1_opt = x_opt[: (N + 1)]
    q2_opt = x_opt[(N + 1) : (2 * (N + 1))]
    dq1_opt = x_opt[(2 * (N + 1)) : (3 * (N + 1))]
    dq2_opt = x_opt[(3 * (N + 1)) : (4 * (N + 1))]
    u_opt = x_opt[(4 * (N + 1)) :]

    # 绘制结果
    plot_results(t, q1_opt, q2_opt, u_opt)

    # 保存结果
    results = {
        "t": t,
        "q1_opt": q1_opt,
        "q2_opt": q2_opt,
        "dq1_opt": dq1_opt,
        "dq2_opt": dq2_opt,
        "u_opt": u_opt,
        "objective_value": result.fun,
        "success": result.success,
        "iterations": result.nit,
    }

    with open("cartpole_results.pkl", "wb") as f:
        pickle.dump(results, f)

    print("结果已保存到 cartpole_results.pkl")

    return results


def cartpole_dynamics(x, u, m1, m2, g, l):
    """
    倒立摆动力学方程
    """
    q1 = x[0]  # 小车位置
    q2 = x[1]  # 摆杆角度
    dq1 = x[2]  # 小车速度
    dq2 = x[3]  # 摆杆角速度

    # 计算加速度
    den = m1 + m2 * (1 - np.cos(q2) ** 2)

    # 小车加速度
    ddq1 = (l * m2 * np.sin(q2) * dq2**2 + u + m2 * g * np.cos(q2) * np.sin(q2)) / den

    # 摆杆角加速度
    ddq2 = -(l * m2 * np.cos(q2) * np.sin(q2) * dq2**2 + u * np.cos(q2) + (m1 + m2) * g * np.sin(q2)) / (l * den)

    return np.array([dq1, dq2, ddq1, ddq2])


def plot_results(t, q1, q2, u):
    """
    绘制优化结果
    """
    plt.rcParams["font.sans-serif"] = ["SimHei"]  # 用来正常显示中文标签
    plt.rcParams["axes.unicode_minus"] = False  # 用来正常显示负号

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

    # 位置图
    ax1.plot(t, q1, "b-", linewidth=2, label="cubic spline (state)")
    ax1.plot(t, q1, "ko", markersize=6, markerfacecolor="k", label="knot-points")
    ax1.set_xlabel("时间 (s)")
    ax1.set_ylabel("位置 (m)")
    ax1.set_title("\n小车位置")
    ax1.grid(True)
    ax1.set_ylim([0, 1.5])
    ax1.legend()

    # 角度图
    ax2.plot(t, q2, "b-", linewidth=2, label="cubic spline (state)")
    ax2.plot(t, q2, "ko", markersize=6, markerfacecolor="k", label="knot-points")
    ax2.set_xlabel("时间 (s)")
    ax2.set_ylabel("角度 (rad)")
    ax2.set_title("摆杆角度")
    ax2.grid(True)
    ax2.set_ylim([-2, 4])
    ax2.legend()

    # 控制力图
    ax3.plot(t, u, "m-", linewidth=2, label="quadratic spline (control)")
    ax3.plot(t, u, "ko", markersize=6, markerfacecolor="k", label="knot-points")
    ax3.set_xlabel("时间 (s)")
    ax3.set_ylabel("力 (N)")
    ax3.set_title("控制力")
    ax3.grid(True)
    ax3.set_ylim([-20, 10])
    ax3.legend()

    # 调整子图间距
    plt.tight_layout()
    plt.suptitle("倒立摆-小车最优控制结果", fontsize=14, fontweight="bold", y=0.98)

    # 保存图片
    plt.savefig("cartpole_optimization_results.png", dpi=300, bbox_inches="tight")
    print("结果图片已保存为 cartpole_optimization_results.png")

    plt.show()


def run_optimization():
    """
    运行脚本
    """
    print("=== 倒立摆-小车最优控制问题 ===")
    print("系统参数:")
    print("  总时间: 2秒")
    print("  分段数: 30")
    print("  目标位置: 1米")
    print("  摆杆目标角度: π弧度")
    print("  控制力约束: ±20N")
    print("  位置约束: ±2米")
    print()

    # 运行优化
    results = cartpole_optimize()

    print("\n=== 优化完成 ===")
    print("请查看生成的图片和pkl文件以获取详细结果")

    return results


if __name__ == "__main__":
    results = run_optimization()
