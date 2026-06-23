# pyright:standard
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def move_block_scipy():
    """
    方块移动最优控制问题 - Python版本
    使用scipy.optimize.minimize求解器
    """

    # 系统参数
    T = 1.0  # 总时间
    N = 20  # 分段数
    dt = T / N  # 时间步长

    # 时间网格
    t = np.linspace(0, T, N + 1)

    # 优化变量维度
    # x = [x_0, x_1, ..., x_N, v_0, v_1, ..., v_N, u_0, u_1, ..., u_N]
    n_states = 2  # 状态变量数 (位置, 速度)
    n_controls = 1  # 控制变量数 (力)
    n_vars = (N + 1) * n_states + (N + 1) * n_controls

    # 初始猜测值
    x0 = np.zeros(n_vars)

    # 为状态变量设置合理的初始猜测
    for k in range(N + 1):
        # 位置：从0到1的线性插值
        x0[k] = k / N

        # 速度：简单的正弦波
        # x0[N+1+k] = 0.5 * np.sin(2*np.pi*k/N)
        x0[N + 1 + k] = 0
        # 控制力：简单的正弦波
        x0[2 * (N + 1) + k] = 0  # 2 * np.sin(2*np.pi*k/N)

    # 边界约束
    bounds = []

    # 位置和速度约束（无限制）
    for i in range(2 * (N + 1)):
        bounds.append((None, None))

    # 控制力约束
    u_max = 10.0
    for i in range(N + 1):
        bounds.append((-u_max, u_max))

    # 约束函数
    constraints = []

    # 动力学约束
    def dynamics_constraint(x):
        """动力学约束：梯形积分"""
        x_pos = x[: (N + 1)]
        x_vel = x[(N + 1) : (2 * (N + 1))]
        u_force = x[(2 * (N + 1)) :]

        constraints = []
        for k in range(N):
            # 位置约束：x_{k+1} = x_k + (dt/2) * (v_{k+1} + v_k)
            constraint_pos = x_pos[k + 1] - x_pos[k] - dt / 2 * (x_vel[k + 1] + x_vel[k])
            constraints.append(constraint_pos)

            # 速度约束：v_{k+1} = v_k + (dt/2) * (u_{k+1} + u_k)
            constraint_vel = x_vel[k + 1] - x_vel[k] - dt / 2 * (u_force[k + 1] + u_force[k])
            constraints.append(constraint_vel)

        return np.array(constraints)

    # 边界条件约束
    def boundary_constraint(x):
        """边界条件约束"""
        x_pos = x[: (N + 1)]
        x_vel = x[(N + 1) : (2 * (N + 1))]

        # 初始条件：x(0) = 0, v(0) = 0
        # 终端条件：x(T) = 1, v(T) = 0
        constraints = [x_pos[0], x_vel[0], x_pos[-1] - 1, x_vel[-1]]  # x(0) = 0  # v(0) = 0  # x(T) = 1  # v(T) = 0

        return np.array(constraints)

    # 添加约束
    constraints.append({"type": "eq", "fun": dynamics_constraint})

    constraints.append({"type": "eq", "fun": boundary_constraint})

    # 目标函数
    def objective_function(x):
        """目标函数：控制力平方和的积分"""
        u_start = 2 * (N + 1)
        u_force = x[u_start:]

        J = 0
        for k in range(N):
            J += dt / 2 * (u_force[k] ** 2 + u_force[k + 1] ** 2)

        return J

    # 优化选项
    options = {"maxiter": 1000, "disp": True, "ftol": 1e-6, "xtol": 1e-6}

    # 运行优化
    print("开始优化...")
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
    x_pos = x_opt[: (N + 1)]
    x_vel = x_opt[(N + 1) : (2 * (N + 1))]
    u_force = x_opt[(2 * (N + 1)) :]

    # 绘制结果
    plot_results(t, x_pos, x_vel, u_force)

    # 保存结果
    results = {
        "t": t,
        "x_pos": x_pos,
        "x_vel": x_vel,
        "u_force": u_force,
        "objective_value": result.fun,
        "success": result.success,
        "iterations": result.nit,
    }

    with open("move_block_results_scipy.pkl", "wb") as f:
        pickle.dump(results, f)

    print("结果已保存到 move_block_results_scipy.pkl")

    return results


def plot_results(t, x_pos, x_vel, u_force):
    """
    绘制优化结果
    """
    plt.rcParams["font.sans-serif"] = ["SimHei"]  # 用来正常显示中文标签
    plt.rcParams["axes.unicode_minus"] = False  # 用来正常显示负号

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

    # 位置图
    ax1.plot(t, x_pos, "b-", linewidth=2, label="位置轨迹")
    ax1.plot(t, x_pos, "ko", markersize=6, markerfacecolor="k", label="节点点")
    ax1.set_xlabel("时间 (s)")
    ax1.set_ylabel("位置 (m)")
    ax1.set_title("\n方块位置")
    ax1.grid(True)
    ax1.set_ylim([0, 1.2])
    ax1.legend()

    # 速度图
    ax2.plot(t, x_vel, "r-", linewidth=2, label="速度轨迹")
    ax2.plot(t, x_vel, "ko", markersize=6, markerfacecolor="k", label="节点点")
    ax2.set_xlabel("时间 (s)")
    ax2.set_ylabel("速度 (m/s)")
    ax2.set_title("方块速度")
    ax2.grid(True)
    ax2.legend()

    # 控制力图
    ax3.plot(t, u_force, "g-", linewidth=2, label="控制力")
    ax3.plot(t, u_force, "ko", markersize=6, markerfacecolor="k", label="节点点")
    ax3.set_xlabel("时间 (s)")
    ax3.set_ylabel("力 (N)")
    ax3.set_title("控制力")
    ax3.grid(True)
    ax3.legend()

    # 调整子图间距
    plt.tight_layout()
    plt.suptitle("控制结果 (Python scipy)", fontsize=14, fontweight="bold", y=0.98)

    # 保存图片
    plt.savefig("move_block_results_scipy.png", dpi=300, bbox_inches="tight")
    print("结果图片已保存为 move_block_results_scipy.png")

    plt.show()


def run_optimization():
    """
    运行脚本
    """
    print("=== 方块移动最优控制问题 (Python scipy版本) ===")
    print("系统参数:")
    print("  总时间: 1秒")
    print("  分段数: 20")
    print("  初始位置: 0米")
    print("  目标位置: 1米")
    print("  控制力约束: ±10N")
    print("  求解器: SLSQP")
    print()

    # 运行优化
    results = move_block_scipy()

    print("\n=== 优化完成 ===")
    print("请查看生成的图片和pkl文件以获取详细结果")

    return results


if __name__ == "__main__":
    results = run_optimization()
