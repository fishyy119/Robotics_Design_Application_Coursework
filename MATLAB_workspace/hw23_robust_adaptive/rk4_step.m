function x_next = rk4_step(x, tau, Ts, params)
% 用 RK4 在一个采样周期内积分机械臂动力学。

k1 = manipulator_rhs(x, tau, params);
k2 = manipulator_rhs(x + 0.5 * Ts * k1, tau, params);
k3 = manipulator_rhs(x + 0.5 * Ts * k2, tau, params);
k4 = manipulator_rhs(x + Ts * k3, tau, params);
x_next = x + Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
