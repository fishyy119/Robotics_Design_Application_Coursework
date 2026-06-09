%% 8.14 关节空间自适应控制纯 MATLAB 仿真
% 基于 hw23 参考程序中的线性参数化形式，
% 对末端附加 10 kg 集中负载的平面机械臂进行离散时域仿真。

clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 参数
Ts = 1e-3;
t_move = 1.0;
t_sim = 4;
t = 0:Ts:t_sim;
q_i = [0; pi / 4];
q_f = [pi / 2; pi / 2];

params = build_model(10.0);
traj = build_joint_trajectory(struct( ...
    'q_i', q_i, ...
    'q_f', q_f, ...
    't_move', t_move, ...
    'Ts', Ts, ...
    't_total', t_sim, ...
    't_acc_ratio', 0.25));

Lam = 5 * eye(2);
K_d = 750 * eye(2);
K_ml = 0.01;

%% 初始条件
N = numel(t);
q = zeros(2, N);
dq = zeros(2, N);
ddq = zeros(2, N);
tau = zeros(2, N);
m_hat = zeros(1, N);

q(:, 1) = q_i;
m_hat(1) = 0;

%% 离散仿真
for k = 1:N - 1
    q_err = traj.q(:, k) - q(:, k);
    dq_err = traj.dq(:, k) - dq(:, k);

    dq_r = Lam * q_err + traj.dq(:, k);
    ddq_r = Lam * dq_err + traj.ddq(:, k);
    sig = dq_r - dq(:, k); % Lam * q_err + dq_err;

    [Y_r, Y_l] = build_regressor(q(:, k), dq(:, k), dq_r, ddq_r, params);
    tau(:, k) = Y_r * params.pi_m + params.F_v * dq_r + Y_l * m_hat(k) + K_d * sig;

    x_k = [q(:, k); dq(:, k)];
    x_kp1 = rk4_step(x_k, tau(:, k), Ts, params);
    q(:, k+1) = x_kp1(1:2);
    dq(:, k+1) = x_kp1(3:4);

    x_dot = manipulator_rhs(x_k, tau(:, k), params);
    ddq(:, k) = x_dot(3:4);

    dm_hat = (Y_l' * sig) / K_ml;
    % 对质量估计加入非负投影。
    m_hat(k+1) = max(0, m_hat(k)+Ts*dm_hat);
end

tau(:, N) = tau(:, N-1);
ddq(:, N) = ddq(:, N-1);

%% 后处理
[p_d_x, p_d_y] = forward_kinematics(traj.q, params.a);
[p_x, p_y] = forward_kinematics(q, params.a);
position_error_norm = sqrt((p_d_x - p_x).^2+(p_d_y - p_y).^2);

fprintf('8.14 adaptive control final end-effector error norm: %.6e m\n', position_error_norm(end));
fprintf('8.14 adaptive control final load estimate: %.6f kg\n', m_hat(end));

%% 绘图
fig1 = utils.createFigureA4(struct('Name', '8.14 Adaptive Control', 'Width', 18, 'AspectRatio', 0.70));
set(fig1, 'Color', 'w');

subplot(2, 2, 1);
plot(t, traj.q(1, :), '--', t, q(1, :), '-');
grid on;
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Joint 1 Angle');
legend({'Ref', 'Act'}, 'Location', 'best');

subplot(2, 2, 2);
plot(t, traj.q(2, :), '--', t, q(2, :), '-');
grid on;
xlabel('Time (s)');
ylabel('Angle (rad)');
title('Joint 2 Angle');
legend({'Ref', 'Act'}, 'Location', 'best');

subplot(2, 2, 3);
plot(t, tau(1, :), '-', t, tau(2, :), '--');
grid on;
xlabel('Time (s)');
ylabel('Torque (N m)');
title('Joint Torque');
legend({'Tau 1', 'Tau 2'}, 'Location', 'best');

subplot(2, 2, 4);
plot(t, position_error_norm, '-');
grid on;
xlabel('Time (s)');
ylabel('Error Norm (m)');
title('Norm of Position Error');

fig2 = utils.createFigureA4(struct('Name', '8.14 Load Estimate', 'Width', 14, 'AspectRatio', 0.72));
set(fig2, 'Color', 'w');
plot(t, m_hat, '-', 'LineWidth', 1.2);
hold on;
plot(t, 10*ones(size(t)), '--');
grid on;
xlabel('Time (s)');
ylabel('Load Estimate (kg)');
title('Estimated Load Mass');
legend({'Estimate', 'Actual'}, 'Location', 'best');

%% 局部函数
function [Y_r, Y_l] = build_regressor(q, dq, dq_r, ddq_r, params)
% 构造关节空间自适应控制回归矩阵。

c1 = cos(q(1));
c2 = cos(q(2));
s2 = sin(q(2));
c12 = cos(q(1)+q(2));

Y_r = zeros(2, 5);
Y_r(1, 1) = params.a(1) * ddq_r(1) + params.g * c1;
Y_r(1, 2) = ddq_r(1);
Y_r(1, 4) = ddq_r(1) + ddq_r(2);
Y_r(1, 3) = params.a(1) * c2 * (ddq_r(1) + Y_r(1, 4)) ...
    -params.a(1) * s2 * (dq(2) * dq_r(1) + dq_r(2) * (dq(1) + dq(2))) ...
    +params.a(2) * Y_r(1, 4) + params.g * c12;
Y_r(1, 5) = params.k_r2 * ddq_r(2) + Y_r(1, 2);
Y_r(2, 3) = params.a(2) * Y_r(1, 4) + params.a(1) * c2 * ddq_r(1) ...
    +params.a(1) * s2 * dq(1) * dq_r(1) + params.g * c12;
Y_r(2, 4) = Y_r(1, 4);
Y_r(2, 5) = params.k_r2 * Y_r(1, 5);

Y_l = [ ...
    params.a(1) * Y_r(1, 1) + params.a(2) * Y_r(1, 3); ...
    params.a(1) * Y_r(2, 1) + params.a(2) * Y_r(2, 3)];
end
