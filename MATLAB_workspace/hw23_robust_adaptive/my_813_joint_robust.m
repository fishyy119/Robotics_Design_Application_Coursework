%% 8.13 关节空间鲁棒控制纯 MATLAB 仿真
% 基于 hw23 参考程序中的两连杆动力学参数与鲁棒控制定义，
% 对末端附加 10 kg 集中负载的平面机械臂进行离散时域仿真。

clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 参数

Ts = 1e-3;
t_move = 1.0;
t_sim = 2;
t = 0:Ts:t_sim;
q_i = [0; pi / 4];
q_f = [pi / 2; pi / 2];

% 控制器只掌握无负载标称模型，但已知末端负载不超过 15 kg。
params = build_model(10.0, 15.0);
traj = build_joint_trajectory(struct( ...
    'q_i', q_i, ...
    'q_f', q_f, ...
    't_move', t_move, ...
    'Ts', Ts, ...
    't_total', t_sim, ...
    't_acc_ratio', 0.25));

K_p = 25 * eye(2);
K_d = 5 * eye(2);

rho = 40;
ep = 0.008;

D = [zeros(2, 2); eye(2, 2)];
H = [D'; -K_p, -K_d];
Q = lyap(H', eye(4));
DQ = D' * Q;

% 按已知负载上界构造惯性矩阵上界，供鲁棒项做保守设计。
% 这里使用 pi_bar 而不是实际对象参数 pi_l
B_h = diag([ ...
    params.a(1) * params.pi_bar(1) + params.pi_bar(2) + ...
    (params.a(2) + 2 * params.a(1)) * params.pi_bar(3) + params.pi_bar(4); ...
    params.a(2) * params.pi_bar(3) + params.pi_bar(4) + params.k_r2^2 * params.pi_bar(5)]);

%% 初始条件
N = numel(t);
q = zeros(2, N);
dq = zeros(2, N);
ddq = zeros(2, N);
tau = zeros(2, N);

q(:, 1) = q_i;

%% 离散仿真
for k = 1:N - 1
    q_err = traj.q(:, k) - q(:, k);
    dq_err = traj.dq(:, k) - dq(:, k);

    xi = [q_err; dq_err];
    z = DQ * xi;
    z_norm = norm(z);
    if z_norm >= ep
        w = rho / z_norm * z;
    else
        w = rho / ep * z;
    end

    y = w + traj.ddq(:, k) + K_p * q_err + K_d * dq_err;
    tau(:, k) = B_h * y + params.F_v * dq(:, k) + ...
        gravity_vector(q(:, k), params.pi_m, params.g);

    x_k = [q(:, k); dq(:, k)];
    x_kp1 = rk4_step(x_k, tau(:, k), Ts, params);
    q(:, k+1) = x_kp1(1:2);
    dq(:, k+1) = x_kp1(3:4);

    x_dot = manipulator_rhs(x_k, tau(:, k), params);
    ddq(:, k) = x_dot(3:4);
end

tau(:, N) = tau(:, N-1);
ddq(:, N) = ddq(:, N-1);

%% 后处理
[p_d_x, p_d_y] = forward_kinematics(traj.q, params.a);
[p_x, p_y] = forward_kinematics(q, params.a);
position_error_norm = sqrt((p_d_x - p_x).^2+(p_d_y - p_y).^2);
%% 绘图
fig1 = utils.createFigureA4(struct('Name', '8.13 Robust Control', 'Width', 18, 'AspectRatio', 0.70));
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
