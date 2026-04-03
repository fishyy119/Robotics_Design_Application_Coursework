%% 平面二连杆 PD 控制末端轨迹跟踪示例
% 轨迹在 t < t1 时为直线，在 t >= t1 时为圆弧。

clearvars;
close all;
clc;

%% 参数
m1 = 1;
m2 = 1;
model.L1 = 0.5;
model.L2 = 0.5;
l1 = model.L1;
l2 = model.L2;
g = 9.81;

%% 初始条件
q1 = -pi / 2;
u1 = 0;
q2 = 0;
u2 = 0.0001;
z0 = [q1; u1; q2; u2];

%% 时间设置
T = 8;
dt = 0.05;
tspan = 0:dt:T;
fs = 1 / dt;

%% PD 控制参数
Kp = diag([500, 500]);
Kd = diag([10, 10]);

global tau_record t_record
tau_record = [];
t_record = [];

%% 生成期望轨迹
w = pi / 2;
t1 = 2;

qd_traj = zeros(2, length(tspan));
qd_dot_traj = zeros(2, length(tspan));
x_traj = zeros(1, length(tspan));
y_traj = zeros(1, length(tspan));

for k = 1:length(tspan)
    t = tspan(k);
    if t < t1
        [x, ~, ~] = TSpline(0, 0, 0, 0, t1 / 2, 0, 0.25 * w, t1, t);
        [y, ~, ~] = TSpline(-1, 0, 0, -0.9, t1 / 2, -0.75, 0, t1, t);
    else
        x = 0.25 * sin(w * (t - t1));
        y = -0.5 - 0.25 * cos(w * (t - t1));
    end

    x_traj(k) = x;
    y_traj(k) = y;

    [q1d, q2d] = inverse_kinematics(x, y, model);
    qd_traj(:, k) = [q1d; q2d];

    if k == 1
        qd_dot_traj(:, k) = [0; 0];
    else
        qd_dot_traj(:, k) = (qd_traj(:, k) - qd_traj(:, k - 1)) / dt;
    end
end

qd_fun = @(t) [interp1(tspan, qd_traj(1, :), t); interp1(tspan, qd_traj(2, :), t)];
qd_dot_fun = @(t) [interp1(tspan, qd_dot_traj(1, :), t); interp1(tspan, qd_dot_traj(2, :), t)];

%% 求解
[t, z] = ode45(@(t, z) rhs_planar_PD(t, z, m1, m2, model.L1, model.L2, g, ...
    qd_fun, qd_dot_fun, Kp, Kd), tspan, z0);

%% 末端位置
x1 = l1 * cos(z(:, 1));
y1 = l1 * sin(z(:, 1));
x2 = x1 + l2 * cos(z(:, 1) + z(:, 3));
y2 = y1 + l2 * sin(z(:, 1) + z(:, 3));

%% 动画准备
all_x = [x1; x2];
all_y = [y1; y2];
x_margin = 0.05 * (max(all_x(:)) - min(all_x(:)));
y_margin = 0.05 * (max(all_y(:)) - min(all_y(:)));
x_lim = [min(all_x(:)) - x_margin, max(all_x(:)) + x_margin];
y_lim = [min(all_y(:)) - y_margin, max(all_y(:)) + y_margin];

figure('Color', 'w');
axis equal;
axis([x_lim, y_lim]);
axis off;
hold on;

% 初始化连杆与轨迹图元。
h_link1 = plot([0, x1(1)], [0, y1(1)], 'r', 'LineWidth', 2);
h_link2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'b', 'LineWidth', 2);
traj1 = animatedline('Color', [1, 0, 0, 0.3], 'LineWidth', 1, 'MaximumNumPoints', 1000);
traj2 = animatedline('Color', [0, 0, 1, 0.9], 'LineWidth', 1, 'MaximumNumPoints', 1000);
plot(x_traj, y_traj, 'k--', 'LineWidth', 1);

% 创建视频写入器。
video = VideoWriter('two_link_PD.mp4', 'MPEG-4');
video.FrameRate = fs;
open(video);

tic
for i = 1:length(t)
    set(h_link1, 'XData', [0, x1(i)], 'YData', [0, y1(i)]);
    set(h_link2, 'XData', [x1(i), x2(i)], 'YData', [y1(i), y2(i)]);
    addpoints(traj1, x1(i), y1(i));
    addpoints(traj2, x2(i), y2(i));
    drawnow;

    % 按仿真时间同步播放动画。
    elapsed = toc;
    target_time = t(i);
    if elapsed < target_time
        pause(target_time - elapsed);
    end

    % 写入当前视频帧。
    frame = getframe(gcf);
    writeVideo(video, frame);
end

close(video);

%% 关节力矩曲线
utils.createFigureA4();
utils.setDefaultGraphics;

plot(t_record, tau_record(1, :), 'r', 'LineWidth', 1.5);
hold on;
plot(t_record, tau_record(2, :), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Torque (N·m)');
legend('$\tau_1$', '$\tau_2$');
grid on;
title('Actual Joint Torques');

%% 局部函数
function [q1, q2] = inverse_kinematics(x, y, model)
% 计算平面二连杆的逆运动学解。

L1 = model.L1;
L2 = model.L2;
r = x^2 + y^2;
c2 = (r - L1^2 - L2^2) / (2 * L1 * L2);
q2 = acos(max(min(c2, 1), -1)); % 防止数值误差导致超出 [-1, 1]。
q1 = atan2(y, x) - atan2(L2 * sin(q2), L1 + L2 * cos(q2));
end
