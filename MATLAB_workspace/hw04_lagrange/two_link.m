clearvars;
close all;
clc;
%% 参数
m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
g = 9.81;
%% 初始条件
q1 = pi / 3;
u1 = 0;
q2 = pi / 2;
u2 = 0;

z0 = [q1; u1; q2; u2];
%% 时间
T = 5;
fs = 100;
tspan = linspace(0, T, T*fs);
%% 求解
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

[t, z] = ode45(@(t, z) rhs_planar(t, z, m1, m2, l1, l2, g), tspan, z0, options);
%% 数据整理（后续绘图用）
q1 = z(:, 1);
u1 = z(:, 2);
q2 = z(:, 3);
u2 = z(:, 4);

% 末端位置（可用于之后绘图）
x1 = l1 * cos(q1);
y1 = l1 * sin(q1);

x2 = x1 + l2 * cos(q1+q2);
y2 = y1 + l2 * sin(q1+q2);
%% ===== 动能（刚体） =====

a1 = l1 / 2;
a2 = l2 / 2;

I1 = (1 / 12) * m1 * l1^2;
I2 = (1 / 12) * m2 * l2^2;

% 质心速度平方
v1_sq = (a1 * u1).^2;

v2_sq = (l1 * u1).^2 + (a2 * (u1 + u2)).^2 + ...
    2 * l1 * a2 * u1 .* (u1 + u2) .* cos(q2);

% 平动动能
KE_trans = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq;

% 转动动能
KE_rot = 0.5 * I1 * u1.^2 + 0.5 * I2 * (u1 + u2).^2;

KE = KE_trans + KE_rot;
%% ===== 势能 =====

yc1 = a1 * sin(q1);
yc2 = l1 * sin(q1) + a2 * sin(q1+q2);

PE = m1 * g * yc1 + m2 * g * yc2;

% 总能量
TE = KE + PE;

% 能量差分
TE_diff = diff(TE);
t_diff = t(1:end-1);
%% ===== 末端轨迹 =====
all_x = [x1; x2];
all_y = [y1; y2];

x_min = min(all_x);
x_max = max(all_x);
y_min = min(all_y);
y_max = max(all_y);

% 留 5% margin
x_margin = 0.05 * (x_max - x_min);
y_margin = 0.05 * (y_max - y_min);

x_lim = [x_min - x_margin, x_max + x_margin];
y_lim = [y_min - y_margin, y_max + y_margin];

figure('Color', 'w'); % 白色背景
axis equal;
axis([x_lim, y_lim]);
axis off;
hold on;

% ===== 视频写入器 =====
video = VideoWriter('two_link.mp4', 'MPEG-4');
video.FrameRate = fs / 3;
open(video);

% ===== 杆 =====
h_link1 = plot([0, x1(1)], [0, y1(1)], 'r', 'LineWidth', 2);
h_link2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'b', 'LineWidth', 2);

% ===== 轨迹（半透明）=====
traj1 = animatedline('Color', [1, 0, 0, 0.3], 'LineWidth', 1);
traj2 = animatedline('Color', [0, 0, 1, 0.3], 'LineWidth', 1);

% （可选）限制轨迹长度
traj1.MaximumNumPoints = 1000;
traj2.MaximumNumPoints = 1000;

tic
for i = 1:length(t)

    % 更新杆
    set(h_link1, 'XData', [0, x1(i)], 'YData', [0, y1(i)]);
    set(h_link2, 'XData', [x1(i), x2(i)], 'YData', [y1(i), y2(i)]);

    % 更新轨迹
    addpoints(traj1, x1(i), y1(i));
    addpoints(traj2, x2(i), y2(i));

    drawnow;

    % ===== 时间同步（关键）=====
    elapsed = toc;
    target_time = t(i);
    if elapsed < target_time
        pause(target_time-elapsed);
    end

    % ===== 写入视频帧 =====
    frame = getframe(gcf);
    writeVideo(video, frame);

end

close(video);
%% ===== 能量守恒检查 =====
utils.createFigureA4();
utils.setDefaultGraphics;

plot(t_diff, TE_diff, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('$\Delta$ Energy');
title('Energy Difference');
grid on;
