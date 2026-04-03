clearvars;
close all;
clc;

%% 路径设置
this_dir = fileparts(mfilename('fullpath'));
matlab_root = fileparts(this_dir);
addpath(this_dir);
addpath(matlab_root);

% 若缺少自动生成文件，则先执行符号推导脚本。
if exist(fullfile(this_dir, 'rhs_planar_contact.m'), 'file') ~= 2 || ...
        exist(fullfile(this_dir, 'energy_planar_contact.m'), 'file') ~= 2
    derive_planar_contact();
end

%% 模型参数
m1 = 1.0;
m2 = 0.8;
l1 = 1.0;
l2 = 0.8;
a1 = l1 / 2;
a2 = l2 / 2;
J1 = (1 / 12) * m1 * l1^2;
J2 = (1 / 12) * m2 * l2^2;
g = 9.81;

%% 地面接触参数
ground_y = -0.85;
k_ground = 1500;
c_ground = 35;

%% 初始状态
q1 = pi / 2 - 0.25;
u1 = 0;
q2 = -pi / 6;
u2 = 0;
z0 = [q1; u1; q2; u2];

%% 仿真设置
T = 6;
fs = 200;
tspan = linspace(0, T, T * fs + 1);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9, 'MaxStep', 1 / fs);

%% 动力学仿真
[t, z] = ode45(@(t, z) rhs_planar_contact(t, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ...
    ground_y, k_ground, c_ground), tspan, z0, options);

%% 运动学后处理
q1 = z(:, 1);
u1 = z(:, 2);
q2 = z(:, 3);
u2 = z(:, 4);

x1 = l1 * cos(q1);
y1 = l1 * sin(q1);
x2 = x1 + l2 * cos(q1 + q2);
y2 = y1 + l2 * sin(q1 + q2);

%% 能量诊断
n = numel(t);
KE = zeros(n, 1);
PEg = zeros(n, 1);
PEc = zeros(n, 1);
Pdamp = zeros(n, 1);
Fy = zeros(n, 1);
penetration = zeros(n, 1);
tip_y = zeros(n, 1);
tip_y_dot = zeros(n, 1);

for idx = 1:n
    [KE(idx), PEg(idx), PEc(idx), Pdamp(idx), Fy(idx), penetration(idx), ...
        tip_y(idx), tip_y_dot(idx)] = energy_planar_contact(t(idx), z(idx, :), ...
        m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground);
end

E_mech = KE + PEg;
E_diss = cumtrapz(t, Pdamp);
E_balance = E_mech + PEc + E_diss;
E_balance_drift = E_balance - E_balance(1);

%% 诊断绘图
utils.setDefaultGraphics;
utils.createFigureA4(struct('Name', 'Planar Contact Diagnostics', 'Width', 18, 'AspectRatio', 0.9));
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, q1, 'r');
hold on;
plot(t, q2, 'b');
grid on;
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('$q_1$', '$q_2$', 'Location', 'best');
title('Joint Angles');

nexttile;
yyaxis left;
plot(t, Fy, 'k');
ylabel('Contact Force (N)');
yyaxis right;
plot(t, penetration, 'm');
ylabel('Penetration (m)');
grid on;
xlabel('Time (s)');
title('Ground Contact');

nexttile;
plot(t, E_mech, 'LineWidth', 1.2);
hold on;
plot(t, PEc, 'LineWidth', 1.2);
plot(t, E_diss, 'LineWidth', 1.2);
plot(t, E_balance, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Energy (J)');
legend('Mechanical', 'Contact Spring', 'Damping Loss', 'Balance Sum', 'Location', 'best');
title('Energy Components');

nexttile;
plot(t, E_balance_drift, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Energy Drift (J)');
title('Energy Balance Drift');

%% 动画
all_x = [0; x1; x2];
all_y = [ground_y; y1; y2];
x_margin = 0.1 * (max(all_x(:)) - min(all_x(:)) + eps);
y_margin = 0.1 * (max(all_y(:)) - min(all_y(:)) + eps);
x_lim = [min(all_x(:)) - x_margin, max(all_x(:)) + x_margin];
y_lim = [min(all_y(:)) - y_margin, max(all_y(:)) + y_margin];

figure('Color', 'w', 'Name', 'Planar Two-Link Contact Animation');
axis equal;
axis([x_lim, y_lim]);
grid on;
hold on;
xlabel('x (m)');
ylabel('y (m)');
title('Planar Two-Link with Ground Contact');

plot(x_lim, [ground_y, ground_y], 'k--', 'LineWidth', 1.0);
h_link1 = plot([0, x1(1)], [0, y1(1)], 'r', 'LineWidth', 2);
h_link2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'b', 'LineWidth', 2);
h_tip = plot(x2(1), y2(1), 'ko', 'MarkerFaceColor', 'k');
traj_tip = animatedline('Color', [0.1, 0.4, 0.8], 'LineWidth', 1.0, 'MaximumNumPoints', 1200);

frame_stride = max(floor(n / 500), 1);
for idx = 1:frame_stride:n
    set(h_link1, 'XData', [0, x1(idx)], 'YData', [0, y1(idx)]);
    set(h_link2, 'XData', [x1(idx), x2(idx)], 'YData', [y1(idx), y2(idx)]);
    set(h_tip, 'XData', x2(idx), 'YData', y2(idx));

    % 触地时用红色标记末端。
    if penetration(idx) > 0
        set(h_tip, 'Color', [0.85, 0.2, 0.2], 'MarkerFaceColor', [0.85, 0.2, 0.2]);
    else
        set(h_tip, 'Color', 'k', 'MarkerFaceColor', 'k');
    end

    addpoints(traj_tip, x2(idx), y2(idx));
    title(sprintf('Planar Two-Link with Ground Contact, t = %.2f s, F_n = %.2f N', t(idx), Fy(idx)));
    drawnow limitrate;
end
