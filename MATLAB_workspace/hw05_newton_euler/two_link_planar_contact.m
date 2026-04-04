clearvars;
close all;
clc;
utils.setDefaultGraphics;

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
m1 = 1.05;
m2 = 0.8;
l1 = 1.0;
l2 = 0.8;
a1 = l1 / 2;
a2 = l2 / 2;
J1 = (1 / 12) * m1 * l1^2;
J2 = (1 / 12) * m2 * l2^2;
g = 9.81;

%% 地面接触参数
ground_y = -1.05;
k_ground = 1500;
c_ground = 25;
mu_ground = 0.3;
v_smooth = 0.02; % 使用tanh平滑sign函数，这个是参考速度

%% 初始状态
q1 = pi / 2 - 0.25;
u1 = 0;
q2 = -pi / 6;
u2 = 0;
z0 = [q1; u1; q2; u2; 0; 0];

%% 仿真设置
T = 6;
fs = 200;
tspan = linspace(0, T, T*fs+1);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-9, 'MaxStep', 1/fs);
%% 动画窗口
anim_window_pos = [120, 80];
anim_window_size = [300, 300];
%% 动力学仿真
% 接触刚度与库仑摩擦会引入刚性，采用刚性求解器提高积分效率。
[t, z] = ode15s(@(t, z) rhs_planar_contact(t, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ...
    ground_y, k_ground, c_ground, mu_ground, v_smooth), tspan, z0, options);

%% 运动学后处理
q1 = z(:, 1);
u1 = z(:, 2);
q2 = z(:, 3);
u2 = z(:, 4);
E_diss_damping = z(:, 5);
E_diss_friction = z(:, 6);

x1 = l1 * cos(q1);
y1 = l1 * sin(q1);
x2 = x1 + l2 * cos(q1+q2);
y2 = y1 + l2 * sin(q1+q2);

%% 能量诊断
n = numel(t);
KE = zeros(n, 1);
PEg = zeros(n, 1);
PEc = zeros(n, 1);
Fy = zeros(n, 1);
penetration = zeros(n, 1);
tip_y = zeros(n, 1);
tip_y_dot = zeros(n, 1);

for idx = 1:n
    [KE(idx), PEg(idx), PEc(idx), ~, ~, Fy(idx), penetration(idx), ...
        tip_y(idx), tip_y_dot(idx)] = energy_planar_contact(t(idx), z(idx, :), ...
        m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground, mu_ground, v_smooth);
end

E_mech = KE + PEg;
E_balance = E_mech + PEc + E_diss_damping + E_diss_friction;

%% 诊断绘图
utils.createFigureA4();
plot(t, E_mech, 'LineWidth', 1.2);
hold on;
plot(t, PEc, 'LineWidth', 1.2);
plot(t, E_diss_damping, 'LineWidth', 1.2);
plot(t, E_diss_friction, 'LineWidth', 1.2);
plot(t, E_balance, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Energy (J)');
legend('Mechanical', 'Spring', 'Damping Loss', 'Friction Loss', 'Balance', 'Location', 'best');
title('Energy Components');

%% 动画
all_x = [0; x1; x2];
all_y = [ground_y; y1; y2];
x_margin = 0.1 * (max(all_x(:)) - min(all_x(:)) + eps);
y_margin = 0.1 * (max(all_y(:)) - min(all_y(:)) + eps);
x_lim = [min(all_x(:)) - x_margin, max(all_x(:)) + x_margin];
y_lim = [min(all_y(:)) - y_margin, max(all_y(:)) + y_margin];

fig_anim = figure( ...
    'Color', 'w', ...
    'Name', 'Planar Two-Link Contact Animation', ...
    'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [anim_window_pos, anim_window_size], ...
    'Resize', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none');
axis equal;
axis([x_lim, y_lim]);
grid on;
hold on;
xlabel('x (m)');
ylabel('y (m)');
title_fmt = '$k = %.0f\\,\\mathrm{N/m}, c_n = %.1f\\,\\mathrm{N\\cdot s/m}, \\mu = %.2f \\\\ t = %.2f\\,\\mathrm{s}, F_n = %.2f\\,\\mathrm{N}$';
h_title = title(sprintf(title_fmt, k_ground, c_ground, mu_ground, t(1), Fy(1)), 'Interpreter', 'latex');

plot(x_lim, [ground_y, ground_y], 'k--', 'LineWidth', 1.0);
h_link1 = plot([0, x1(1)], [0, y1(1)], 'r', 'LineWidth', 2);
h_link2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'b', 'LineWidth', 2);
h_tip = plot(x2(1), y2(1), 'ko', 'MarkerFaceColor', 'k');
traj_tip = animatedline('Color', [0.1, 0.4, 0.8], 'LineWidth', 1.0, 'MaximumNumPoints', 1200);

% 保存动画视频。
video_path = fullfile(this_dir, 'two_link_planar_contact.mp4');
video = VideoWriter(video_path, 'MPEG-4');
video.FrameRate = 50;
frame_height = [];
frame_width = [];
open(video);
video_is_open = true;

try
    frame_stride = max(floor(n / 500), 1);
    tic;
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
        set(h_title, 'String', sprintf(title_fmt, k_ground, c_ground, mu_ground, t(idx), Fy(idx)));
        drawnow;

        % 按仿真时间同步播放动画。
        elapsed = toc;
        target_time = t(idx) - t(1);
        if elapsed < target_time
            pause(target_time-elapsed);
        end

        frame = getframe(fig_anim);
        if isempty(frame_height)
            [frame_height, frame_width, ~] = size(frame.cdata);
        else
            frame = normalizeVideoFrame(frame, frame_height, frame_width);
        end
        writeVideo(video, frame);
    end

    % 显式关闭写入器，确保 MP4 尾信息落盘。
    close(video);
    video_is_open = false;
catch ME
    if video_is_open
        close(video);
    end
    rethrow(ME);
end

function frame = normalizeVideoFrame(frame, target_height, target_width)
% 通过裁剪或补边保证视频帧尺寸恒定。

cdata = frame.cdata;
[height, width, channels] = size(cdata);

if height == target_height && width == target_width
    return;
end

normalized_cdata = 255 * ones(target_height, target_width, channels, 'uint8');
copy_height = min(height, target_height);
copy_width = min(width, target_width);

src_row = floor((height - copy_height)/2) + 1;
src_col = floor((width - copy_width)/2) + 1;
dst_row = floor((target_height - copy_height)/2) + 1;
dst_col = floor((target_width - copy_width)/2) + 1;

normalized_cdata(dst_row:(dst_row + copy_height - 1), dst_col:(dst_col + copy_width - 1), :) = ...
    cdata(src_row:(src_row + copy_height - 1), src_col:(src_col + copy_width - 1), :);

frame.cdata = normalized_cdata;
end
