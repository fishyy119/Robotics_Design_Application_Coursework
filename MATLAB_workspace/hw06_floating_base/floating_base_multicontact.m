clearvars;
close all;
clc;
utils.setDefaultGraphics;
%% 路径设置
this_dir = fileparts(mfilename('fullpath'));
matlab_root = fileparts(this_dir);
addpath(this_dir);
addpath(matlab_root);

regenerate_symbolic = true; % 调试阶段强制重生成功能文件，避免继续使用旧展开式
if regenerate_symbolic || ...
        exist(fullfile(this_dir, 'rhs_floating_base_multicontact.m'), 'file') ~= 2 || ...
        exist(fullfile(this_dir, 'energy_floating_base_multicontact.m'), 'file') ~= 2 || ...
        exist(fullfile(this_dir, 'draw_floating_base_multicontact.m'), 'file') ~= 2
    derive_floating_base_multicontact();
    clear rhs_floating_base_multicontact energy_floating_base_multicontact draw_floating_base_multicontact
    rehash;
end
%% 模型参数
m = 1.6;
lCuboid = 0.48;
wCuboid = 0.32;
hCuboid = 0.24;

xG = 0.04;
yG = -0.03;
zG = 0.02;

I_body = diag([ ...
    m * (wCuboid^2 + hCuboid^2) / 12, ...
    m * (lCuboid^2 + hCuboid^2) / 12, ...
    m * (lCuboid^2 + wCuboid^2) / 12]);

g = 5;
%% 地面接触参数
ground_z = 0.0;
k_ground = 1800;
c_ground = 10;
%% 初始条件
q1 = -0.20;
u1 = 0.35;
q2 = 0.15;
u2 = -0.18;
q3 = 0.95;
u3 = 0.00;
q4 = 0.32;
u4 = 0.90;
q5 = -0.28;
u5 = -0.65;
q6 = 0.45;
u6 = 0.75;
E_diss_damping0 = 0;

z0 = [q1; u1; q2; u2; q3; u3; q4; u4; q5; u5; q6; u6; E_diss_damping0];
%% 数值积分
T = 8.0;
fs = 600;
tspan = linspace(0, T, T*fs+1);
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-11, 'MaxStep', 1/fs);

% 接触弹簧阻尼会引入刚性，采用刚性求解器提高积分稳定性。
[t, z] = ode15s(@(t, z) rhs_floating_base_multicontact(t, z, m, xG, yG, zG, I_body, g, ...
    lCuboid, wCuboid, hCuboid, ground_z, k_ground, c_ground), tspan, z0, options);
%% 后处理
n_step = numel(t);
r_B_all = zeros(3, n_step);
r_G_all = zeros(3, n_step);
B_axes_all = zeros(3, 3, n_step);
G_axes_all = zeros(3, 3, n_step);
vertices_all = zeros(3, 8, n_step);
penetration_all = zeros(8, n_step);
vertex_z_all = zeros(8, n_step);

KE = zeros(n_step, 1);
PEg = zeros(n_step, 1);
PEc = zeros(n_step, 1);
E_diss_damping = z(:, 13);

for idx = 1:n_step
    [r_B_all(:, idx), B_axes_all(:, :, idx), r_G_all(:, idx), G_axes_all(:, :, idx), ...
        vertices_all(:, :, idx)] = draw_floating_base_multicontact(t(idx), z(idx, :), ...
        xG, yG, zG, lCuboid, wCuboid, hCuboid);

    [KE(idx), PEg(idx), PEc(idx), ~, penetration_all(:, idx), vertex_z_all(:, idx), ~] = ...
        energy_floating_base_multicontact(t(idx), z(idx, :), m, xG, yG, zG, I_body, g, ...
        lCuboid, wCuboid, hCuboid, ground_z, k_ground, c_ground);
end

contact_count = sum(penetration_all > 0, 1).';
E_balance = KE + PEg + PEc + E_diss_damping;
%% 能量诊断图
fig_energy = utils.createFigureA4(struct('Name', 'Floating Base Energy', 'Width', 18, 'AspectRatio', 0.75));
figure(fig_energy);
plot(t, KE, 'LineWidth', 1.2);
hold on;
plot(t, PEg, 'LineWidth', 1.2);
plot(t, PEc, 'LineWidth', 1.2);
plot(t, E_diss_damping, 'LineWidth', 1.2);
plot(t, E_balance, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Energy (J)');
legend('Kinetic', 'Gravity', 'Spring', 'Damping Loss', 'Balance', 'Location', 'best');
title('Energy Components');
%% 动画保存
% 动画显示选项
display_options = struct( ...
    'show_world_frame', true, ...
    'show_body_com_frame', true, ...
    'show_body_ref_frame', false);


faces = [ ...
    1, 2, 3, 4; ...
    1, 4, 8, 5; ...
    3, 4, 8, 7; ...
    5, 6, 7, 8; ...
    2, 3, 7, 6; ...
    1, 2, 6, 5];

face_colors = [ ...
    0.80, 0.27, 0.25; ...
    0.26, 0.55, 0.78; ...
    0.20, 0.58, 0.41; ...
    0.85, 0.57, 0.18; ...
    0.56, 0.44, 0.74; ...
    0.37, 0.37, 0.37];

all_x = reshape(vertices_all(1, :, :), [], 1);
all_y = reshape(vertices_all(2, :, :), [], 1);
all_z = reshape(vertices_all(3, :, :), [], 1);

x_margin = 0.12 * (max(all_x) - min(all_x) + eps);
y_margin = 0.12 * (max(all_y) - min(all_y) + eps);
z_margin = 0.12 * (max(all_z) - min(all_z) + eps);

x_lim = [min(all_x) - x_margin, max(all_x) + x_margin];
y_lim = [min(all_y) - y_margin, max(all_y) + y_margin];
z_lim = [min(ground_z, min(all_z)) - z_margin, max(all_z) + z_margin];

fig_anim = figure( ...
    'Color', 'w', ...
    'Name', 'Floating Base Multicontact Animation', ...
    'NumberTitle', 'off', ...
    'Units', 'pixels', ...
    'Position', [60, 60, 980, 760], ...
    'Renderer', 'opengl');

ax = axes('Parent', fig_anim);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
axis(ax, [x_lim, y_lim, z_lim]);
view(ax, 32, 22);
xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
zlabel(ax, 'z (m)');
title_handle = title(ax, sprintf('t = %.2f s, Active Contacts = %d', t(1), contact_count(1)));
axis_colors = [0.80, 0.15, 0.15; 0.15, 0.35, 0.85; 0.12, 0.60, 0.25];
world_axis_length = 0.35;

visibility = struct( ...
    'world', logical_to_on_off(display_options.show_world_frame), ...
    'body_com', logical_to_on_off(display_options.show_body_com_frame), ...
    'body_ref', logical_to_on_off(display_options.show_body_ref_frame));

ground_patch = patch(ax, ...
    'XData', [x_lim(1), x_lim(2), x_lim(2), x_lim(1)], ...
    'YData', [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
    'ZData', ground_z*ones(1, 4), ...
    'FaceColor', [0.88, 0.88, 0.88], ...
    'FaceAlpha', 0.45, ...
    'EdgeColor', [0.65, 0.65, 0.65], ...
    'LineWidth', 1.0); %#ok<NASGU>

h_world_axes = gobjects(3, 1);
for axis_idx = 1:3
    world_end = zeros(3, 1);
    world_end(axis_idx) = world_axis_length;
    h_world_axes(axis_idx) = plot3(ax, ...
        [0, world_end(1)], ...
        [0, world_end(2)], ...
        [0, world_end(3)], ...
        'Color', axis_colors(axis_idx, :), ...
        'LineWidth', 1.6, ...
        'Visible', visibility.world);
end

h_faces = gobjects(size(faces, 1), 1);
for face_idx = 1:size(faces, 1)
    verts = vertices_all(:, faces(face_idx, :), 1);
    h_faces(face_idx) = patch(ax, ...
        'XData', verts(1, :), ...
        'YData', verts(2, :), ...
        'ZData', verts(3, :), ...
        'FaceColor', face_colors(face_idx, :), ...
        'FaceAlpha', 0.78, ...
        'EdgeColor', [0.15, 0.15, 0.15], ...
        'LineWidth', 1.0);
end

h_B_axes = gobjects(3, 1);
h_G_axes = gobjects(3, 1);

for axis_idx = 1:3
    h_B_axes(axis_idx) = plot3(ax, ...
        [r_B_all(1, 1), B_axes_all(1, axis_idx, 1)], ...
        [r_B_all(2, 1), B_axes_all(2, axis_idx, 1)], ...
        [r_B_all(3, 1), B_axes_all(3, axis_idx, 1)], ...
        '--', 'Color', axis_colors(axis_idx, :), 'LineWidth', 1.5, ...
        'Visible', visibility.body_ref);

    h_G_axes(axis_idx) = plot3(ax, ...
        [r_G_all(1, 1), G_axes_all(1, axis_idx, 1)], ...
        [r_G_all(2, 1), G_axes_all(2, axis_idx, 1)], ...
        [r_G_all(3, 1), G_axes_all(3, axis_idx, 1)], ...
        '-', 'Color', axis_colors(axis_idx, :), 'LineWidth', 2.1, ...
        'Visible', visibility.body_com);
end

h_B_point = plot3(ax, r_B_all(1, 1), r_B_all(2, 1), r_B_all(3, 1), 'ko', ...
    'MarkerFaceColor', [0.12, 0.12, 0.12], 'MarkerSize', 6, 'Visible', visibility.body_ref);
h_G_point = plot3(ax, r_G_all(1, 1), r_G_all(2, 1), r_G_all(3, 1), 's', ...
    'MarkerFaceColor', [0.98, 0.82, 0.18], 'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
    'MarkerSize', 7, 'Visible', visibility.body_com);

vertex_colors = repmat([0.15, 0.15, 0.15], 8, 1);
vertex_colors(penetration_all(:, 1) > 0, :) = repmat([0.86, 0.18, 0.18], nnz(penetration_all(:, 1) > 0), 1);
h_vertices = scatter3(ax, vertices_all(1, :, 1), vertices_all(2, :, 1), vertices_all(3, :, 1), ...
    42, vertex_colors, 'filled', 'MarkerEdgeColor', [0.05, 0.05, 0.05]);

lighting(ax, 'gouraud');
camlight(ax, 'headlight');

video_path = fullfile(this_dir, 'floating_base_multicontact.mp4');
video = VideoWriter(video_path, 'MPEG-4');
video.FrameRate = 30;
open(video);
video_is_open = true;
frame_height = [];
frame_width = [];

scene = struct();
scene.faces = faces;
scene.handles = struct( ...
    'world_axes', h_world_axes, ...
    'faces', h_faces, ...
    'body_ref_axes', h_B_axes, ...
    'body_com_axes', h_G_axes, ...
    'body_ref_point', h_B_point, ...
    'body_com_point', h_G_point, ...
    'vertices', h_vertices, ...
    'title', title_handle);

try
    frame_stride = max(round(fs/video.FrameRate), 1);
    for idx = 1:frame_stride:n_step
        frame_data = struct( ...
            'vertices', vertices_all(:, :, idx), ...
            'r_B', r_B_all(:, idx), ...
            'B_axes', B_axes_all(:, :, idx), ...
            'r_G', r_G_all(:, idx), ...
            'G_axes', G_axes_all(:, :, idx), ...
            'penetration', penetration_all(:, idx), ...
            'time', t(idx), ...
            'contact_count', contact_count(idx));
        update_cuboid_scene(scene, frame_data);

        drawnow;
        frame = getframe(fig_anim);
        if isempty(frame_height)
            [frame_height, frame_width, ~] = size(frame.cdata);
        else
            frame = normalize_video_frame(frame, frame_height, frame_width);
        end
        writeVideo(video, frame);
    end

    close(video);
    video_is_open = false;
catch ME
    if video_is_open
        close(video);
    end
    rethrow(ME);
end

function update_cuboid_scene(scene, frame_data)
% 更新三维动画中的几何对象。

for face_idx = 1:size(scene.faces, 1)
    verts = frame_data.vertices(:, scene.faces(face_idx, :));
    set(scene.handles.faces(face_idx), ...
        'XData', verts(1, :), ...
        'YData', verts(2, :), ...
        'ZData', verts(3, :));
end

for axis_idx = 1:3
    set(scene.handles.body_ref_axes(axis_idx), ...
        'XData', [frame_data.r_B(1), frame_data.B_axes(1, axis_idx)], ...
        'YData', [frame_data.r_B(2), frame_data.B_axes(2, axis_idx)], ...
        'ZData', [frame_data.r_B(3), frame_data.B_axes(3, axis_idx)]);

    set(scene.handles.body_com_axes(axis_idx), ...
        'XData', [frame_data.r_G(1), frame_data.G_axes(1, axis_idx)], ...
        'YData', [frame_data.r_G(2), frame_data.G_axes(2, axis_idx)], ...
        'ZData', [frame_data.r_G(3), frame_data.G_axes(3, axis_idx)]);
end

set(scene.handles.body_ref_point, ...
    'XData', frame_data.r_B(1), 'YData', frame_data.r_B(2), 'ZData', frame_data.r_B(3));
set(scene.handles.body_com_point, ...
    'XData', frame_data.r_G(1), 'YData', frame_data.r_G(2), 'ZData', frame_data.r_G(3));

vertex_colors = repmat([0.15, 0.15, 0.15], 8, 1);
vertex_colors(frame_data.penetration > 0, :) = repmat([0.86, 0.18, 0.18], nnz(frame_data.penetration > 0), 1);
set(scene.handles.vertices, ...
    'XData', frame_data.vertices(1, :), ...
    'YData', frame_data.vertices(2, :), ...
    'ZData', frame_data.vertices(3, :), ...
    'CData', vertex_colors);

set(scene.handles.title, 'String', ...
    sprintf('t = %.2f s, Active Contacts = %d', frame_data.time, frame_data.contact_count));
end

function value = logical_to_on_off(flag)
% 将逻辑量转换为 MATLAB 图形对象可见性字符串。

if flag
    value = 'on';
else
    value = 'off';
end
end

function frame = normalize_video_frame(frame, target_height, target_width)
% 通过裁剪或补边保持视频帧尺寸不变。

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
