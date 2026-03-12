clear;
clc;

%% 参数
model.L1 = 0.5;
model.L2 = 0.5;

% 轨迹关键点
p0 = [-0.5; -0.5];
p1 = [0; -0.25];
p2 = [0.5; -0.5];

t0 = 0;
t1 = 2;
t2 = 4;
dt = 0.05;

% 开关：true -> 先关节角插值, false -> 先末端插值
joint_interp_mode = true;

%% 数据收集
t_snap = t0:dt:t2;
n_snap = length(t_snap);

snapshots = struct('y1', [], 'z1', [], 'y2', [], 'z2', [], 'color', []);
% cmap = utils.viridis(n_snap);
cmap = parula(n_snap);

if joint_interp_mode
    [q10, q20] = inverse_kinematics(p0(1), p0(2), model);
    [q11, q21] = inverse_kinematics(p1(1), p1(2), model);
    [q12, q22] = inverse_kinematics(p2(1), p2(2), model);
end
for k = 1:n_snap
    t = t_snap(k);
    if joint_interp_mode
        % 关节角插值
        q1 = TSpline(q10, 0, t0, q11, t1, q12, 0, t2, t);
        q2 = TSpline(q20, 0, t0, q21, t1, q22, 0, t2, t);
        [y1, z1, y2, z2] = forward_kinematics(q1, q2, model);
    else
        % 末端插值
        [y, yd, ydd] = TSpline(p0(1), 0, t0, p1(1), t1, p2(1), 0, t2, t);
        [z, zd, zdd] = TSpline(p0(2), 0, t0, p1(2), t1, p2(2), 0, t2, t);
        [q1, q2] = inverse_kinematics(y, z, model);
        [y1, z1, y2, z2] = forward_kinematics(q1, q2, model);
    end
    % 保存快照
    snapshots(k).y1 = y1;
    snapshots(k).z1 = z1;
    snapshots(k).y2 = y2;
    snapshots(k).z2 = z2;
    snapshots(k).color = cmap(k, :);
end
%% 绘图
utils.createFigureA4();
utils.setDefaultGraphics;

hold on
grid on
axis equal
xlim([-0.6, 0.6])
ylim([-0.6, 0.2])

if joint_interp_mode
    title('Trajectory (Joint Angle Interpolation)')
else
    title('Trajectory (End-Effector Interpolation)')
end

% 绘制所有连杆
for k = 1:n_snap
    if ~isempty(snapshots(k).y1)
        plot([0, snapshots(k).y1], [0, snapshots(k).z1], '-', 'Color', snapshots(k).color, 'LineWidth', 0.5);
        plot([snapshots(k).y1, snapshots(k).y2], [snapshots(k).z1, snapshots(k).z2], '-', 'Color', snapshots(k).color, 'LineWidth', 0.5);
    end
end

% 绘制所有端点在最上层
for k = 1:n_snap
    if ~isempty(snapshots(k).y1)
        plot(snapshots(k).y1, snapshots(k).z1, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
        plot(snapshots(k).y2, snapshots(k).z2, '.', 'MarkerSize', 8, 'Color', 'k');
    end
end

%%
function [q1, q2] = inverse_kinematics(y, z, model)
L1 = model.L1;
L2 = model.L2;

r2 = y^2 + z^2;
q2 = acos((r2 - L1^2 - L2^2)/(2 * L1 * L2));
q1 = atan2(y, -z) - acos((L1^2 + r2 - L2^2)/(2 * L1 * sqrt(r2)));
end

function [y1, z1, y2, z2] = forward_kinematics(q1, q2, model)
y1 = model.L1 * sin(q1);
z1 = -model.L1 * cos(q1);
y2 = y1 + model.L2 * sin(q1+q2);
z2 = z1 - model.L2 * cos(q1+q2);
end
