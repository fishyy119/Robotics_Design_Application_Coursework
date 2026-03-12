clc;
clear;
close all;

% 定义不同仿真条件列表
simParams = [; ...
    struct('m', 1.0, 'g', 9.8, 'k_collision', 0.6, 'k_friction', 0.4, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    ];

% 循环执行仿真
for i = 1:numel(simParams)
    params = simParams(i);
    [t, z] = simulateBall(params);
    plotBallResults(t, z, params, i);
end
%% 仿真函数
function [t_total, z_total] = simulateBall(params)
% 参数解构
m = params.m;
g = params.g;
k_collision = params.k_collision;
k_friction = params.k_friction;
x0 = params.x0;
y0 = params.y0;
xd0 = params.xd0;
yd0 = params.yd0;

z0 = [x0, y0, xd0, yd0];
t0 = 0;
dt = 10;
N_time = 10000;

t_total = t0;
z_total = z0;

% 高精度 ODE 选项
options = odeset('abstol', 2.25e-14, 'reltol', 2.25e-14, 'events', @(t, z) collisionEvent(t, z));

%% 弹跳阶段
while true
    tspan = linspace(t0, t0+dt, N_time);
    [t_temp, z_temp, ~] = ode113(@(t, z) flyingBall(t, z, g), tspan, z0, options);

    % 碰撞处理
    z0 = collisionBall(z_temp(end, :), k_collision);

    t0 = t_temp(end) + dt / N_time;
    t_total = [t_total; t_temp(2:end); t0]; %#ok<AGROW>
    z_total = [z_total; z_temp(2:end, :); z0]; %#ok<AGROW>

    if abs(z_total(end, 4)) < 1e-4
        break;
    end
end

%% 滑动阶段
t0 = t_total(end) + dt / N_time;
z0 = z_total(end, :);
tspan = linspace(t0, t0+dt, N_time);
options = odeset('abstol', 2.25e-14, 'reltol', 2.25e-14);

[t_temp, z_temp] = ode113(@(t, z) slipBall(t, z, g, k_friction), tspan, z0, options);
t_total = [t_total; t_temp(2:end)];
z_total = [z_total; z_temp(2:end, :)];
end
%% 飞行阶段动力学
function zdot = flyingBall(~, z, g)
xd = z(3);
yd = z(4);
xdd = 0;
ydd = -g;
zdot = [xd; yd; xdd; ydd];
end

%% 碰撞处理
function zplus = collisionBall(z, k_collision)
zplus = [z(1), z(2), z(3), -k_collision * z(4)];
end

%% 滑动阶段动力学
function zdot = slipBall(~, z, g, k_friction)
xd = z(3);
yd = z(4);
if abs(xd) > 1e-6
    xdd = -sign(xd) * k_friction * g;
else
    xdd = 0;
end
zdot = [xd; yd; xdd; 0];
end

%% 碰撞事件函数
function [gstop, isterminal, direction] = collisionEvent(~, z)
gstop = z(2); % y = 0 碰撞触发
isterminal = 1; % 停止 ODE
direction = -1; % 仅下降触发
end

%% 绘图函数
function plotBallResults(t, z, params, idx)
utils.createFigureA4(struct('Name', sprintf('Simulation %d', idx), 'AspectRatio', 0.5))
utils.setDefaultGraphics;

% 水平、垂直速度
vx = z(:, 3);
vy = z(:, 4);

% x-y轨迹
subplot(2, 1, 1)
plot(z(:, 1), z(:, 2), 'b-')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
title(sprintf('Ball trajectory ($k_c$=%.2f, $k_f$=%.2f)', params.k_collision, params.k_friction))
grid on

% 速度随时间（x、y分量）
subplot(2, 1, 2)
plot(t, vx, 'r--'); hold on
plot(t, vy, 'b-')
xlabel('$t$ (s)')
ylabel('Velocity (m/s)')
legend({'$v_x$', '$v_y$'}, 'Location', 'best')
title('Velocity components vs time')
grid on
end
