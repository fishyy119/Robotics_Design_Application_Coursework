clc;
clear;
close all;

%% 仿真参数
simParams = [ ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 1000, 'c_damp', 10, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 1000, 'c_damp', 5, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 250, 'c_damp', 5, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0) ...
    ];

nSim = numel(simParams);

%% 循环仿真
allResults = cell(nSim, 1);
for i = 1:nSim
    params = simParams(i);
    [t, z] = simulateSpringBounce(params);
    allResults{i} = struct('t', t, 'z', z, 'params', params);
end

plotSpringBounceResults(allResults);

%% 局部函数
function [t_total, z_total] = simulateSpringBounce(params)
% 计算弹簧阻尼接触模型下的小球运动。

m = params.m;
g = params.g;
k = params.k_spring;
c = params.c_damp;
x0 = params.x0;
y0 = params.y0;
xd0 = params.xd0;
yd0 = params.yd0;

z0 = [x0; y0; xd0; yd0];
t0 = 0;
dt = 12;
n_time = 10000;

t_total = t0;
z_total = z0';

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
tspan = linspace(t0, t0 + dt, n_time);
[t_temp, z_temp] = ode113(@(t, z) springBounceODE(t, z, m, g, k, c), tspan, z0, options);

t_total = [t_total; t_temp(2:end)];
z_total = [z_total; z_temp(2:end, :)];
end

function zdot = springBounceODE(~, z, m, g, k, c)
% 计算弹簧阻尼接触模型的状态导数。

y = z(2);
xd = z(3);
yd = z(4);

% 仅在触地后施加竖直方向弹簧阻尼力。
if y <= 0
    F_spring = -k * y; % y < 0 表示弹簧被压缩。
    F_damp = -c * yd;
else
    F_spring = 0;
    F_damp = 0;
end

ydd = -g + (F_spring + F_damp) / m;
xdd = 0; % 水平方向保持匀速。
zdot = [xd; yd; xdd; ydd];
end

function plotSpringBounceResults(results)
% 绘制多组参数下的小球轨迹。

nSim = numel(results);
utils.createFigureA4(struct('AspectRatio', 1));

for i = 1:nSim
    res = results{i};
    z = res.z;
    params = res.params;

    % 绘制 x-y 轨迹。
    subplot(nSim, 1, i)
    plot(z(:, 1), z(:, 2), 'b-')
    xlabel('$x$ (m)')
    ylabel('$y$ (m)')
    title(sprintf('Ball trajectory ($k_s$=%.1f N/m, $c$=%.2f kg/s)', params.k_spring, params.c_damp))
    xlim([z(1, 1), z(end, 1)])
    grid on
end
end
