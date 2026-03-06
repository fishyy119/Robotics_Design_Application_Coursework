clc;
clear;
close all;

%% 仿真参数列表（弹簧-阻尼碰撞）
simParams = [ ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 1000, 'c_damp', 10, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 1000, 'c_damp', 5, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    struct('m', 1.0, 'g', 9.8, 'k_spring', 250, 'c_damp', 5, 'x0', 0, 'y0', 10, 'xd0', 5, 'yd0', 0), ...
    ];

nSim = numel(simParams);

%% 循环执行仿真
allResults = cell(nSim, 1);
for i = 1:nSim
    params = simParams(i);
    [t, z] = simulateSpringBounce(params);
    allResults{i} = struct('t', t, 'z', z, 'params', params);
end
plotSpringBounceResults(allResults)
%% 仿真函数（弹簧-阻尼碰撞）
function [t_total, z_total] = simulateSpringBounce(params)
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
N_time = 10000;


t_total = t0;
z_total = z0';

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

tspan = linspace(t0, t0+dt, N_time);
[t_temp, z_temp] = ode113(@(t, z) springBounceODE(t, z, m, g, k, c), tspan, z0, options);

t_total = [t_total; t_temp(2:end)];
z_total = [z_total; z_temp(2:end, :)];

end

%% 弹簧-阻尼动力学
function zdot = springBounceODE(~, z, m, g, k, c)
x = z(1);
y = z(2);
xd = z(3);
yd = z(4);

% 仅垂直方向弹簧-阻尼，触地后受力
if y <= 0
    F_spring = -k * y; % y<0 表示压缩弹簧
    F_damp = -c * yd; % 阻尼
else
    F_spring = 0;
    F_damp = 0;
end

ydd = -g + (F_spring + F_damp) / m;
xdd = 0; % 水平匀速
zdot = [xd; yd; xdd; ydd];
end

function plotSpringBounceResults(results)

nSim = numel(results);
utils.createFigureA4(struct('AspectRatio', 1))

for i = 1:nSim
    res = results{i};
    t = res.t;
    z = res.z;
    params = res.params;

    % x-y轨迹
    subplot(nSim, 1, i)
    plot(z(:, 1), z(:, 2), 'b-')
    xlabel('$x$ (m)'); ylabel('$y$ (m)')
    title(sprintf('Ball trajectory ($k_s$=%.1f N/m, $c$=%.2f kg/s)', params.k_spring, params.c_damp))
    xlim([z(1, 1), z(end, 1)])
    grid on

end

end
