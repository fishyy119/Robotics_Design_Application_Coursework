clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 电机参数
motor = getMotorParams('RE40 150W');

%% 负载设置
light_load.name = 'Light load';
light_load.J_load = 6.0e-4;
light_load.B = 1.2e-3;
light_load.tau_c = 0.012;
light_load.w_smooth = 1.0;
light_load.J_total = motor.J + light_load.J_load;

medium_load.name = 'Medium load';
medium_load.J_load = 1.2e-3;
medium_load.B = 1.5e-3;
medium_load.tau_c = 0.020;
medium_load.w_smooth = 1.0;
medium_load.J_total = motor.J + medium_load.J_load;

heavy_load.name = 'Heavy load';
heavy_load.J_load = 2.0e-3;
heavy_load.B = 1.8e-3;
heavy_load.tau_c = 0.030;
heavy_load.w_smooth = 1.0;
heavy_load.J_total = motor.J + heavy_load.J_load;

load_cases = [light_load, medium_load, heavy_load];
controller.mode = 'open_loop';
controller.voltageLimit = motor.U;
controller.voltageFcn = @(t, x) motor.U * double(t >= 0);
command = [];

%% 数值积分
x0 = [0; 0; 0];
tspan = linspace(0, 2.0, 2001);
results = cell(numel(load_cases), 1);

for idx = 1:numel(load_cases)
    results{idx} = simulateMotorResponse(tspan, x0, motor, load_cases(idx), controller, command);
end

%% 绘图
fig_params.Name = 'Open Loop Responses';
fig_params.Width = 18;
fig_params.AspectRatio = 0.90;
utils.createFigureA4(fig_params);

subplot(3, 1, 1);
hold on;
for idx = 1:numel(results)
    plot(results{idx}.t, results{idx}.signals.speedRpm, 'Color', utils.tab10Color(idx), 'LineWidth', 1.3);
end
xlabel('Time (s)');
ylabel('Speed (rpm)');
title('Speed Response under Different Loads');
legend({load_cases.name}, 'Location', 'best');
grid on;

subplot(3, 1, 2);
hold on;
for idx = 1:numel(results)
    plot(results{idx}.t, results{idx}.signals.current, 'Color', utils.tab10Color(idx), 'LineWidth', 1.3);
end
xlabel('Time (s)');
ylabel('Current (A)');
title('Current Response under Different Loads');
legend({load_cases.name}, 'Location', 'best');
grid on;

subplot(3, 1, 3);
hold on;
for idx = 1:numel(results)
    plot(results{idx}.t, results{idx}.signals.tauLoad, 'Color', utils.tab10Color(idx), 'LineWidth', 1.3);
end
xlabel('Time (s)');
ylabel('Torque (N m)');
title('Friction Load Torque under Different Loads');
legend({load_cases.name}, 'Location', 'best');
grid on;
