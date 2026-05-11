clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 电机参数
motor = getMotorParams('RE40 150W');

%% 负载与控制器
load_case.name = 'Speed-loop load';
load_case.J_load = 1.8e-3;
load_case.B = 1.8e-3;
load_case.tau_c = 0.020;
load_case.w_smooth = 1.0;
load_case.J_total = motor.J + load_case.J_load;
controller.mode = 'speed_loop';
controller.voltageLimit = motor.U;
controller.currentPI.Kp = 16.0;
controller.currentPI.Ki = 3.8e4;
controller.currentPI.limit = motor.U;
controller.speedPI.Kp = 0.26;
controller.speedPI.Ki = 3.2;
controller.speedPI.limit = 1.70;
command.omegaRefFcn = @(t) 55.0 * double(t >= 0.05) - 25.0 * double(t >= 1.00);

%% 数值积分
x0 = [0; 0; 0; 0; 0];
tspan = linspace(0, 2.0, 1601);
result = simulateMotorResponse(tspan, x0, motor, load_case, controller, command);

%% 绘图
fig_params.Name = 'Speed Loop Control';
fig_params.Width = 18;
fig_params.AspectRatio = 0.6;
utils.createFigureA4(fig_params);

subplot(2, 2, 1);
hold on;
h_measured = plot(result.t, result.signals.speedRpm, 'Color', utils.tab10Color(1, 0.65), 'LineWidth', 1.3);
h_ref = plot(result.t, result.signals.omegaRef * 30 / pi, 'k--', 'LineWidth', 1.1);
xlabel('Time (s)');
ylabel('Speed (rpm)');
title('Speed Loop Tracking');
legend([h_ref, h_measured], {'Reference', 'Measured'}, 'Location', 'best');
grid on;

subplot(2, 2, 2);
hold on;
h_measured = plot(result.t, result.signals.current, 'Color', utils.tab10Color(2, 0.65), 'LineWidth', 1.3);
h_ref = plot(result.t, result.signals.currentRef, 'k--', 'LineWidth', 1.1);
xlabel('Time (s)');
ylabel('Current (A)');
title('Inner-Loop Current Response');
legend([h_ref, h_measured], {'Current reference', 'Measured current'}, 'Location', 'best');
grid on;

subplot(2, 2, 3);
plot(result.t, result.signals.voltage, 'Color', utils.tab10Color(3), 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Applied Voltage');
grid on;

subplot(2, 2, 4);
hold on;
plot(result.t, result.signals.tauLoad, 'Color', utils.tab10Color(4), 'LineWidth', 1.3);
plot(result.t, result.signals.tauMotor, 'Color', utils.tab10Color(5), 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Torque (N m)');
title('Load Torque and Motor Torque');
legend({'Load torque', 'Motor torque'}, 'Location', 'best');
grid on;
