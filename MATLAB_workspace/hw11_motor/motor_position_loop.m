clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 电机参数
motor = getMotorParams('RE40 150W');

%% 负载与控制器
load_case.name = 'Position-loop load';
load_case.J_load = 1.6e-3;
load_case.B = 1.5e-3;
load_case.tau_c = 0.018;
load_case.w_smooth = 1.0;
load_case.J_total = motor.J + load_case.J_load;
controller.mode = 'position_loop';
controller.voltageLimit = motor.U;
controller.currentPI.Kp = 16.0;
controller.currentPI.Ki = 3.8e4;
controller.currentPI.limit = motor.U;
controller.speedPI.Kp = 0.24;
controller.speedPI.Ki = 2.8;
controller.speedPI.limit = 1.60;
controller.positionP.Kp = 18.0;
controller.positionP.limit = 90.0;
command.thetaRefFcn = @(t) 0.50 * pi * sin(2 * pi * 0.25 * t);

%% 数值积分
x0 = [0; 0; 0; 0; 0];
tspan = linspace(0, 8.00, 4001);
result = simulateMotorResponse(tspan, x0, motor, load_case, controller, command);

%% 绘图
fig_params.Name = 'Position Loop Control';
fig_params.Width = 18;
fig_params.AspectRatio = 0.85;
utils.createFigureA4(fig_params);

subplot(3, 2, 1);
hold on;
h_measured = plot(result.t, result.signals.theta, 'Color', utils.tab10Color(1, 0.65), 'LineWidth', 1.3);
h_ref = plot(result.t, result.signals.thetaRef, 'k--', 'LineWidth', 1.1);
xlabel('Time (s)');
ylabel('Position (rad)');
title('Position Loop Tracking');
legend([h_ref, h_measured], {'Reference', 'Measured'}, 'Location', 'best');
grid on;

subplot(3, 2, 2);
hold on;
h_measured = plot(result.t, result.signals.speedRpm, 'Color', utils.tab10Color(2, 0.65), 'LineWidth', 1.3);
h_ref = plot(result.t, result.signals.omegaRef * 30 / pi, 'k--', 'LineWidth', 1.1);
xlabel('Time (s)');
ylabel('Speed (rpm)');
title('Generated Speed Command and Measured Speed');
legend([h_ref, h_measured], {'Speed reference', 'Measured speed'}, 'Location', 'best');
grid on;

subplot(3, 2, 3);
plot(result.t, result.signals.current, 'Color', utils.tab10Color(3), 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Current (A)');
title('Inner Current Response');
grid on;

subplot(3, 2, 4);
plot(result.t, result.signals.voltage, 'Color', utils.tab10Color(4), 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Applied Voltage');
grid on;

subplot(3, 2, [5, 6]);
hold on;
plot(result.t, result.signals.tauLoad, 'Color', utils.tab10Color(5), 'LineWidth', 1.3);
plot(result.t, result.signals.tauMotor, 'Color', utils.tab10Color(6), 'LineWidth', 1.3);
xlabel('Time (s)');
ylabel('Torque (N m)');
title('Load Torque and Motor Torque');
legend({'Load torque', 'Motor torque'}, 'Location', 'best');
grid on;
