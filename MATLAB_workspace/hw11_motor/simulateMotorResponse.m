function result = simulateMotorResponse(tspan, x0, motor, load_case, controller, command)
% 对给定负载与控制器进行电机数值仿真。

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', tspan(2) - tspan(1));
[t, x] = ode45(@(t, x) motorRhs(t, x, motor, load_case, controller, command), tspan, x0, options);

signals = evaluateSignals(t, x, motor, load_case, controller, command);

result.t = t;
result.x = x;
result.signals = signals;
result.load = load_case;
result.controller = controller;
end

function dx = motorRhs(t, x, motor, load_case, controller, command)
% 计算直流电机在不同控制模式下的状态导数。

omega = x(2);
current = x(3);

control = evaluateController(t, x, motor, controller, command);
tau_load = computeLoadTorque(omega, load_case);

dtheta = omega;
domega = (motor.Kt * current - tau_load) / load_case.J_total;
dcurrent = (control.voltage - motor.R * current - motor.Ke * omega) / motor.L;

dx = [dtheta; domega; dcurrent; control.integratorDot];
end

function control = evaluateController(t, x, motor, controller, command)
% 根据控制模式计算电压指令和积分器导数。

theta = x(1);
omega = x(2);
current = x(3);

control.thetaRef = 0;
control.omegaRef = 0;
control.currentRef = 0;
control.voltage = 0;
control.integratorDot = [];

switch controller.mode
    case 'open_loop'
        control.voltage = clamp(controller.voltageFcn(t, x), controller.voltageLimit);

    case 'current_loop'
        control.currentRef = command.currentRefFcn(t);
        xi_i = x(4);
        [control.voltage, control.integratorDot] = saturatingPI( ...
            control.currentRef - current, xi_i, controller.currentPI.Kp, ...
            controller.currentPI.Ki, controller.currentPI.limit);

    case 'speed_loop'
        control.omegaRef = command.omegaRefFcn(t);
        xi_i = x(4);
        xi_w = x(5);
        [control.currentRef, dxi_w] = saturatingPI( ...
            control.omegaRef - omega, xi_w, controller.speedPI.Kp, ...
            controller.speedPI.Ki, controller.speedPI.limit);
        [control.voltage, dxi_i] = saturatingPI( ...
            control.currentRef - current, xi_i, controller.currentPI.Kp, ...
            controller.currentPI.Ki, controller.currentPI.limit);
        control.integratorDot = [dxi_i; dxi_w];

    case 'position_loop'
        control.thetaRef = command.thetaRefFcn(t);
        xi_i = x(4);
        xi_w = x(5);
        control.omegaRef = clamp(controller.positionP.Kp * (control.thetaRef - theta), controller.positionP.limit);
        [control.currentRef, dxi_w] = saturatingPI( ...
            control.omegaRef - omega, xi_w, controller.speedPI.Kp, ...
            controller.speedPI.Ki, controller.speedPI.limit);
        [control.voltage, dxi_i] = saturatingPI( ...
            control.currentRef - current, xi_i, controller.currentPI.Kp, ...
            controller.currentPI.Ki, controller.currentPI.limit);
        control.integratorDot = [dxi_i; dxi_w];

    otherwise
        error('Unsupported control mode: %s', controller.mode);
end

control.voltage = clamp(control.voltage, motor.U);
end

function [u, dxi] = saturatingPI(error_value, xi, Kp, Ki, limit)
% 计算带饱和与简单抗积分饱和的 PI 输出。

u_unsat = Kp * error_value + Ki * xi;
u = clamp(u_unsat, limit);

dxi = error_value;
if abs(u_unsat) > limit && sign(u_unsat) == sign(error_value)
    dxi = 0;
end
end

function tau_load = computeLoadTorque(omega, load_case)
% 计算负载力矩与平滑摩擦力矩。

direction = tanh(omega / load_case.w_smooth);
tau_load = load_case.B * omega + load_case.tau_c * direction;
end

function signals = evaluateSignals(t, x, motor, load_case, controller, command)
% 后处理计算参考量、控制量和力矩。

n = numel(t);
signals.theta = x(:, 1);
signals.omega = x(:, 2);
signals.current = x(:, 3);
signals.speedRpm = x(:, 2) * 30 / pi;
signals.voltage = zeros(n, 1);
signals.thetaRef = zeros(n, 1);
signals.omegaRef = zeros(n, 1);
signals.currentRef = zeros(n, 1);
signals.tauLoad = zeros(n, 1);
signals.tauMotor = motor.Kt * x(:, 3);

for idx = 1:n
    control = evaluateController(t(idx), x(idx, :).', motor, controller, command);
    signals.voltage(idx) = control.voltage;
    signals.thetaRef(idx) = control.thetaRef;
    signals.omegaRef(idx) = control.omegaRef;
    signals.currentRef(idx) = control.currentRef;
    signals.tauLoad(idx) = computeLoadTorque(x(idx, 2), load_case);
end
end

function y = clamp(u, limit)
% 对称限幅。

y = min(max(u, -limit), limit);
end
