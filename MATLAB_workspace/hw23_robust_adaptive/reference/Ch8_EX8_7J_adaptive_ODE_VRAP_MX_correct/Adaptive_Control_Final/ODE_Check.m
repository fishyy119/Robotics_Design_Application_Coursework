% Adaptive Control——ODE Check Uses ODE113 and some other numerical
% integration methods(such as rk4) to check whether the dynamics derivation
% is right.
clear; close all; clc

global a k_r1 k_r2 pi_m pi_l F_v g Lam K_ml K_d

% load manipulator dynamic parameters without load mass
params;
% load mass parameters
load_m;

% gravity acceleration
g = 9.81;

% transmission ratio and friction matrix
K_r = diag([k_r1 k_r2]);
F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
Tc = 0.001;

% controller gains(refer to Scheme.J pp.349)
Lam = 5*diag([1 1]);
K_d = 750*diag([1 1]);
K_ml = 0.01;

% initial estimate of load mass
m_l0 = 0;

% duration of simulation
t_d = 4.0;
  
% trajectory generation
tra_generation_ee; % obtain q_d dq_d ddq_d 

% tspan, options and initial state vector for ode113
tspan = [0, Tc];
options = odeset('abstol', 1e-9, 'reltol', 1e-9);
% state vector: z = [q1; dq1; q2; dq2; m_l]
z0 = [q_d(1, 1) dq_d(1, 1) q_d(1, 2) dq_d(1, 2) m_l0]';

z_record = zeros(n, 5);

% Dynamics function handle
DynFunc_mx = @(z, q_ref, dq_ref, ddq_ref)DynFunc([], z, q_ref, dq_ref, ddq_ref);

NumIntgMethod = 'euler'; % you can choose different methods: {'ode113', 'ode45', 'rk4', 'euler', 'midpoint'}

for i = 1:n
    
    q_ref = q_d(i, :);
    dq_ref = dq_d(i, :);
    ddq_ref = ddq_d(i, :);
    
    switch NumIntgMethod
        case 'ode113'
            [~, z_temp] = ode113(@DynFunc, tspan, z0, options, q_ref, dq_ref, ddq_ref);
            z0_temp = z_temp(end, :)';
        case 'ode45'
            [~, z_temp] = ode45(@DynFunc, tspan,  z0, options, q_ref, dq_ref, ddq_ref);
            z0_temp = z_temp(end, :)';
        case 'rk4'
            k1 = feval(DynFunc_mx,                z0, q_ref, dq_ref, ddq_ref);
            k2 = feval(DynFunc_mx, z0 + 0.5 * Tc * k1, q_ref, dq_ref, ddq_ref);
            k3 = feval(DynFunc_mx, z0 + 0.5 * Tc * k2, q_ref, dq_ref, ddq_ref);
            k4 = feval(DynFunc_mx,       z0 + Tc * k3, q_ref, dq_ref, ddq_ref);
            z0_temp = z0 + Tc*(1/6)*(k1 + 2*k2 + 2*k3 + k4);
        case 'euler'
            k1 = feval(DynFunc_mx,                 z0, q_ref, dq_ref, ddq_ref);
            z0_temp = z0 + Tc * k1;
        case 'midpoint'
            k1 = feval(DynFunc_mx,                 z0, q_ref, dq_ref, ddq_ref);
            k2 = feval(DynFunc_mx, z0 + 0.5 * Tc * k1, q_ref, dq_ref, ddq_ref);
            z0_temp = z0 + Tc * k2;
        otherwise
            error('Invalid method!');
    end
    
    z0 = z0_temp;
    z_record(i, :) = z0';
end

% joint angle comparison
figure(1); clf;
subplot(1, 2, 1); grid on; hold on;
plot(T, q_d(:, 1), 'LineWidth', 2, 'DisplayName', 'q1-ref');
hold on
plot(T, z_record(:, 1), 'r--', 'LineWidth', 2, 'DisplayName', 'q1-act');
xlabel('[s]');  ylabel('[rad]');
title('关节1角度');
legend
subplot(1, 2, 2); grid on; hold on;
plot(T, q_d(:, 2), 'LineWidth', 2, 'DisplayName', 'q2-ref');
hold on
plot(T, z_record(:, 3), 'r--', 'LineWidth', 2, 'DisplayName', 'q2-act');
xlabel('[s]');  ylabel('[rad]');
title('关节2角度');
legend

figure(2); clf;
plot(T, q_d(:, 1) - z_record(:, 1), 'b-', 'DisplayName', 'q1');
hold on
plot(T, q_d(:, 2) - z_record(:, 3), 'r--', 'DisplayName', 'q2');
grid on
xlabel('[s]');  ylabel('[rad]');
title('关节角度误差');
legend

figure(3); clf;
x_des_mx = a(1) .* cos(q_d(:, 1)) + a(2) .* cos(q_d(:, 1) + q_d(:, 2));
y_des_mx = a(1) .* sin(q_d(:, 1)) + a(2) .* sin(q_d(:, 1) + q_d(:, 2));
x_rel_mx = a(1) .* cos(z_record(:, 1)) + a(2) .* cos(z_record(:, 1) + z_record(:, 3));
y_rel_mx = a(1) .* sin(z_record(:, 1)) + a(2) .* sin(z_record(:, 1) + z_record(:, 3));
err_pos_ee = [x_des_mx, y_des_mx] - [x_rel_mx, y_rel_mx];
Nt = length(err_pos_ee);     Norm_p_mx = zeros(Nt, 1);
for i = 1:Nt
    Norm_p_mx(i) = norm(err_pos_ee(i, :));
end
plot(T, Norm_p_mx, 'b-');
grid on
xlabel('[s]');  ylabel('[m]');
title('机械臂末端位置误差2范数');