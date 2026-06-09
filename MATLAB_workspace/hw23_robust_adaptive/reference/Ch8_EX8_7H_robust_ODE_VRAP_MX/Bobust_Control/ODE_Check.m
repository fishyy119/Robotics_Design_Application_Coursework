% Robust Control ODE check
clear; close all; clc

params;
load_m;
pi_m = pi_l;

% gravity acceleration
g = 9.81;

% friction matrix
K_r = diag([k_r1 k_r2]);
F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
Tc = 0.001;

% controller gains
K_d = 5*diag([1 1]);
K_p = 25*diag([1 1]);

B_h = diag([a(1)*pi_l(1) + pi_l(2) + a(2) + pi_l(4);
             a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5)]);

D  = [zeros(2,2);eye(2,2)];
H = [D';-K_p -K_d];
Q = lyap2(H',eye(4,4));
%X = lyap2(A,C) solves the special form of the Lyapunov matrix equation:
% A*X + X*A' = -C
DQ = D'*Q;
rho = 70;
ep = 0.004;

% duration of simulation
t_d = 4.0; % 2.5
  
% trajectory generation
tra_generation_ee;

% tspan, options and initial state vector for ode113
tspan = [0, Tc];
options = odeset('abstol',1e-9,'reltol',1e-9);
z0 = [q_d(1, 1) dq_d(1, 1) q_d(1, 2) dq_d(1, 2)]';
z_total = zeros(n, 4);

z_ode_mx = [];

for i = 1:n
    t = (i - 1) * Tc;
    q_ref = q_d(i, :);
    dq_ref = dq_d(i, :);
    ddq_ref = ddq_d(i, :);
    [t_temp, z_temp] = ode45(@DynFunc, tspan, z0, options, q_ref, dq_ref, ddq_ref, g, pi_m, pi_l, F_v, K_p, K_d, B_h, DQ, rho, ep, a, k_r2);
    
    z_ode_mx = [z_ode_mx; z_temp];
    
    z0 = z_temp(end,:);
    z_total(i, :) = z0;
end

% joint angle comparison
figure(1); clf;
subplot(1, 2, 1); grid on; hold on;
plot(T, q_d(:, 1), 'LineWidth', 2, 'DisplayName', 'q1-ref');
hold on
plot(T, z_total(:, 1), 'LineWidth', 2, 'DisplayName', 'q1-act');
xlabel('[s]');  ylabel('[rad]');
title('关节1角度');
legend
subplot(1, 2, 2); grid on; hold on;
plot(T, q_d(:, 2), 'LineWidth', 2, 'DisplayName', 'q2-ref');
hold on
plot(T, z_total(:, 3), 'LineWidth', 2, 'DisplayName', 'q2-act');
xlabel('[s]');  ylabel('[rad]');
title('关节2角度');
legend

figure(2); clf;
plot(T, q_d(:, 1) - z_total(:, 1), 'b-', 'DisplayName', 'q1');
hold on
grid on
plot(T, q_d(:, 2) - z_total(:, 3), 'r--', 'DisplayName', 'q2');
xlabel('[s]');  ylabel('[rad]');
title('关节角度误差');
legend