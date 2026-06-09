% ============= Robust Control ================= %
%MAIN_Programs.m
% This script runs a simulation of two dofs manipulator with a load mass
% implementing a robust control law. Refer to pp.332 in Bruno's book.

clear; close all; clc

addpath ../Simulation_Lib

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

% B_hat
B_h = diag([a(1)*pi_l(1) + pi_l(2) + a(2) + pi_l(4);
             a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5)]);

D  = [zeros(2,2);eye(2,2)]; % 公式：8.71
H = [D';-K_p -K_d]; % H_tilt(公式：8.78)
Q = lyap2(H',eye(4,4)); % 见公式： 8.81（pp.335）
%X = lyap2(A,C) solves the special form of the Lyapunov matrix equation:
% A*X + X*A' = -C
DQ = D'*Q; % 公式：8.82
rho = 70; % 公式：8.83
ep = 0.004; % 公式：8.89 中的 epsilon

% duration of simulation
t_d = 4.0;
  
% trajectory generation
tra_generation_ee;

% simulation environment config
mod = MatlabVrepTemplate(Tc * 4000);
mod = mod.init();
disp('Simulation Start!!!');
mod.go();

Vel_Limit = 9999999;
TAU = zeros(n, 2);
dir = zeros(2, 1);
Q_rel_last = zeros(1, 2);
Q_rel_r = zeros(n, 2);

for i = 1:n
    % get real joint position
    Q_rel = mod.read_joint();
    if i == 1
        Q_rel = q_d(1, :);
    end
    Q_rel_r(i, :) = Q_rel;
    
    % obtain real joint velocity
    if i == 1
        dQ_rel = [0.0, 0.0];
    else
        dQ_rel = (Q_rel - Q_rel_last) / Tc;
    end
    Q_rel_last = Q_rel;
    
    % calculate the errors of joint pos and vel
    q_err = q_d(i, :) - Q_rel;
    dq_err = dq_d(i, :) - dQ_rel;
    Xi = [q_err dq_err]';% 公式：8.69
    
    % get robust component w
    z = DQ * Xi; % 2*1 matrix
    if norm(z) >= ep % 公式：8.89
        w = rho / norm(z) * z;
    else
        w = rho / ep * z;
    end
    
    % get the inertia torque
    y = w + ddq_d(i, :)' + K_p * q_err' + K_d * dq_err';
    Inertia_Tau = B_h * y;
    
    % get the friction and gravity compensation
    Tau_FGC = F_v * dQ_rel' + [g*pi_m(1)*cos(Q_rel(1)) + g*pi_m(3)*cos(Q_rel(1)+Q_rel(2)); g*pi_m(3)*cos(Q_rel(1)+Q_rel(2))];
    
    % integrate Inertia_Tau and Tau_FGC to get total torque
    Tau = Inertia_Tau + Tau_FGC;
    TAU(i, :) = Tau';
    
    for j = 1:2
        if Tau(j) > 0
            dir(j) = Vel_Limit;
        else
            dir(j) = - Vel_Limit;
        end        
    end

    mod.enable_joint_torque(); % necessary setting : switch from position control to torque control

    for j = 1:2
        mod.Main.simxSetJointTargetVelocity(mod.ClientID, mod.Joint_Handle(j), dir(j),mod.Main.simx_opmode_oneshot);
        mod.Main.simxSetJointForce(mod.ClientID, mod.Joint_Handle(j), abs(Tau(j)), mod.Main.simx_opmode_oneshot);        
    end
    
    mod.trigger();
end

mod.stop();