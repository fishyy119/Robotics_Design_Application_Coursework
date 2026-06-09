% ============= Adaptive Control ================= %
%MAIN_Programs.m
% This script runs a simulation of two dofs manipulator with a load mass
% implementing an adaptive control law. Refer to pp.339 in Bruno's book.

clear; close all; clc

addpath ..\Simulation_Lib

global a k_r1 k_r2 pi_m pi_l F_v g

% load manipulator dynamic parameters without load mass
params;
load_m;

% gravity acceleration
g = 9.81;

% friction matrix
K_r = diag([k_r1 k_r2]);
F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
Tc = 0.001;

% controller gains
Lam = 5*diag([1 1]); % Lamda matrix 见公式：8.92
K_d = 750*diag([1 1]);
K_ml = 0.01;

% initial estimate of load mass
m_l0 = 0;
m_l = 0;
dm_l = 0;

% duration of simulation
t_d = 4.0;
  
% trajectory generation
tra_generation_ee;

% simulation environment configs
mod = MatlabVrepTemplate(Tc * 1000);
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

    dq_r = Lam * q_err' + dq_d(i, :)'; % 公式：8.92
    ddq_r = Lam * dq_err' + ddq_d(i, :)';
    
    sig = dq_r - dQ_rel'; % 公式：8.93
    
    % update parameter vector estimation
    if i == 1
        m_l = m_l0;
    else
        m_l = m_l + dm_l / K_ml * Tc;
    end
    
    % utilize the adaptive controller to calculate the joint torque and
    % dot_Pi
    [Tau, dm_l] = Adaptive_Controller(Q_rel, dQ_rel, ddq_r, dq_r, sig, m_l);
    TAU(i, :) = (Tau + K_d * sig)';    
    
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