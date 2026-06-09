%Variables initialization 
%Example Scheme J (adaptive control) in Section 8.7 
%TextBook: Bruno Siciliano, et al. Robotics: Modelling, Planning and Control, Springer,2009

clc;
clear all;
global a k_r1 k_r2 pi_m pi_l F_v g

c=0;
% load manipulator dynamic parameters with load mass
  param
  load_m;
% gravity acceleration
  g = 9.81;
% friction matrix
  K_r = diag([k_r1 k_r2]);
  F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
  Tc = 0.001;
 
% controller gains
  Lam = 5*diag([1 1]);
  K_d = 750*diag([1 1]);
  K_ml = 0.01;  %refer to K_pai in textbook

% initial estimate of load mass
  m_l0 = 0;


% duration of simulation
  t_d = 4;

% trajectory generation
  GenTraject;

% sample time for plots
  Ts = Tc;
% for the initial value of the integral in the model of robot 
  q_i=[-1.4706; 2.9413];