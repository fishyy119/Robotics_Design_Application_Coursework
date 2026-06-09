%Variables initialization 
%Example Scheme A in Section 8.7 
%TextBook: Bruno Siciliano, et al. Robotics: Modelling, Planning and Control, Springer,2009

global a k_r1 k_r2 pi_m pi_l

% load manipulator dynamic parameters with load mass
  param
  load_m;
%  pi_m = pi_l; %pi_l:带载10kg与不带载
% gravity acceleration
  g = 9.81;
% friction matrix
  K_r = diag([k_r1 k_r2]);
  F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
  Tc = 0.001;

% controller gains
  K_p = 25*diag([1 1]);
  K_d = 5*diag([1 1]);


  B_h = diag([a(1)*pi_l(1) + pi_l(2) + a(2) + pi_l(4);
             a(2)*pi_l(3) + pi_l(4) + k_r2*k_r2*pi_l(5)]);

  D  = [zeros(2,2);eye(2,2)];
  H = [D';-K_p -K_d];
  Q = lyap2(H',eye(4,4)); %P=eye(4)
%X = lyap2(A,C) solves the special form of the Lyapunov matrix equation:
% A*X + X*A' = -C
  DQ = D'*Q;
  rho = 70;
  ep = 0.004;


% duration of simulation
  t_d = 4;

% trajectory generation
  GenTraject;

% sample time for plots
  Ts = Tc;
% for the initial value of the integral in the model of robot 
  q_i=[-1.4706; 2.9413]; 
  c=0
