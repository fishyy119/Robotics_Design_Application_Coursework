%TRA_GENERATION  Generation of trajectory.
% clear;
% Tc = 0.001;  t_d = 4.0; a = [1, 1];
% initial joint configuration
  p_i = [0.2; 0];

% final joint configuration
  p_f = [1.8; 0];

% time base vector
  T = (0:Tc:t_d)';

% final time
  t_f = 3;

% trapezoidal velocity profile trajectory for end-effector 
  dp_c = 0.8;               % maximum velocity
  [T_1,px,~,~,err] = Trapezoid(p_i(1),p_f(1),dp_c,t_f,Tc);

% joint space trajectory of t_d sec duration
  n = size(T,1);
  m = size(T_1,1);
  q_d = zeros(n,2);
  dq_d = q_d;
  ddq_d = q_d;
  
  p_d = zeros(n, 2);
  p_d(1:m, :) = [px, zeros(size(px))];
  p_d(m+1:n,:) = ones(n-m,1)*p_f';
  
  for i = 1:n
      x = p_d(i, 1);    y = p_d(i, 2);
      [q1,q2] = inv_kinematics(x,y,a);
      q_d(i, :) = [q1, q2];
  end
  
  dq_d(2:n, :) = diff(q_d)/Tc;
  ddq_d(2:n, :) = diff(dq_d)/Tc;
  
%   clear T_1 q_1 dq_1 ddq_1 T_2 q_2 dq_2 ddq_2

