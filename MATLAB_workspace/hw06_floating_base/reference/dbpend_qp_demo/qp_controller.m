function [T] = qp_controller(z, m1, m2, I1, I2, l, a, g)
%QP_CONTROLLER 此处显示有关此函数的摘要
%   此处显示详细说明
q1 = z(1);                          
u1 = z(2);                          
q2 = z(3);                        
u2 = z(4);  

% Reference
q1_ref = pi;
u1_ref = 0;
q2_ref = 0;
u2_ref = 0;                     

% [M][alpha] = RHS + T
M11 = I1 + I2 + a^2*m1 + a^2*m2 + l^2*m2 + 2*a*l*m2*cos(q2);%-I1-a^2*m1-m2*l^2;
M12 = m2*a^2 + l*m2*cos(q2)*a + I2;%-cos(-q1+q2)*l*a*m2;
M21 = m2*a^2 + l*m2*cos(q2)*a + I2;%-cos(-q1+q2)*l*a*m2;
M22 = m2*a^2 + I2;%-m2*a^2-I2;
RHS1 = a*l*m2*sin(q2)*u2^2 + 2*a*l*m2*u1*sin(q2)*u2 - a*g*m2*sin(q1 + q2) - a*g*m1*sin(q1) - g*l*m2*sin(q1);%m2*g*l*sin(q1)-m2*l*u2^2*a*sin(-q1+q2)+a*sin(q1)*m1*g-T1+T2;
RHS2 = - a*l*m2*sin(q2)*u1^2 - a*g*m2*sin(q1 + q2);%a*sin(q2)*m2*g-T2+m2*a*u1^2*l*sin(-q1+q2);
M    = [M11 M12; M21 M22];
RHS  = [RHS1;RHS2];

% Simple PD controller
kp = 100;
kd = 30;
% T1 = kp*(q1_ref-q1)+kd*(u1_ref-u1); %zero torques for unactuated system
% T2 = kp*(q2_ref-q2)+kd*(u2_ref-u2); 
% T1 = max(T1,-100); T1 = min(T1,100);
% T2 = max(T2,-100); T2 = min(T2,100);

% QP controller
% min 0.5*|Aw*x-bw|^2 ==> min x'*Aw'*Aw*x + (-Aw'*bw)'*x
%                         min x'*   H  *x +     f'    *x
%    Aw*x - bw = [w1*(x1-b1)] = [w1*[1 0 0 0]*x1 - w1*b1] = [w1*[1 0 0 0]] * [x1] - [w1*b1] 
%                [w2*(x2-b2)]   [w2*[0 1 0 0]*x2 - w2*b2]   [w2*[0 1 0 0]]   [x2]   [w2*b2]  
%                [w3*(x3-b3)]   [w3*[0 0 1 0]*x3 - w3*b3]   [w3*[0 0 1 0]]   [x3]   [w3*b3] 
%                [w4*(x4-b4)]   [w4*[0 0 0 1]*x4 - w4*b4]   [w4*[0 0 0 1]]   [x4]   [w4*b4] 
% s.t. Aeq *x  = beq
%      Aiep*x <= bieq
% x = [ud1,ud2,T1,T2]'

% Cost function
% Reference ud
u1dot_ref = kp*(q1_ref-q1)+kd*(u1_ref-u1);
u2dot_ref = kp*(q2_ref-q2)+kd*(u2_ref-u2); 
w = [10,10,0,0];    % Weights
Aw = diag([1,1,1,1].*w);
bw = [u1dot_ref;u2dot_ref;0;0].*w';
H =  Aw'*Aw;
f = -Aw'*bw;

% Eq. constraints
% M*ud - T = RHS ==> [M11, M12, -1, 0 ] * [ud1,ud2,T1,T2]' = [RHS1]
%                    [M21, M22,  0, -1]                      [RHS2]
Aeq = [M,[-1 0;0 -1]];
beq = RHS;

% Ineq. constraints
% T1 <=  T1_max
%-T1 <= -T1_min
% T2 <=  T2_max
%-T2 <= -T2_min
% ud1 <=  ud1_max
%-ud1 <= -ud1_min
% ud2 <=  ud2_max
%-ud2 <= -ud2_min
Aieq = [0  0  1  0;
        0  0 -1  0;
        0  0  0  1;
        0  0  0 -1;
        1  0  0  0;
       -1  0  0  0;
        0  1  0  0;
        0 -1  0  0;];
Bieq = [100;100;100;100;50;50;50;50];

options = optimoptions('quadprog','Display','off');
x = quadprog(H,f,Aieq,Bieq,Aeq,beq,[],[],[],options);
T1 = x(3); T2 = x(4);

T = [T1;T2];

end

