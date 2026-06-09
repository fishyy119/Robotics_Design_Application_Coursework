function [Tau, dot_pi_hat] = Adaptive_Controller(q, dq, ddq_r, dq_r, sig, pi_hat)
%ADAPTIVE_CONTROLLER  Utilize adaptive control law to obtain joint torque.
%                     [Tau, dm_l] = Adaptive_Controller(q, dq, ddq_r, dq_r, sig, m_l)
%                     returns joint torque and the rate of m_l.

global a k_r2 g pi_m F_v

c_1   = cos(q(1));
c_2   = cos(q(2));
s_2   = sin(q(2));
c_12  = cos(q(1) + q(2));
Y_r = zeros(2, 5);

% evaluates complete regressor 公式：8.99 Y
Y_r(1,1) = a(1)*ddq_r(1) + g*c_1;
Y_r(1,2) = ddq_r(1);
Y_r(1,4) = ddq_r(1) + ddq_r(2);
Y_r(1,3) = a(1)*c_2*(ddq_r(1) + Y_r(1,4)) - a(1)*s_2*(dq(2)*dq_r(1) +...
           dq_r(2)*(dq(1) + dq(2))) + a(2)*Y_r(1,4) + g*c_12;
% Y_r(1,5) = k_r2*ddq_r(2);
Y_r(1,5) = k_r2*ddq_r(2) + Y_r(1,2);
Y_r(2,3) = a(2)*Y_r(1,4) + a(1)*c_2*ddq_r(1) + a(1)*s_2*dq(1)*dq_r(1) +...
           g*c_12;
Y_r(2,4) = Y_r(1,4);
Y_r(2,5) = k_r2*Y_r(1,5);

% builds load regressor
Y_l = [a(1)*Y_r(1,1) + a(2)*Y_r(1,3);
       a(1)*Y_r(2,1) + a(2)*Y_r(2,3)];

dot_pi_hat = Y_l'*sig; % 公式：8.104 --- 更新参数向量估计

Tau =  Y_r*pi_m' + F_v*dq_r + Y_l*pi_hat;  % 公式：8.91

end
