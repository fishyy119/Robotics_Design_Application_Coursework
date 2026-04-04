function [KE, PEg, PEc, Pdiss_damping, Pdiss_friction, Fy, penetration, tip_y, tip_y_dot] = energy_planar_contact(~, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground, mu_ground, v_smooth)
% 本文件由 derive_planar_contact.m 自动生成。
% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。

% 状态量
q1 = z(1); u1 = z(2);
q2 = z(3); u2 = z(4);

% 接触量
tip_x_dot = - l1*u1*sin(q1) - l2*u1*cos(q1)*sin(q2) - l2*u1*cos(q2)*sin(q1) - l2*u2*cos(q1)*sin(q2) - l2*u2*cos(q2)*sin(q1);
tip_y = l1*sin(q1) + l2*cos(q1)*sin(q2) + l2*cos(q2)*sin(q1);
tip_y_dot = l1*u1*cos(q1) - l2*u1*sin(q1)*sin(q2) - l2*u2*sin(q1)*sin(q2) + l2*u1*cos(q1)*cos(q2) + l2*u2*cos(q1)*cos(q2);
penetration = max(ground_y - tip_y, 0);
approach_speed = max(-tip_y_dot, 0);
Fy = 0;
Pdiss_damping = 0;
Pdiss_friction = 0;
if penetration > 0
    Fy = k_ground * penetration + c_ground * approach_speed;
    Fx = -mu_ground * Fy * tanh(tip_x_dot / max(v_smooth, 1e-6));
    Pdiss_damping = c_ground * approach_speed^2;
    Pdiss_friction = -Fx * tip_x_dot;
end
PEc = 0.5 * k_ground * penetration^2;

% 能量表达式
KE = (J1*u1^2)/2 + (J2*u1^2)/2 + (J2*u2^2)/2 + J2*u1*u2 + (a1^2*m1*u1^2*cos(q1)^2)/2 + (l1^2*m2*u1^2*cos(q1)^2)/2 + (a1^2*m1*u1^2*sin(q1)^2)/2 + (l1^2*m2*u1^2*sin(q1)^2)/2 + (a2^2*m2*u1^2*cos(q1)^2*cos(q2)^2)/2 + (a2^2*m2*u2^2*cos(q1)^2*cos(q2)^2)/2 + (a2^2*m2*u1^2*cos(q1)^2*sin(q2)^2)/2 + (a2^2*m2*u1^2*cos(q2)^2*sin(q1)^2)/2 + (a2^2*m2*u2^2*cos(q1)^2*sin(q2)^2)/2 + (a2^2*m2*u2^2*cos(q2)^2*sin(q1)^2)/2 + (a2^2*m2*u1^2*sin(q1)^2*sin(q2)^2)/2 + (a2^2*m2*u2^2*sin(q1)^2*sin(q2)^2)/2 + a2*l1*m2*u1^2*cos(q1)^2*cos(q2) + a2*l1*m2*u1^2*cos(q2)*sin(q1)^2 + a2^2*m2*u1*u2*cos(q1)^2*cos(q2)^2 + a2^2*m2*u1*u2*cos(q1)^2*sin(q2)^2 + a2^2*m2*u1*u2*cos(q2)^2*sin(q1)^2 + a2^2*m2*u1*u2*sin(q1)^2*sin(q2)^2 + a2*l1*m2*u1*u2*cos(q1)^2*cos(q2) + a2*l1*m2*u1*u2*cos(q2)*sin(q1)^2;
PEg = g*l1*m1 + g*l1*m2 + g*l2*m1 + g*l2*m2 + a1*g*m1*sin(q1) + g*l1*m2*sin(q1) + a2*g*m2*cos(q1)*sin(q2) + a2*g*m2*cos(q2)*sin(q1);
