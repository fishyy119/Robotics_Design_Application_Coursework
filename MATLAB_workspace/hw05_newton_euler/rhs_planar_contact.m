function zdot = rhs_planar_contact(~, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground)
% 本文件由 derive_planar_contact.m 自动生成。
% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。

%% 状态量
q1 = z(1); u1 = z(2);
q2 = z(3); u2 = z(4);

%% 接触量
tip_y = l1*sin(q1) + l2*cos(q1)*sin(q2) + l2*cos(q2)*sin(q1);
tip_y_dot = l1*u1*cos(q1) - l2*u1*sin(q1)*sin(q2) - l2*u2*sin(q1)*sin(q2) + l2*u1*cos(q1)*cos(q2) + l2*u2*cos(q1)*cos(q2);
penetration = max(ground_y - tip_y, 0);
approach_speed = max(-tip_y_dot, 0);

%% 接触力
Fx = 0;
Fy = 0;
if penetration > 0
    Fy = k_ground * penetration + c_ground * approach_speed;
end

%% 外加关节力矩
T1 = 0;
T2 = 0;

%% 动力学方程
M11 = - J1 - J2 - a1^2*m1*cos(q1)^2 - l1^2*m2*cos(q1)^2 - a1^2*m1*sin(q1)^2 - l1^2*m2*sin(q1)^2 - a2^2*m2*cos(q1)^2*sin(q2)^2 - a2^2*m2*cos(q2)^2*sin(q1)^2 - a2^2*m2*sin(q1)^2*sin(q2)^2 - a2^2*m2*cos(q1)^2*cos(q2)^2 - 2*a2*l1*m2*cos(q1)^2*cos(q2) - 2*a2*l1*m2*cos(q2)*sin(q1)^2;
M12 = - J2 - a2^2*m2*cos(q1)^2*sin(q2)^2 - a2^2*m2*cos(q2)^2*sin(q1)^2 - a2^2*m2*sin(q1)^2*sin(q2)^2 - a2^2*m2*cos(q1)^2*cos(q2)^2 - a2*l1*m2*cos(q1)^2*cos(q2) - a2*l1*m2*cos(q2)*sin(q1)^2;
M21 = - J2 - a2^2*m2*cos(q1)^2*sin(q2)^2 - a2^2*m2*cos(q2)^2*sin(q1)^2 - a2^2*m2*sin(q1)^2*sin(q2)^2 - a2^2*m2*cos(q1)^2*cos(q2)^2 - a2*l1*m2*cos(q1)^2*cos(q2) - a2*l1*m2*cos(q2)*sin(q1)^2;
M22 = - J2 - a2^2*m2*cos(q1)^2*sin(q2)^2 - a2^2*m2*cos(q2)^2*sin(q1)^2 - a2^2*m2*sin(q1)^2*sin(q2)^2 - a2^2*m2*cos(q1)^2*cos(q2)^2;

RHS1 = Fx*l1*sin(q1) - Fy*l1*cos(q1) - T1 + a1*g*m1*cos(q1) + g*l1*m2*cos(q1) - Fy*l2*cos(q1)*cos(q2) + Fx*l2*cos(q1)*sin(q2) + Fx*l2*cos(q2)*sin(q1) + Fy*l2*sin(q1)*sin(q2) + a2*g*m2*cos(q1)*cos(q2) - a2*g*m2*sin(q1)*sin(q2) - a2*l1*m2*u2^2*cos(q1)^2*sin(q2) - a2*l1*m2*u2^2*sin(q1)^2*sin(q2) - 2*a2*l1*m2*u1*u2*cos(q1)^2*sin(q2) - 2*a2*l1*m2*u1*u2*sin(q1)^2*sin(q2);
RHS2 = Fx*l2*cos(q1)*sin(q2) - Fy*l2*cos(q1)*cos(q2) - T2 + Fx*l2*cos(q2)*sin(q1) + Fy*l2*sin(q1)*sin(q2) + a2*g*m2*cos(q1)*cos(q2) - a2*g*m2*sin(q1)*sin(q2) + a2*l1*m2*u1^2*cos(q1)^2*sin(q2) + a2*l1*m2*u1^2*sin(q1)^2*sin(q2);

MM = [M11, M12; M21, M22];
qdd = MM \ [RHS1; RHS2];

zdot = [u1; qdd(1); u2; qdd(2)];
