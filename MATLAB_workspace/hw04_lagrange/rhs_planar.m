function zdot = rhs_planar(~, z, m1, m2, l1, l2, g)

q1 = z(1);
u1 = z(2);
q2 = z(3);
u2 = z(4);

a1 = l1 / 2;
a2 = l2 / 2;

I1 = (1 / 12) * m1 * l1^2;
I2 = (1 / 12) * m2 * l2^2;
%% ===== 质量矩阵 M(q) =====
M11 = I1 + I2 + m1 * a1^2 + m2 * (l1^2 + a2^2 + 2 * l1 * a2 * cos(q2));
M12 = I2 + m2 * (a2^2 + l1 * a2 * cos(q2));
M21 = M12;
M22 = I2 + m2 * a2^2;

M = [M11, M12; ...
    M21, M22];
%% ===== 科里奥利/离心项 =====
h = -m2 * l1 * a2 * sin(q2);

C1 = h * (2 * u1 * u2 + u2^2);
C2 = -h * u1^2;
%% ===== 重力项（作用在质心） =====
G1 = (m1 * a1 + m2 * l1) * g * cos(q1) + m2 * a2 * g * cos(q1+q2);
G2 = m2 * a2 * g * cos(q1+q2);
%% 无控制
tau = [0; 0];

RHS = tau - [C1; C2] - [G1; G2];

qdd = M \ RHS;

zdot = [u1; qdd(1); u2; qdd(2)];

end
