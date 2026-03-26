function zdot = rhs_planar_PD(t, z, m1, m2, l1, l2, g, qd_fun, qd_dot_fun, Kp, Kd)

q = z([1,3]);
qd  = qd_fun(t);
qd_dot = qd_dot_fun(t);
u = z([2,4]);

%% 刚体参数
a1 = l1/2; a2 = l2/2;
I1 = (1/12)*m1*l1^2; I2 = (1/12)*m2*l2^2;

%% 质量矩阵
M11 = I1 + I2 + m1*a1^2 + m2*(l1^2 + a2^2 + 2*l1*a2*cos(q(2)));
M12 = I2 + m2*(a2^2 + l1*a2*cos(q(2)));
M21 = M12;
M22 = I2 + m2*a2^2;
M = [M11 M12; M21 M22];

%% Coriolis/离心项
h = -m2*l1*a2*sin(q(2));
C = [h*(2*u(1)*u(2) + u(2)^2);
     -h*u(1)^2];

%% 重力
G = [(m1*a1 + m2*l1)*g*cos(q(1)) + m2*a2*g*cos(q(1)+q(2));
     m2*a2*g*cos(q(1)+q(2))];

%% PD 控制
tau = Kp*(qd - q) + Kd*(qd_dot - u);

% 记录力矩到全局变量
global tau_record t_record
t_record(end+1) = t;
tau_record(:,end+1) = tau;

%% 计算加速度
qdd = M \ (tau - C - G);

zdot = [u(1); qdd(1); u(2); qdd(2)];
end