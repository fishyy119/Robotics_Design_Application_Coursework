function x_dot = manipulator_rhs(x, tau, params)
% 计算真实带载两连杆系统的状态导数。

q = x(1:2);
dq = x(3:4);
pi_dyn = params.pi_l;
B = inertia_matrix(q(2), pi_dyn, params);
C = coriolis_vector(q, dq, pi_dyn, params);
G = gravity_vector(q, pi_dyn, params.g);
ddq = B \ (tau - params.F_v * dq - C - G);
x_dot = [dq; ddq];
end
