function C = coriolis_vector(q, dq, pi_dyn, params)
% 计算离心与科氏力项。

s2 = sin(q(2));
h = -params.a(1) * pi_dyn(3) * s2;
C = [h * dq(2) * (2 * dq(1) + dq(2)); -h * dq(1)^2];
end
