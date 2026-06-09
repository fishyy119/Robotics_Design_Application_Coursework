function G = gravity_vector(q, pi_dyn, g)
% 计算关节空间重力项。

G = [ ...
    g * pi_dyn(1) * cos(q(1)) + g * pi_dyn(3) * cos(q(1) + q(2)); ...
    g * pi_dyn(3) * cos(q(1) + q(2))];
end
