function params = build_model(load_mass, load_mass_bar)
% 构造平面两连杆机械臂的标称、真实与上界模型参数。
% load_mass 为真实末端负载质量，load_mass_bar 为已知负载上界。

if nargin < 2
    load_mass_bar = load_mass;
end

a = [1; 1];
m_l1 = 50;
m_l2 = 50;
l_1 = 0.5;
l_2 = 0.5;
I_l1 = 10;
I_l2 = 10;
m_m1 = 5;
m_m2 = 5;
I_m1 = 0.01;
I_m2 = 0.01;
k_r1 = 100;
k_r2 = 100;

m_1 = m_l1 + m_m2;
m_2 = m_l2;
m1_lC1 = m_l1 * (l_1 - a(1));
m2_lC2 = m_l2 * (l_2 - a(2));
I_1 = I_l1 + m_l1 * (l_1 - a(1))^2 + I_m2;
I_2 = I_l2 + m_l2 * (l_2 - a(2))^2;

pi_m = [ ...
    a(1) * m_1 + m1_lC1 + a(1) * m_2; ...
    a(1) * m1_lC1 + I_1 + k_r1^2 * I_m1; ...
    a(2) * m_2 + m2_lC2; ...
    a(2) * m2_lC2 + I_2; ...
    I_m2];

% 末端集中负载建模为附着在连杆 2 坐标系原点的点质量。
% 对 8.7 节这套参数化，连杆 2 坐标系原点位于末端，因此点载荷
% 相对该坐标系的一阶矩和转动惯量附加项都为 0。
m_2_load = m_2 + load_mass;
m2_lC2_load = m2_lC2;
I_2_load = I_2;

pi_l = [ ...
    a(1) * m_1 + m1_lC1 + a(1) * m_2_load; ...
    a(1) * m1_lC1 + I_1 + k_r1^2 * I_m1; ...
    a(2) * m_2_load + m2_lC2_load; ...
    a(2) * m2_lC2_load + I_2_load; ...
    I_m2];

% pi_bar 不是控制器知道的真实参数，而是鲁棒设计中用于构造
% 惯性上界的保守包络参数，其负载上界由 load_mass_bar 给定。
m_2_bar = m_2 + load_mass_bar;
m2_lC2_bar = m2_lC2;
I_2_bar = I_2;

pi_bar = [ ...
    a(1) * m_1 + m1_lC1 + a(1) * m_2_bar; ...
    a(1) * m1_lC1 + I_1 + k_r1^2 * I_m1; ...
    a(2) * m_2_bar + m2_lC2_bar; ...
    a(2) * m2_lC2_bar + I_2_bar; ...
    I_m2];

K_r = diag([k_r1, k_r2]);

params.a = a;
params.g = 9.81;
params.k_r1 = k_r1;
params.k_r2 = k_r2;
params.pi_m = pi_m;
params.pi_l = pi_l;
params.pi_bar = pi_bar;
params.load_mass = load_mass;
params.load_mass_bar = load_mass_bar;
params.F_v = K_r * diag([0.01, 0.01]) * K_r;
end
