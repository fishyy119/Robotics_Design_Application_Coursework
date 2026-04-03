function [S, Sd, Sdd] = TSpline(p0, v0, t0, p1, t1, p2, v2, t2, t)
% 计算过三点约束的分段三次样条及其一阶、二阶导数。

%% 构造线性方程组
R = [1, t0, t0 * t0, t0 * t0 * t0, 0, 0, 0, 0; ...
    0, 1, 2 * t0, 3 * t0 * t0, 0, 0, 0, 0; ...
    1, t1, t1 * t1, t1 * t1 * t1, 0, 0, 0, 0; ...
    0, 0, 0, 0, 1, t2, t2 * t2, t2 * t2 * t2; ...
    0, 0, 0, 0, 0, 1, 2 * t2, 3 * t2 * t2; ...
    0, 0, 0, 0, 1, t1, t1 * t1, t1 * t1 * t1; ...
    0, 1, 2 * t1, 3 * t1 * t1, 0, -1, -2 * t1, -3 * t1 * t1; ...
    0, 0, 2, 6 * t1, 0, 0, -2, -6 * t1];

rhs = [p0, v0, p1, p2, v2, p1, 0, 0]';
coeff = R \ rhs;
ti = t;

%% 计算样条值
if t <= t1
    S = coeff(1) + coeff(2) * ti + coeff(3) * ti * ti + coeff(4) * ti * ti * ti;
    Sd = coeff(2) + 2 * coeff(3) * ti + 3 * coeff(4) * ti * ti;
    Sdd = 2 * coeff(3) + 6 * coeff(4) * ti;
else
    S = coeff(5) + coeff(6) * ti + coeff(7) * ti * ti + coeff(8) * ti * ti * ti;
    Sd = coeff(6) + 2 * coeff(7) * ti + 3 * coeff(8) * ti * ti;
    Sdd = 2 * coeff(7) + 6 * coeff(8) * ti;
end
