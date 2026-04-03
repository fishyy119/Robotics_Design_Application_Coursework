function derive_planar_contact()
% 通过符号推导生成平面二连杆接触模型所需的动力学函数。

clc;

%% 符号变量
syms l1 l2 a1 a2 m1 m2 J1 J2 g T1 T2 Fx Fy real
syms q1 q2 u1 u2 ud1 ud2 real

q = [q1; q2];
u = [u1; u2];
ud = [ud1; ud2];

%% 参考坐标系
i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

R1 = [cos(q1), -sin(q1), 0;
    sin(q1), cos(q1), 0;
    0, 0, 1];
R2 = [cos(q1 + q2), -sin(q1 + q2), 0;
    sin(q1 + q2), cos(q1 + q2), 0;
    0, 0, 1];

e1 = R1(:, 1);
e2 = R2(:, 1);

%% 位置、速度与加速度
r_O_G1 = a1 * e1;
r_O_E = l1 * e1;
r_E_G2 = a2 * e2;
r_E_F = l2 * e2;
r_O_G2 = r_O_E + r_E_G2;
r_O_F = r_O_E + r_E_F;

om1 = u1 * k;
om2 = (u1 + u2) * k;
al1 = ud1 * k;
al2 = (ud1 + ud2) * k;

v_O = sym(zeros(3, 1));
v_G1 = v_O + cross(om1, r_O_G1, 1);
v_E = v_O + cross(om1, r_O_E, 1);
v_G2 = v_E + cross(om2, r_E_G2, 1);
v_F = v_E + cross(om2, r_E_F, 1);

a_O = sym(zeros(3, 1));
a_G1 = a_O + cross(om1, cross(om1, r_O_G1, 1), 1) + cross(al1, r_O_G1, 1);
a_E = a_O + cross(om1, cross(om1, r_O_E, 1), 1) + cross(al1, r_O_E, 1);
a_G2 = a_E + cross(om2, cross(om2, r_E_G2, 1), 1) + cross(al2, r_E_G2, 1);

%% Newton-Euler 力矩平衡
F_tip = [Fx; Fy; 0];
F_g1 = -m1 * g * j;
F_g2 = -m2 * g * j;

M_O = dot(cross(r_O_G1, F_g1, 1) + cross(r_O_G2, F_g2, 1) + cross(r_O_F, F_tip, 1) + T1 * k, k);
Hdot_O = dot(cross(r_O_G1, m1 * a_G1, 1) + J1 * al1 + ...
    cross(r_O_G2, m2 * a_G2, 1) + J2 * al2, k);

M_E = dot(cross(r_E_G2, F_g2, 1) + cross(r_E_F, F_tip, 1) + T2 * k, k);
Hdot_E = dot(cross(r_E_G2, m2 * a_G2, 1) + J2 * al2, k);

eqn1 = collect(expand(M_O - Hdot_O), [ud1, ud2]);
eqn2 = collect(expand(M_E - Hdot_E), [ud1, ud2]);

RHS1 = expand(-subs(eqn1, [ud1, ud2], [0, 0]));
M11 = expand(subs(eqn1, [ud1, ud2], [1, 0]) + RHS1);
M12 = expand(subs(eqn1, [ud1, ud2], [0, 1]) + RHS1);

RHS2 = expand(-subs(eqn2, [ud1, ud2], [0, 0]));
M21 = expand(subs(eqn2, [ud1, ud2], [1, 0]) + RHS2);
M22 = expand(subs(eqn2, [ud1, ud2], [0, 1]) + RHS2);

M = [M11, M12; M21, M22];
RHS = [RHS1; RHS2];

%% 能量表达式
KE = expand(0.5 * m1 * dot(v_G1, v_G1) + 0.5 * J1 * u1^2 + ...
    0.5 * m2 * dot(v_G2, v_G2) + 0.5 * J2 * (u1 + u2)^2);
PEg = expand(m1 * g * r_O_G1(2) + m2 * g * r_O_G2(2));
y_tip = expand(r_O_F(2));
y_tip_dot = expand(v_F(2));

%% 写入自动生成文件
this_dir = fileparts(mfilename('fullpath'));
write_rhs_file(fullfile(this_dir, 'rhs_planar_contact.m'), M, RHS, y_tip, y_tip_dot);
write_energy_file(fullfile(this_dir, 'energy_planar_contact.m'), KE, PEg, y_tip, y_tip_dot);

disp('已生成 rhs_planar_contact.m 和 energy_planar_contact.m');
end

function write_rhs_file(filename, M, RHS, y_tip, y_tip_dot)
% 写入状态导数函数文件。

fid = fopen(filename, 'w');
assert(fid ~= -1, 'Failed to open %s for writing.', filename);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'function zdot = rhs_planar_contact(~, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground)\n');
fprintf(fid, '%% 本文件由 derive_planar_contact.m 自动生成。\n');
fprintf(fid, '%% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。\n\n');
fprintf(fid, '%% 状态量\n');
fprintf(fid, 'q1 = z(1); u1 = z(2);\n');
fprintf(fid, 'q2 = z(3); u2 = z(4);\n\n');
fprintf(fid, '%% 接触量\n');
fprintf(fid, 'tip_y = %s;\n', char(y_tip));
fprintf(fid, 'tip_y_dot = %s;\n', char(y_tip_dot));
fprintf(fid, 'penetration = max(ground_y - tip_y, 0);\n');
fprintf(fid, 'approach_speed = max(-tip_y_dot, 0);\n\n');
fprintf(fid, '%% 接触力\n');
fprintf(fid, 'Fx = 0;\n');
fprintf(fid, 'Fy = 0;\n');
fprintf(fid, 'if penetration > 0\n');
fprintf(fid, '    Fy = k_ground * penetration + c_ground * approach_speed;\n');
fprintf(fid, 'end\n\n');
fprintf(fid, '%% 外加关节力矩\n');
fprintf(fid, 'T1 = 0;\n');
fprintf(fid, 'T2 = 0;\n\n');
fprintf(fid, '%% 动力学方程\n');
fprintf(fid, 'M11 = %s;\n', char(M(1, 1)));
fprintf(fid, 'M12 = %s;\n', char(M(1, 2)));
fprintf(fid, 'M21 = %s;\n', char(M(2, 1)));
fprintf(fid, 'M22 = %s;\n\n', char(M(2, 2)));
fprintf(fid, 'RHS1 = %s;\n', char(RHS(1)));
fprintf(fid, 'RHS2 = %s;\n\n', char(RHS(2)));
fprintf(fid, 'MM = [M11, M12; M21, M22];\n');
fprintf(fid, 'qdd = MM \\ [RHS1; RHS2];\n\n');
fprintf(fid, 'zdot = [u1; qdd(1); u2; qdd(2)];\n');
end

function write_energy_file(filename, KE, PEg, y_tip, y_tip_dot)
% 写入能量诊断函数文件。

fid = fopen(filename, 'w');
assert(fid ~= -1, 'Failed to open %s for writing.', filename);
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'function [KE, PEg, PEc, Pdamp, Fy, penetration, tip_y, tip_y_dot] = energy_planar_contact(~, z, m1, m2, l1, l2, a1, a2, J1, J2, g, ground_y, k_ground, c_ground)\n');
fprintf(fid, '%% 本文件由 derive_planar_contact.m 自动生成。\n');
fprintf(fid, '%% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。\n\n');
fprintf(fid, '%% 状态量\n');
fprintf(fid, 'q1 = z(1); u1 = z(2);\n');
fprintf(fid, 'q2 = z(3); u2 = z(4);\n\n');
fprintf(fid, '%% 接触量\n');
fprintf(fid, 'tip_y = %s;\n', char(y_tip));
fprintf(fid, 'tip_y_dot = %s;\n', char(y_tip_dot));
fprintf(fid, 'penetration = max(ground_y - tip_y, 0);\n');
fprintf(fid, 'approach_speed = max(-tip_y_dot, 0);\n');
fprintf(fid, 'Fy = 0;\n');
fprintf(fid, 'if penetration > 0\n');
fprintf(fid, '    Fy = k_ground * penetration + c_ground * approach_speed;\n');
fprintf(fid, 'end\n');
fprintf(fid, 'PEc = 0.5 * k_ground * penetration^2;\n');
fprintf(fid, 'Pdamp = c_ground * approach_speed^2;\n\n');
fprintf(fid, '%% 能量表达式\n');
fprintf(fid, 'KE = %s;\n', char(KE));
fprintf(fid, 'PEg = %s;\n', char(PEg));
end
