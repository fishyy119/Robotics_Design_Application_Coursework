function derive_floating_base_multicontact()
% 通过符号推导生成三维浮动基座多顶点接触模型所需的数值函数。

clc;

%% 符号变量
syms m g xG yG zG lCuboid wCuboid hCuboid real
syms Ixx Iyy Izz Ixy Ixz Iyz real
syms q1 q2 q3 q4 q5 q6 u1 u2 u3 u4 u5 u6 ud1 ud2 ud3 ud4 ud5 ud6 real

q = [q1; q2; q3; q4; q5; q6];
u = [u1; u2; u3; u4; u5; u6];
ud = [ud1; ud2; ud3; ud4; ud5; ud6];

I_body = [Ixx, Ixy, Ixz; ...
    Ixy, Iyy, Iyz; ...
    Ixz, Iyz, Izz];

%% 参考坐标系
i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

Rz = rot_axis('z', q6);
Ry = rot_axis('y', q5);
Rx = rot_axis('x', q4);
R = Rz * Ry * Rx;
I_world = R * I_body * R.';

%% 几何关系
r_B = [q1; q2; q3];
r_B_G = R * [xG; yG; zG];
r_G = r_B + r_B_G;

local_vertices = 0.5 * [ ...
    lCuboid, -lCuboid, -lCuboid, lCuboid, lCuboid, -lCuboid, -lCuboid, lCuboid; ...
    wCuboid,  wCuboid, -wCuboid, -wCuboid,  wCuboid,  wCuboid, -wCuboid, -wCuboid; ...
    hCuboid,  hCuboid,  hCuboid,  hCuboid, -hCuboid, -hCuboid, -hCuboid, -hCuboid];

vertices = sym(zeros(3, 8));
for idx = 1:8
    vertices(:, idx) = r_G + R * local_vertices(:, idx);
end

vertex_z = vertices(3, :).';
grad_z = jacobian(vertex_z, q);

B_axes = [r_B + R * i, r_B + R * j, r_B + R * k];
G_axes = [r_G + R * i, r_G + R * j, r_G + R * k];

%% 速度与能量
om = R * (u4 * i) + Rz * Ry * (u5 * j) + Rz * (u6 * k);
v_B = [u1; u2; u3];
v_G = v_B + cross(om, r_B_G, 1);

KE = 0.5 * m * (v_G.' * v_G) + 0.5 * om.' * I_world * om;
PEg = m * g * r_G(3);

%% 拉格朗日动力学
L = KE - PEg;
Ld_u = jacobian(L, u).';
d_dt_Ld_u = jacobian(Ld_u, q) * u + jacobian(Ld_u, u) * ud;
Ld_q = jacobian(L, q).';
eom = d_dt_Ld_u - Ld_q;

M = sym(zeros(6, 6));
RHS_free = sym(zeros(6, 1));
zero_ud = sym(zeros(6, 1));

for row = 1:6
    eqn_row = collect(eom(row), ud);
    RHS_free(row) = -subs(eqn_row, ud, zero_ud);
    for col = 1:6
        basis = sym(zeros(6, 1));
        basis(col) = 1;
        M(row, col) = subs(eqn_row, ud, basis) + RHS_free(row);
    end
end

%% 自动生成数值函数
this_dir = fileparts(mfilename('fullpath'));

write_rhs_file(fullfile(this_dir, 'rhs_floating_base_multicontact.m'), M, RHS_free, vertex_z, grad_z);
write_energy_file(fullfile(this_dir, 'energy_floating_base_multicontact.m'));
write_draw_file(fullfile(this_dir, 'draw_floating_base_multicontact.m'), r_B, B_axes, r_G, G_axes, vertices);

disp('已生成 rhs_floating_base_multicontact.m、energy_floating_base_multicontact.m 和 draw_floating_base_multicontact.m');
end

function write_rhs_file(filename, M, RHS_free, vertex_z, grad_z)
% 写入状态导数函数文件。

fid = fopen(filename, 'w');
assert(fid ~= -1, 'Failed to open %s for writing.', filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'function zdot = rhs_floating_base_multicontact(~, z, m, xG, yG, zG, I_body, g, lCuboid, wCuboid, hCuboid, ground_z, k_ground, c_ground)\n');
fprintf(fid, '%% 本文件由 derive_floating_base_multicontact.m 自动生成。\n');
fprintf(fid, '%% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。\n\n');

fprintf(fid, '%% 状态量\n');
write_state_unpack(fid, 6);
fprintf(fid, '%% z(13) 为累计阻尼耗散能，仅用于能量诊断。\n');
fprintf(fid, '\n');

fprintf(fid, '%% 惯量参数\n');
write_inertia_unpack(fid);
fprintf(fid, '\n');

fprintf(fid, '%% 顶点高度\n');
write_vector_assignments(fid, 'vertex_z', vertex_z);
fprintf(fid, 'penetration = max(ground_z - vertex_z, 0);\n');
fprintf(fid, '\n');

fprintf(fid, '%% 顶点高度对广义坐标的雅可比\n');
write_matrix_assignments(fid, 'grad_z', grad_z);
fprintf(fid, 'u_vec = [u1; u2; u3; u4; u5; u6];\n');
fprintf(fid, 'vertex_z_dot = grad_z * u_vec;\n');
fprintf(fid, 'contact_active = penetration > 0;\n');
fprintf(fid, 'approach_speed = max(-vertex_z_dot, 0) .* contact_active;\n');
fprintf(fid, 'fz = k_ground * penetration + c_ground * approach_speed;\n');
fprintf(fid, 'damping_dissipation_rate = c_ground * sum(approach_speed .^ 2);\n');
fprintf(fid, 'contact_tau = grad_z.'' * fz;\n\n');

fprintf(fid, '%% 接触前自由动力学\n');
write_matrix_assignments(fid, 'MM', M);
write_vector_assignments(fid, 'RHS_free', RHS_free);
fprintf(fid, 'qdd = MM \\ (RHS_free + contact_tau);\n\n');

fprintf(fid, '%% 状态导数\n');
fprintf(fid, 'zdot = [u1; qdd(1); u2; qdd(2); u3; qdd(3); u4; qdd(4); u5; qdd(5); u6; qdd(6); damping_dissipation_rate];\n');
end

function write_energy_file(filename)
% 写入能量诊断函数文件。

fid = fopen(filename, 'w');
assert(fid ~= -1, 'Failed to open %s for writing.', filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'function [KE, PEg, PEc, Pdiss_damping, penetration, vertex_z, vertex_z_dot] = energy_floating_base_multicontact(~, z, m, xG, yG, zG, I_body, g, lCuboid, wCuboid, hCuboid, ground_z, k_ground, c_ground)\n');
fprintf(fid, '%% 本文件由 derive_floating_base_multicontact.m 自动生成。\n');
fprintf(fid, '%% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。\n\n');

fprintf(fid, '%% 状态量\n');
write_state_unpack(fid, 6);
fprintf(fid, '\n');

fprintf(fid, '%% 惯量参数\n');
fprintf(fid, 'I_world = I_body;\n');
fprintf(fid, '\n');

fprintf(fid, '%% 姿态、角速度与质心运动\n');
fprintf(fid, 'Rz = [cos(q6), -sin(q6), 0; sin(q6), cos(q6), 0; 0, 0, 1];\n');
fprintf(fid, 'Ry = [cos(q5), 0, sin(q5); 0, 1, 0; -sin(q5), 0, cos(q5)];\n');
fprintf(fid, 'Rx = [1, 0, 0; 0, cos(q4), -sin(q4); 0, sin(q4), cos(q4)];\n');
fprintf(fid, 'R = Rz * Ry * Rx;\n');
fprintf(fid, 'I_world = R * I_world * R.'';\n');
fprintf(fid, 'r_B = [q1; q2; q3];\n');
fprintf(fid, 'r_BG = R * [xG; yG; zG];\n');
fprintf(fid, 'r_G = r_B + r_BG;\n');
fprintf(fid, 'omega = R * [u4; 0; 0] + Rz * Ry * [0; u5; 0] + Rz * [0; 0; u6];\n');
fprintf(fid, 'v_B = [u1; u2; u3];\n');
fprintf(fid, 'v_G = v_B + cross(omega, r_BG);\n\n');

fprintf(fid, '%% 顶点高度与接触弹簧势能\n');
fprintf(fid, 'local_vertices = 0.5 * [ ...\n');
fprintf(fid, '    lCuboid, -lCuboid, -lCuboid, lCuboid, lCuboid, -lCuboid, -lCuboid, lCuboid; ...\n');
fprintf(fid, '    wCuboid,  wCuboid, -wCuboid, -wCuboid,  wCuboid,  wCuboid, -wCuboid, -wCuboid; ...\n');
fprintf(fid, '    hCuboid,  hCuboid,  hCuboid,  hCuboid, -hCuboid, -hCuboid, -hCuboid, -hCuboid];\n');
fprintf(fid, 'vertex_offsets = R * local_vertices;\n');
fprintf(fid, 'vertices = r_G + vertex_offsets;\n');
fprintf(fid, 'vertex_velocities = repmat(v_G, 1, 8) + cross(repmat(omega, 1, 8), vertex_offsets, 1);\n');
fprintf(fid, 'vertex_z = vertices(3, :).'';\n');
fprintf(fid, 'vertex_z_dot = vertex_velocities(3, :).'';\n');
fprintf(fid, 'penetration = max(ground_z - vertex_z, 0);\n');
fprintf(fid, 'contact_active = penetration > 0;\n');
fprintf(fid, 'approach_speed = max(-vertex_z_dot, 0) .* contact_active;\n');
fprintf(fid, 'PEc = 0.5 * k_ground * sum(penetration .^ 2);\n');
fprintf(fid, 'Pdiss_damping = c_ground * sum(approach_speed .^ 2);\n\n');

fprintf(fid, '%% 能量\n');
fprintf(fid, 'KE = 0.5 * m * dot(v_G, v_G) + 0.5 * omega.'' * I_world * omega;\n');
fprintf(fid, 'PEg = m * g * (r_G(3) - ground_z);\n');
end

function write_draw_file(filename, r_B, B_axes, r_G, G_axes, vertices)
% 写入几何绘制函数文件。

fid = fopen(filename, 'w');
assert(fid ~= -1, 'Failed to open %s for writing.', filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'function [r_B, B_axes, r_G, G_axes, vertices] = draw_floating_base_multicontact(~, z, xG, yG, zG, lCuboid, wCuboid, hCuboid)\n');
fprintf(fid, '%% 本文件由 derive_floating_base_multicontact.m 自动生成。\n');
fprintf(fid, '%% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。\n\n');

fprintf(fid, '%% 状态量\n');
write_state_unpack(fid, 6);
fprintf(fid, '\n');

fprintf(fid, '%% 参考点、质心和顶点位置\n');
write_vector_assignments(fid, 'r_B', r_B);
write_matrix_assignments(fid, 'B_axes', B_axes);
write_vector_assignments(fid, 'r_G', r_G);
write_matrix_assignments(fid, 'G_axes', G_axes);
write_matrix_assignments(fid, 'vertices', vertices);
end

function write_state_unpack(fid, n_q)
% 写入状态展开代码。

for idx = 1:n_q
    fprintf(fid, 'q%d = z(%d); u%d = z(%d);\n', idx, 2 * idx - 1, idx, 2 * idx);
end
end

function write_inertia_unpack(fid)
% 写入惯量矩阵展开代码。

fprintf(fid, 'Ixx = I_body(1, 1);\n');
fprintf(fid, 'Iyy = I_body(2, 2);\n');
fprintf(fid, 'Izz = I_body(3, 3);\n');
fprintf(fid, 'Ixy = I_body(1, 2);\n');
fprintf(fid, 'Ixz = I_body(1, 3);\n');
fprintf(fid, 'Iyz = I_body(2, 3);\n');
end

function write_vector_assignments(fid, name, vec)
% 写入列向量赋值代码。

n = numel(vec);
fprintf(fid, '%s = zeros(%d, 1);\n', name, n);
for idx = 1:n
    fprintf(fid, '%s(%d) = %s;\n', name, idx, expr_to_code(vec(idx)));
end
fprintf(fid, '\n');
end

function write_matrix_assignments(fid, name, mat)
% 写入矩阵赋值代码。

[n_row, n_col] = size(mat);
fprintf(fid, '%s = zeros(%d, %d);\n', name, n_row, n_col);
for row = 1:n_row
    for col = 1:n_col
        fprintf(fid, '%s(%d, %d) = %s;\n', name, row, col, expr_to_code(mat(row, col)));
    end
end
fprintf(fid, '\n');
end

function code = expr_to_code(expr)
% 将符号表达式转换为单行 MATLAB 代码。

code = char(expr);
code = strrep(code, sprintf('\n'), ' ');
end

function R = rot_axis(axis_name, angle)
% 返回绕坐标轴的旋转矩阵。

switch axis_name
    case 'x'
        R = [1, 0, 0; ...
            0, cos(angle), -sin(angle); ...
            0, sin(angle), cos(angle)];
    case 'y'
        R = [cos(angle), 0, sin(angle); ...
            0, 1, 0; ...
            -sin(angle), 0, cos(angle)];
    case 'z'
        R = [cos(angle), -sin(angle), 0; ...
            sin(angle), cos(angle), 0; ...
            0, 0, 1];
    otherwise
        error('Unsupported axis "%s".', axis_name);
end
end
