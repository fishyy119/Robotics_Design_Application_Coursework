function [KE, PEg, PEc, Pdiss_damping, penetration, vertex_z, vertex_z_dot] = energy_floating_base_multicontact(~, z, m, xG, yG, zG, I_body, g, lCuboid, wCuboid, hCuboid, ground_z, k_ground, c_ground)
% 本文件由 derive_floating_base_multicontact.m 自动生成。
% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。

% 状态量
q1 = z(1); u1 = z(2);
q2 = z(3); u2 = z(4);
q3 = z(5); u3 = z(6);
q4 = z(7); u4 = z(8);
q5 = z(9); u5 = z(10);
q6 = z(11); u6 = z(12);

% 惯量参数
I_world = I_body;

% 姿态、角速度与质心运动
Rz = [cos(q6), -sin(q6), 0; sin(q6), cos(q6), 0; 0, 0, 1];
Ry = [cos(q5), 0, sin(q5); 0, 1, 0; -sin(q5), 0, cos(q5)];
Rx = [1, 0, 0; 0, cos(q4), -sin(q4); 0, sin(q4), cos(q4)];
R = Rz * Ry * Rx;
I_world = R * I_world * R.';
r_B = [q1; q2; q3];
r_BG = R * [xG; yG; zG];
r_G = r_B + r_BG;
omega = R * [u4; 0; 0] + Rz * Ry * [0; u5; 0] + Rz * [0; 0; u6];
v_B = [u1; u2; u3];
v_G = v_B + cross(omega, r_BG);

% 顶点高度与接触弹簧势能
local_vertices = 0.5 * [ ...
    lCuboid, -lCuboid, -lCuboid, lCuboid, lCuboid, -lCuboid, -lCuboid, lCuboid; ...
    wCuboid,  wCuboid, -wCuboid, -wCuboid,  wCuboid,  wCuboid, -wCuboid, -wCuboid; ...
    hCuboid,  hCuboid,  hCuboid,  hCuboid, -hCuboid, -hCuboid, -hCuboid, -hCuboid];
vertex_offsets = R * local_vertices;
vertices = r_G + vertex_offsets;
vertex_velocities = repmat(v_G, 1, 8) + cross(repmat(omega, 1, 8), vertex_offsets, 1);
vertex_z = vertices(3, :).';
vertex_z_dot = vertex_velocities(3, :).';
penetration = max(ground_z - vertex_z, 0);
contact_active = penetration > 0;
approach_speed = max(-vertex_z_dot, 0) .* contact_active;
PEc = 0.5 * k_ground * sum(penetration .^ 2);
Pdiss_damping = c_ground * sum(approach_speed .^ 2);

% 能量
KE = 0.5 * m * dot(v_G, v_G) + 0.5 * omega.' * I_world * omega;
PEg = m * g * (r_G(3) - ground_z);
