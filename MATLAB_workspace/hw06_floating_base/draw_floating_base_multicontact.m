function [r_B, B_axes, r_G, G_axes, vertices] = draw_floating_base_multicontact(~, z, xG, yG, zG, lCuboid, wCuboid, hCuboid)
% 本文件由 derive_floating_base_multicontact.m 自动生成。
% 如需修改推导过程，请编辑生成脚本而非直接手工修改本文件。

% 状态量
q1 = z(1); u1 = z(2);
q2 = z(3); u2 = z(4);
q3 = z(5); u3 = z(6);
q4 = z(7); u4 = z(8);
q5 = z(9); u5 = z(10);
q6 = z(11); u6 = z(12);

% 参考点、质心和顶点位置
r_B = zeros(3, 1);
r_B(1) = q1;
r_B(2) = q2;
r_B(3) = q3;

B_axes = zeros(3, 3);
B_axes(1, 1) = q1 + cos(q5)*cos(q6);
B_axes(1, 2) = q1 - cos(q4)*sin(q6) + cos(q6)*sin(q4)*sin(q5);
B_axes(1, 3) = q1 + sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5);
B_axes(2, 1) = q2 + cos(q5)*sin(q6);
B_axes(2, 2) = q2 + cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6);
B_axes(2, 3) = q2 - cos(q6)*sin(q4) + cos(q4)*sin(q5)*sin(q6);
B_axes(3, 1) = q3 - sin(q5);
B_axes(3, 2) = q3 + cos(q5)*sin(q4);
B_axes(3, 3) = q3 + cos(q4)*cos(q5);

r_G = zeros(3, 1);
r_G(1) = q1 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + xG*cos(q5)*cos(q6);
r_G(2) = q2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + xG*cos(q5)*sin(q6);
r_G(3) = q3 - xG*sin(q5) + zG*cos(q4)*cos(q5) + yG*cos(q5)*sin(q4);

G_axes = zeros(3, 3);
G_axes(1, 1) = q1 + cos(q5)*cos(q6) - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + xG*cos(q5)*cos(q6);
G_axes(1, 2) = q1 - cos(q4)*sin(q6) - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + xG*cos(q5)*cos(q6) + cos(q6)*sin(q4)*sin(q5);
G_axes(1, 3) = q1 + sin(q4)*sin(q6) - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + xG*cos(q5)*cos(q6) + cos(q4)*cos(q6)*sin(q5);
G_axes(2, 1) = q2 + cos(q5)*sin(q6) + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + xG*cos(q5)*sin(q6);
G_axes(2, 2) = q2 + cos(q4)*cos(q6) + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + sin(q4)*sin(q5)*sin(q6) + xG*cos(q5)*sin(q6);
G_axes(2, 3) = q2 - cos(q6)*sin(q4) + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + xG*cos(q5)*sin(q6) + cos(q4)*sin(q5)*sin(q6);
G_axes(3, 1) = q3 - sin(q5) - xG*sin(q5) + zG*cos(q4)*cos(q5) + yG*cos(q5)*sin(q4);
G_axes(3, 2) = q3 + cos(q5)*sin(q4) - xG*sin(q5) + zG*cos(q4)*cos(q5) + yG*cos(q5)*sin(q4);
G_axes(3, 3) = q3 + cos(q4)*cos(q5) - xG*sin(q5) + zG*cos(q4)*cos(q5) + yG*cos(q5)*sin(q4);

vertices = zeros(3, 8);
vertices(1, 1) = q1 + (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 - (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 2) = q1 + (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 - (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) - (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 3) = q1 + (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 + (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) - (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 4) = q1 + (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 + (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 5) = q1 - (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 - (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 6) = q1 - (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 - (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) - (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 7) = q1 - (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 + (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) - (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(1, 8) = q1 - (hCuboid*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)))/2 + (wCuboid*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)))/2 - yG*(cos(q4)*sin(q6) - cos(q6)*sin(q4)*sin(q5)) + zG*(sin(q4)*sin(q6) + cos(q4)*cos(q6)*sin(q5)) + (lCuboid*cos(q5)*cos(q6))/2 + xG*cos(q5)*cos(q6);
vertices(2, 1) = q2 - (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 + (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 2) = q2 - (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 + (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) - (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 3) = q2 - (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 - (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) - (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 4) = q2 - (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 - (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 5) = q2 + (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 + (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 6) = q2 + (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 + (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) - (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 7) = q2 + (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 - (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) - (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(2, 8) = q2 + (hCuboid*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)))/2 - (wCuboid*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)))/2 + yG*(cos(q4)*cos(q6) + sin(q4)*sin(q5)*sin(q6)) - zG*(cos(q6)*sin(q4) - cos(q4)*sin(q5)*sin(q6)) + (lCuboid*cos(q5)*sin(q6))/2 + xG*cos(q5)*sin(q6);
vertices(3, 1) = q3 - (lCuboid*sin(q5))/2 - xG*sin(q5) + (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) + (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 2) = q3 + (lCuboid*sin(q5))/2 - xG*sin(q5) + (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) + (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 3) = q3 + (lCuboid*sin(q5))/2 - xG*sin(q5) + (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) - (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 4) = q3 - (lCuboid*sin(q5))/2 - xG*sin(q5) + (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) - (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 5) = q3 - (lCuboid*sin(q5))/2 - xG*sin(q5) - (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) + (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 6) = q3 + (lCuboid*sin(q5))/2 - xG*sin(q5) - (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) + (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 7) = q3 + (lCuboid*sin(q5))/2 - xG*sin(q5) - (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) - (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);
vertices(3, 8) = q3 - (lCuboid*sin(q5))/2 - xG*sin(q5) - (hCuboid*cos(q4)*cos(q5))/2 + zG*cos(q4)*cos(q5) - (wCuboid*cos(q5)*sin(q4))/2 + yG*cos(q5)*sin(q4);

