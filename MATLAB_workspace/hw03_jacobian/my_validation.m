%% 平面二连杆末端速度（使用三维叉乘）
clear;
clc;

% --- 符号变量 ---
syms q1 q2 dq1 dq2 l1 l2 real

% 关节角和速度
q = [q1; q2];
dq = [dq1; dq2];

% 单位向量
i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

% 连杆方向（扩展为3D z=0）
e1 = sin(q1) * i - cos(q1) * j;
e2 = sin(q1 + q2) * i - cos(q1 + q2) * j;
L = [l1; l2];
P0 = [0; 0; 0];
P1 = P0 + L(1) * e1;
P2 = P1 + L(2) * e2;

% --- 速度推导 ---
V0 = [0; 0; 0]; % 基座速度

% 第一连杆末端速度
omega1 = [0; 0; dq1]; % 绕z轴旋转
V1 = V0 + cross(omega1, P1 - P0);

% 第二连杆末端速度（末端执行器）
omega2 = [0; 0; dq2 + dq1]; % 绕z轴旋转
V2 = V1 + cross(omega2, P2 - P1);

% 只取 x,y 分量作为平面速度
V_ee = simplify(V2(1:2))

% --- 与雅可比法对比 ---
J = jacobian(P2(1:2), q); % 末端位置对关节角的雅可比
V_jac = J * dq;
V_jac = simplify(V_jac)

% --- 差值检查 ---
diff = simplify(V_ee - V_jac)
