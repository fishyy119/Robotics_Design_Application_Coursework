function [x, y] = forward_kinematics(q, a)
% 计算末端在平面内的位置轨迹。

x = a(1) * cos(q(1, :)) + a(2) * cos(q(1, :) + q(2, :));
y = a(1) * sin(q(1, :)) + a(2) * sin(q(1, :) + q(2, :));
end
