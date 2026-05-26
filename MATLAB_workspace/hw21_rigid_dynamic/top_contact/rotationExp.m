function R = rotationExp(rotationVector)
% 由旋转向量计算旋转矩阵。

rotationVector = rotationVector(:);
theta = norm(rotationVector);

if theta < 1.0e-12
    R = eye(3) + hw21utils.skew3(rotationVector);
    return;
end

axis = rotationVector / theta;
axisHat = hw21utils.skew3(axis);
R = eye(3) + sin(theta) * axisHat + (1 - cos(theta)) * (axisHat * axisHat);
end
