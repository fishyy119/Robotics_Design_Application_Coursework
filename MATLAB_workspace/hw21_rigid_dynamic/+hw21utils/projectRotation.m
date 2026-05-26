function R = projectRotation(RApprox)
% 将数值积分后的矩阵投影回 SO(3)。

[U, ~, V] = svd(RApprox);
R = U * diag([1, 1, det(U * V.')]) * V.';
end
