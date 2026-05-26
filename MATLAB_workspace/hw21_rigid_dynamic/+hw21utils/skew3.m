function S = skew3(vector)
% 将三维向量转换为反对称矩阵。

vector = vector(:);
S = [
    0, -vector(3),  vector(2);
    vector(3), 0, -vector(1);
    -vector(2), vector(1), 0
];
end
