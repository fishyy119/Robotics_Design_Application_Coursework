function body = createBoxBody(dimensions, mass)
% 创建以质心为原点的长方体刚体模型。

dimensions = dimensions(:);
halfDim = dimensions / 2;

body = struct();
body.name = 'Block';
body.mass = mass;
body.dimensions = dimensions;
body.comBody = [0; 0; 0];
body.inertiaBody = diag([
    1 / 12 * (dimensions(2)^2 + dimensions(3)^2), ...
    1 / 12 * (dimensions(1)^2 + dimensions(3)^2), ...
    1 / 12 * (dimensions(1)^2 + dimensions(2)^2)
]) * mass;
body.tipBody = [];
body.referenceCornerBody = halfDim;

body.vertices = [
    -halfDim(1), -halfDim(1),  halfDim(1),  halfDim(1), -halfDim(1), -halfDim(1),  halfDim(1),  halfDim(1);
    -halfDim(2),  halfDim(2),  halfDim(2), -halfDim(2), -halfDim(2),  halfDim(2),  halfDim(2), -halfDim(2);
    -halfDim(3), -halfDim(3), -halfDim(3), -halfDim(3),  halfDim(3),  halfDim(3),  halfDim(3),  halfDim(3)
];

body.faces = [
    1, 2, 3, 4;
    2, 6, 7, 3;
    4, 3, 7, 8;
    1, 5, 8, 4;
    1, 2, 6, 5;
    5, 6, 7, 8
].';
end
