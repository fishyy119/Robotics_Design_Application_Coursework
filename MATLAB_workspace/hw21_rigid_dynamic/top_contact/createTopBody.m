function body = createTopBody(radius, diskThickness, tipOffset, density, shaftRadius)
% 创建以尖端为局部原点的陀螺刚体模型。

if nargin < 5
    shaftRadius = 0.01;
end

[diskVertices, diskFaces] = createZCylinder([0; 0; tipOffset], radius, diskThickness, 24);
[shaftVertices, shaftFaces] = createZCylinder([0; 0; tipOffset], shaftRadius, 2 * tipOffset, 18);

body = struct();
body.name = 'Top';
body.mass = pi * radius^2 * diskThickness * density;
body.comBody = [0; 0; tipOffset];
body.inertiaBody = diag([
    (diskThickness^2 + 3 * radius^2) / 12, ...
    (diskThickness^2 + 3 * radius^2) / 12, ...
    radius^2 / 2
]) * body.mass;
body.tipBody = [0; 0; 0];
body.topPointBody = [0; 0; 2 * tipOffset];
body.vertices = [diskVertices, shaftVertices];
body.faces = [diskFaces, shaftFaces + size(diskVertices, 2)];
end
