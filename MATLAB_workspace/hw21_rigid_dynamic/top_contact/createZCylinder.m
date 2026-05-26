function [vertices, faces] = createZCylinder(center, radius, bodyLength, sideCount)
% 生成沿 z 轴放置的圆柱网格。

if nargin < 4
    sideCount = 16;
end

theta = (0:sideCount - 1) / sideCount * 2 * pi;
x = radius * cos(theta);
y = radius * sin(theta);
zTop = bodyLength / 2 * ones(1, sideCount);
zBottom = -zTop;

vertices = [
    x, x, 0, 0;
    y, y, 0, 0;
    zTop, zBottom, bodyLength / 2, -bodyLength / 2
];

vertices = vertices + center(:);

sideFaces = [
    1:sideCount;
    sideCount + (1:sideCount);
    sideCount + [2:sideCount, 1];
    [2:sideCount, 1]
];

topFaces = [
    1:sideCount;
    [2:sideCount, 1];
    repmat(2 * sideCount + 1, 1, sideCount);
    repmat(2 * sideCount + 1, 1, sideCount)
];

bottomFaces = [
    sideCount + [2:sideCount, 1];
    sideCount + (1:sideCount);
    repmat(2 * sideCount + 2, 1, sideCount);
    repmat(2 * sideCount + 2, 1, sideCount)
];

faces = [sideFaces, topFaces, bottomFaces];
end
