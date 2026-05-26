function pointsWorld = transformBodyPoints(rotationMatrix, translation, pointsBody)
% 将刚体坐标系中的点转换到世界坐标系。

pointCount = size(pointsBody, 2);
pointsWorld = rotationMatrix * pointsBody + translation(:) * ones(1, pointCount);
end
