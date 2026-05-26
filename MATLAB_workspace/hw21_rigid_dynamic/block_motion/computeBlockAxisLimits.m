function axisLimits = computeBlockAxisLimits(result, bodyRadius, axisPadding)
% 根据质心和关键点轨迹分别构造长方体三轴显示范围。

allPoints = result.rComHistory;
for pointIdx = 1:size(result.trackedPointHistory, 3)
    allPoints = [allPoints, result.trackedPointHistory(:, :, pointIdx)]; %#ok<AGROW>
end

mins = min(allPoints, [], 2) - (bodyRadius + axisPadding);
maxs = max(allPoints, [], 2) + (bodyRadius + axisPadding);
axisLimits = [mins, maxs];
end
