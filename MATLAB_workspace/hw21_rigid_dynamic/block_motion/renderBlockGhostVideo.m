function renderBlockGhostVideo(result, blockBody, caseConfig, videoDir, bodyRadius)
% 渲染固定间隔拖影样式的长方体视频。

frameRate = 30;
frameStride = max(1, round(1 / (frameRate * caseConfig.dt)));
frameIds = unique([1:frameStride:numel(result.time), numel(result.time)]);
axisLimits = computeBlockAxisLimits(result, bodyRadius, caseConfig.axisPadding);
poseTranslationHistory = result.originHistory;
ghostStepSamples = max(1, round(caseConfig.ghostTimeStep / caseConfig.dt));
ghostSampleIds = 1:ghostStepSamples:numel(result.time);
ghostColor = [0.45, 0.55, 0.72];
ghostAlphaNewest = 0.24;
ghostAlphaOldest = 0.08;
ghostEdgeColorScale = 0.65;
ghostLineWidth = 0.7;
currentFaceColor = [0.70, 0.70, 0.78];
currentFaceAlpha = 0.62;
currentEdgeColor = [0.20, 0.20, 0.24];
currentLineWidth = 1.0;
ghostZWindowSpan = max(caseConfig.zWindowSpan, 4 * (bodyRadius + 0.20));

figParams = struct('Name', [result.name, ' ghost'], 'Width', 18, 'AspectRatio', 0.75);
fig = utils.createFigureA4(figParams);
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, caseConfig.viewAngle(1), caseConfig.viewAngle(2));
xlim(ax, axisLimits(1, :));
ylim(ax, axisLimits(2, :));
if caseConfig.followZWindow
    zWindowCenterHistory = result.rComHistory(3, :);
    zWindowSpan = ghostZWindowSpan;
    zlim(ax, zWindowCenterHistory(1) + 0.5 * zWindowSpan * [-1, 1]);
else
    zlim(ax, axisLimits(3, :));
end
xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
zlabel(ax, 'z (m)');

ghostPatchHandles = gobjects(numel(ghostSampleIds), 1);
for ghostIdx = 1:numel(ghostSampleIds)
    ghostPatchHandles(ghostIdx) = patch( ...
        'Parent', ax, ...
        'Faces', blockBody.faces.', ...
        'Vertices', zeros(size(blockBody.vertices, 2), 3), ...
        'FaceColor', ghostColor, ...
        'FaceAlpha', ghostAlphaOldest, ...
        'EdgeColor', min(1.0, ghostColor * ghostEdgeColorScale), ...
        'LineWidth', ghostLineWidth, ...
        'Visible', 'off', ...
        'HandleVisibility', 'off');
end

currentPatchHandle = patch( ...
    'Parent', ax, ...
    'Faces', blockBody.faces.', ...
    'Vertices', zeros(size(blockBody.vertices, 2), 3), ...
    'FaceColor', currentFaceColor, ...
    'FaceAlpha', currentFaceAlpha, ...
    'EdgeColor', currentEdgeColor, ...
    'LineWidth', currentLineWidth, ...
    'HandleVisibility', 'off');
title(ax, [result.name, ' ghost trail']);

videoPath = fullfile(videoDir, caseConfig.videoNames.ghost);
videoWriter = VideoWriter(videoPath, 'MPEG-4');
videoWriter.FrameRate = frameRate;
if isprop(videoWriter, 'Quality')
    videoWriter.Quality = 95;
end
open(videoWriter);

for frameIdx = frameIds
    currentVertices = hw21utils.transformBodyPoints( ...
        result.rotationHistory(:, :, frameIdx), ...
        poseTranslationHistory(:, frameIdx), ...
        blockBody.vertices ...
    );
    if caseConfig.followZWindow
        zlim(ax, zWindowCenterHistory(frameIdx) + 0.5 * zWindowSpan * [-1, 1]);
    end
    set(currentPatchHandle, 'Vertices', currentVertices.');

    visibleGhostIds = ghostSampleIds(ghostSampleIds < frameIdx);
    visibleGhostCount = numel(visibleGhostIds);
    for ghostIdx = 1:numel(ghostSampleIds)
        if ghostIdx > visibleGhostCount
            set(ghostPatchHandles(ghostIdx), 'Visible', 'off');
            continue;
        end

        sampleIdx = visibleGhostIds(ghostIdx);
        ghostVertices = hw21utils.transformBodyPoints( ...
            result.rotationHistory(:, :, sampleIdx), ...
            poseTranslationHistory(:, sampleIdx), ...
            blockBody.vertices ...
        );
        if visibleGhostCount == 1
            ghostAlpha = ghostAlphaNewest;
        else
            alphaRatio = (ghostIdx - 1) / (visibleGhostCount - 1);
            ghostAlpha = ghostAlphaOldest + ...
                (ghostAlphaNewest - ghostAlphaOldest) * alphaRatio;
        end
        set(ghostPatchHandles(ghostIdx), ...
            'Vertices', ghostVertices.', ...
            'FaceAlpha', ghostAlpha, ...
            'Visible', 'on');
    end

    title(ax, sprintf('%s ghost trail, t = %.2f s', result.name, result.time(frameIdx)));
    drawnow;
    writeVideo(videoWriter, getframe(fig));
end

close(videoWriter);
close(fig);
end
