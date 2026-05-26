clearvars;
close all;
clc;
utils.setDefaultGraphics;
%% 路径
scriptDir = fileparts(mfilename('fullpath'));
hw21Dir = fileparts(scriptDir);
videoDir = fullfile(hw21Dir, 'videos');

if ~exist(videoDir, 'dir')
    mkdir(videoDir);
end
%% 刚体参数
topBody = createTopBody(0.20, 0.05, 0.20, 2.7e3, 0.01);
tiltAngle = pi / 50;
R0Top = rotationExp([tiltAngle; 0; 0]);
tipHeight0 = 0.30;
%% 初始条件
caseConfig = struct();
caseConfig.name = 'Top with Ground Contact';
caseConfig.body = topBody;
caseConfig.gravity = [0; 0; -9.81];
caseConfig.contactStiffness = 1.0e4;
caseConfig.contactNormalDamping = 1.0e3;
caseConfig.contactTangentialDamping = 1.0e3;
caseConfig.origin0 = [0; 0; tipHeight0];
caseConfig.vSpatial0 = [0; 0; 0];
caseConfig.R0 = R0Top;
caseConfig.wWorld0 = [0; 0; 50.0];
caseConfig.dt = 0.001;
caseConfig.tEnd = 15.00;
caseConfig.trackedPointsBody = topBody.topPointBody;
caseConfig.trackedPointNames = {'Top end'};
%% 数值积分
result = simulateTopCase(caseConfig);
%% 视频导出
frameRate = 30;
viewAngle = [34, 18];
axisPadding = 0.10;
bodyRadius = max(vecnorm(topBody.vertices, 2, 1));
frameStride = max(1, round(1/(frameRate * caseConfig.dt)));
frameIds = unique([1:frameStride:numel(result.time), numel(result.time)]);
colors = lines(max(1, size(result.trackedPointHistory, 3)));
poseTranslationHistory = result.originHistory;

figParams = struct('Name', result.name, 'Width', 18, 'AspectRatio', 0.75);
fig = utils.createFigureA4(figParams);
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
view(ax, viewAngle(1), viewAngle(2));
xlim(ax, [-0.2, 0.6]);
ylim(ax, [-0.4, 0.4]);
zlim(ax, [0, 0.8]);
xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
zlabel(ax, 'z (m)');

patchHandle = patch( ...
    'Parent', ax, ...
    'Faces', topBody.faces.', ...
    'Vertices', zeros(size(topBody.vertices, 2), 3), ...
    'FaceColor', [0.70, 0.70, 0.78], ...
    'FaceAlpha', 0.45, ...
    'EdgeColor', [0.25, 0.25, 0.30], ...
    'LineWidth', 0.9, ...
    'HandleVisibility', 'off');

comLine = plot3(ax, nan, nan, nan, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Center of mass');
comMarker = plot3(ax, nan, nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');

trackedLineHandles = gobjects(size(result.trackedPointHistory, 3), 1);
trackedMarkerHandles = gobjects(size(result.trackedPointHistory, 3), 1);
for pointIdx = 1:size(result.trackedPointHistory, 3)
    trackedLineHandles(pointIdx) = plot3( ...
        ax, nan, nan, nan, '-', ...
        'Color', colors(pointIdx, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', result.trackedPointNames{pointIdx});
    trackedMarkerHandles(pointIdx) = plot3( ...
        ax, nan, nan, nan, 'o', ...
        'Color', colors(pointIdx, :), ...
        'MarkerFaceColor', colors(pointIdx, :), ...
        'MarkerSize', 5, ...
        'HandleVisibility', 'off');
end

legend(ax, 'Location', 'best');
title(ax, result.name);

videoPath = fullfile(videoDir, 'hw21_top_contact.mp4');
videoWriter = VideoWriter(videoPath, 'MPEG-4');
videoWriter.FrameRate = frameRate;
if isprop(videoWriter, 'Quality')
    videoWriter.Quality = 95;
end
open(videoWriter);

for frameIdx = frameIds
    worldVertices = hw21utils.transformBodyPoints( ...
        result.rotationHistory(:, :, frameIdx), ...
        poseTranslationHistory(:, frameIdx), ...
        topBody.vertices ...
        );
    set(patchHandle, 'Vertices', worldVertices.');
    set(comLine, ...
        'XData', result.rComHistory(1, 1:frameIdx), ...
        'YData', result.rComHistory(2, 1:frameIdx), ...
        'ZData', result.rComHistory(3, 1:frameIdx));
    set(comMarker, ...
        'XData', result.rComHistory(1, frameIdx), ...
        'YData', result.rComHistory(2, frameIdx), ...
        'ZData', result.rComHistory(3, frameIdx));

    for pointIdx = 1:size(result.trackedPointHistory, 3)
        pointHistory = result.trackedPointHistory(:, :, pointIdx);
        set(trackedLineHandles(pointIdx), ...
            'XData', pointHistory(1, 1:frameIdx), ...
            'YData', pointHistory(2, 1:frameIdx), ...
            'ZData', pointHistory(3, 1:frameIdx));
        set(trackedMarkerHandles(pointIdx), ...
            'XData', pointHistory(1, frameIdx), ...
            'YData', pointHistory(2, frameIdx), ...
            'ZData', pointHistory(3, frameIdx));
    end

    title(ax, sprintf('%s, t = %.2f s', result.name, result.time(frameIdx)));
    drawnow;
    writeVideo(videoWriter, getframe(fig));
end

close(videoWriter);
close(fig);
