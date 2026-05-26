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
blockBody = createBoxBody([0.10; 0.40; 0.90], 36.0);
bodyRadius = max(vecnorm(blockBody.vertices, 2, 1));

%% 运动学配置
caseConfig = struct();
caseConfig.name = 'Block with Constant Spatial Velocity';
caseConfig.videoTypes = {'ghost'};
caseConfig.videoNames = struct( ...
    'ghost', 'hw21_block_spatial_kinematics.mp4');
caseConfig.origin0 = [0; 0; 0];
caseConfig.vSpatial0 = [0.3; 0; 1];
caseConfig.R0 = eye(3);
caseConfig.wWorld0 = [1; 0.0; 0.0];
caseConfig.dt = 0.01;
caseConfig.tEnd = 10.0;
caseConfig.viewAngle = [-34, 18];
caseConfig.axisPadding = 0.20;
caseConfig.followZWindow = false;
caseConfig.zWindowSpan = 0.0;
caseConfig.ghostTimeStep = 0.20;
caseConfig.body = blockBody;
caseConfig.trackedPointsBody = blockBody.referenceCornerBody;
caseConfig.trackedPointNames = {'Upper corner'};

%% 运动学积分
result = simulateBlockConstantSpatialKinematics(caseConfig);

%% 视频导出
for videoTypeIdx = 1:numel(caseConfig.videoTypes)
    videoType = caseConfig.videoTypes{videoTypeIdx};

    switch videoType
        case 'trajectory'
            renderBlockTrajectoryVideo(result, blockBody, caseConfig, videoDir, bodyRadius);
        case 'ghost'
            renderBlockGhostVideo(result, blockBody, caseConfig, videoDir, bodyRadius);
        otherwise
            error('Unsupported block video type: %s', videoType);
    end
end
