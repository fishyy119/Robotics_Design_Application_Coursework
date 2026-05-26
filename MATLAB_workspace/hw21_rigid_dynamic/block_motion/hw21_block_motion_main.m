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

%% 算例配置
caseConfigs = getBlockCaseConfigs(bodyRadius);

%% 数值积分与视频导出
for caseIdx = 1:numel(caseConfigs)
    caseConfig = caseConfigs(caseIdx);
    caseConfig.body = blockBody;
    caseConfig.trackedPointsBody = blockBody.referenceCornerBody;
    caseConfig.trackedPointNames = {'Upper corner'};

    result = simulateBlockSpatialCase(caseConfig);

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
end
