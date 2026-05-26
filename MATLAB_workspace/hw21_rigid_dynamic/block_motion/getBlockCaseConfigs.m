function caseConfigs = getBlockCaseConfigs(bodyRadius)
% 返回长方体仿真的算例与绘图配置。

caseConfigs = repmat(struct(), 2, 1);

caseConfigs(1).name = 'Block under Gravity';
caseConfigs(1).videoTypes = {'trajectory'};
caseConfigs(1).videoNames = struct( ...
    'trajectory', 'hw21_block_gravity.mp4');
caseConfigs(1).gravity = [0; 0; -9.81];
caseConfigs(1).origin0 = [0; 0; 0];
caseConfigs(1).vSpatial0 = [0; 0; 0];
caseConfigs(1).R0 = eye(3);
caseConfigs(1).wWorld0 = [1.0; 1.0; 1.0];
caseConfigs(1).dt = 0.005;
caseConfigs(1).tEnd = 5.0;
caseConfigs(1).viewAngle = [36, 20];
caseConfigs(1).axisPadding = 0.15;
caseConfigs(1).followZWindow = true;
caseConfigs(1).zWindowSpan = 2 * (bodyRadius + 0.20);
caseConfigs(1).ghostTimeStep = 0.04;

caseConfigs(2).name = 'Block in Zero Gravity';
caseConfigs(2).videoTypes = {'ghost'};
caseConfigs(2).videoNames = struct( ...
    'ghost', 'hw21_block_zero_gravity2.mp4');
caseConfigs(2).gravity = [0; 0; 0];
caseConfigs(2).origin0 = [0; 0; 0];
caseConfigs(2).vSpatial0 = [0.50; 0.10; 0.0];
caseConfigs(2).R0 = eye(3);
caseConfigs(2).wWorld0 = [1.0; 0.0; 1.0];
caseConfigs(2).dt = 0.01;
caseConfigs(2).tEnd = 5.0;
caseConfigs(2).viewAngle = [34, 18];
caseConfigs(2).axisPadding = 0.20;
caseConfigs(2).followZWindow = false;
caseConfigs(2).zWindowSpan = 0.0;
caseConfigs(2).ghostTimeStep = 0.3;


end
