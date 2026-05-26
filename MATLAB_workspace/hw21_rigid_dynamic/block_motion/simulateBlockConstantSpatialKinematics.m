function result = simulateBlockConstantSpatialKinematics(config)
% 在恒定空间速度下生成长方体的运动学轨迹。

body = config.body;
time = 0:config.dt:config.tEnd;
stepCount = numel(time);

trackedPointsBody = [];
trackedPointNames = {};

if isfield(config, 'trackedPointsBody') && ~isempty(config.trackedPointsBody)
    trackedPointsBody = config.trackedPointsBody;
    trackedPointCount = size(trackedPointsBody, 2);
    trackedPointHistory = zeros(3, stepCount, trackedPointCount);

    if isfield(config, 'trackedPointNames') && ~isempty(config.trackedPointNames)
        trackedPointNames = config.trackedPointNames;
    else
        trackedPointNames = arrayfun(@(idx) sprintf('Point %d', idx), 1:trackedPointCount, 'UniformOutput', false);
    end
else
    trackedPointCount = 0;
    trackedPointHistory = zeros(3, stepCount, 0);
end

originPosition = config.origin0(:);
rotationMatrix = config.R0;
vSpatial = config.vSpatial0(:);
wWorld = config.wWorld0(:);

originHistory = zeros(3, stepCount);
rComHistory = zeros(3, stepCount);
vOriginHistory = zeros(3, stepCount);
vComHistory = zeros(3, stepCount);
vSpatialHistory = zeros(3, stepCount);
wWorldHistory = zeros(3, stepCount);
angularMomentumHistory = zeros(3, stepCount);
linearMomentumHistory = zeros(3, stepCount);
rotationHistory = zeros(3, 3, stepCount);
kineticTransHistory = zeros(1, stepCount);
kineticRotHistory = zeros(1, stepCount);
potentialHistory = zeros(1, stepCount);

for idx = 1:stepCount
    rCom = originPosition + rotationMatrix * body.comBody;
    vCom = vSpatial + cross(wWorld, rCom);
    state = struct( ...
        'rCom', rCom, ...
        'vCom', vCom, ...
        'R', rotationMatrix, ...
        'wWorld', wWorld);
    diagnostics = hw21utils.computeComStateDiagnostics(state, body);

    originHistory(:, idx) = diagnostics.originPosition;
    rComHistory(:, idx) = diagnostics.rCom;
    vOriginHistory(:, idx) = diagnostics.vOrigin;
    vComHistory(:, idx) = diagnostics.vCom;
    vSpatialHistory(:, idx) = diagnostics.vSpatial;
    wWorldHistory(:, idx) = wWorld;
    rotationHistory(:, :, idx) = rotationMatrix;
    linearMomentumHistory(:, idx) = diagnostics.linearMomentum;
    angularMomentumHistory(:, idx) = diagnostics.angularMomentum;
    if trackedPointCount > 0
        trackedPointHistory(:, idx, :) = reshape( ...
            hw21utils.transformBodyPoints(rotationMatrix, diagnostics.originPosition, trackedPointsBody), ...
            3, 1, trackedPointCount ...
        );
    end
    kineticTransHistory(idx) = 0.5 * body.mass * dot(diagnostics.vCom, diagnostics.vCom);
    kineticRotHistory(idx) = 0.5 * wWorld.' * diagnostics.worldInertiaCom * wWorld;

    if idx == stepCount
        continue;
    end

    [originPosition, rotationMatrix] = advanceConstantSpatialPose( ...
        originPosition, rotationMatrix, vSpatial, wWorld, config.dt);
end

result = struct();
result.name = config.name;
result.body = body;
result.gravity = [0; 0; 0];
result.time = time;
result.originHistory = originHistory;
result.rComHistory = rComHistory;
result.vOriginHistory = vOriginHistory;
result.vComHistory = vComHistory;
result.vSpatialHistory = vSpatialHistory;
result.wWorldHistory = wWorldHistory;
result.linearMomentumHistory = linearMomentumHistory;
result.angularMomentumHistory = angularMomentumHistory;
result.rotationHistory = rotationHistory;
result.kineticTransHistory = kineticTransHistory;
result.kineticRotHistory = kineticRotHistory;
result.potentialHistory = potentialHistory;
result.trackedPointsBody = trackedPointsBody;
result.trackedPointNames = trackedPointNames;
result.trackedPointHistory = trackedPointHistory;
result.velocityDescription = 'spatial';
end

function [originNext, rotationNext] = advanceConstantSpatialPose(originPosition, rotationMatrix, vSpatial, wWorld, dt)
% 对恒定空间速度进行一步精确位姿更新。

omegaNorm = norm(wWorld);
skewOmega = hw21utils.skew3(wWorld);

if omegaNorm < 1.0e-12
    deltaRotation = eye(3);
    deltaTranslation = dt * vSpatial;
else
    angle = omegaNorm * dt;
    skewOmega2 = skewOmega * skewOmega;
    sinOverOmega = sin(angle) / omegaNorm;
    oneMinusCosOverOmega2 = (1 - cos(angle)) / (omegaNorm^2);
    dtMinusSinOverOmegaOverOmega2 = (dt - sin(angle) / omegaNorm) / (omegaNorm^2);

    deltaRotation = eye(3) + ...
        sinOverOmega * skewOmega + ...
        oneMinusCosOverOmega2 * skewOmega2;
    deltaTranslation = ( ...
        dt * eye(3) + ...
        oneMinusCosOverOmega2 * skewOmega + ...
        dtMinusSinOverOmegaOverOmega2 * skewOmega2 ...
    ) * vSpatial;
end

originNext = deltaRotation * originPosition + deltaTranslation;
rotationNext = hw21utils.projectRotation(deltaRotation * rotationMatrix);
end
