function result = simulateBlockSpatialCase(config)
% 输入使用空间速度，内部按质心平动与转动递推长方体运动。

body = config.body;
time = 0:config.dt:config.tEnd;
stepCount = numel(time);
gravity = config.gravity(:);

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

state = hw21utils.createComStateFromSpatialInit(config, body);

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
    diagnostics = hw21utils.computeComStateDiagnostics(state, body);

    originHistory(:, idx) = diagnostics.originPosition;
    rComHistory(:, idx) = diagnostics.rCom;
    vOriginHistory(:, idx) = diagnostics.vOrigin;
    vComHistory(:, idx) = diagnostics.vCom;
    vSpatialHistory(:, idx) = diagnostics.vSpatial;
    wWorldHistory(:, idx) = state.wWorld;
    rotationHistory(:, :, idx) = state.R;
    linearMomentumHistory(:, idx) = diagnostics.linearMomentum;
    angularMomentumHistory(:, idx) = diagnostics.angularMomentum;
    if trackedPointCount > 0
        trackedPointHistory(:, idx, :) = reshape( ...
            hw21utils.transformBodyPoints(state.R, diagnostics.originPosition, trackedPointsBody), ...
            3, 1, trackedPointCount ...
        );
    end
    kineticTransHistory(idx) = 0.5 * body.mass * dot(diagnostics.vCom, diagnostics.vCom);
    kineticRotHistory(idx) = 0.5 * state.wWorld.' * diagnostics.worldInertiaCom * state.wWorld;
    potentialHistory(idx) = -body.mass * dot(gravity, diagnostics.rCom);

    if idx == stepCount
        continue;
    end

    forceWorld = body.mass * gravity;
    torqueWorld = cross(diagnostics.rCom, forceWorld);
    state = hw21utils.advanceComState(state, body, forceWorld, torqueWorld, config.dt);
end

result = struct();
result.name = config.name;
result.body = body;
result.gravity = gravity;
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
