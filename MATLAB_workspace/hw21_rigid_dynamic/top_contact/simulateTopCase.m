function result = simulateTopCase(config)
% 输入使用空间速度，内部按质心平动与转动递推陀螺运动。

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

rOriginHistory = zeros(3, stepCount);
rComHistory = zeros(3, stepCount);
vOriginHistory = zeros(3, stepCount);
vComHistory = zeros(3, stepCount);
wWorldHistory = zeros(3, stepCount);
rotationHistory = zeros(3, 3, stepCount);
tipPositionHistory = zeros(3, stepCount);
tipVelocityHistory = zeros(3, stepCount);
contactForceHistory = zeros(3, stepCount);
totalForceHistory = zeros(3, stepCount);
penetrationHistory = zeros(1, stepCount);
normalForceHistory = zeros(1, stepCount);
kineticTransHistory = zeros(1, stepCount);
kineticRotHistory = zeros(1, stepCount);
potentialHistory = zeros(1, stepCount);
linearMomentumHistory = zeros(3, stepCount);
angularMomentumHistory = zeros(3, stepCount);
vSpatialHistory = zeros(3, stepCount);

for idx = 1:stepCount
    [forceWorld, torqueWorld, contactInfo] = evaluateTopContact(state, body, config);
    diagnostics = hw21utils.computeComStateDiagnostics(state, body);

    rOriginHistory(:, idx) = diagnostics.originPosition;
    rComHistory(:, idx) = diagnostics.rCom;
    vOriginHistory(:, idx) = diagnostics.vOrigin;
    vComHistory(:, idx) = diagnostics.vCom;
    vSpatialHistory(:, idx) = diagnostics.vSpatial;
    wWorldHistory(:, idx) = state.wWorld;
    rotationHistory(:, :, idx) = state.R;
    if trackedPointCount > 0
        trackedPointHistory(:, idx, :) = reshape( ...
            hw21utils.transformBodyPoints(state.R, diagnostics.originPosition, trackedPointsBody), ...
            3, 1, trackedPointCount ...
        );
    end
    tipPositionHistory(:, idx) = contactInfo.tipPosition;
    tipVelocityHistory(:, idx) = contactInfo.tipVelocity;
    contactForceHistory(:, idx) = contactInfo.contactForce;
    totalForceHistory(:, idx) = forceWorld;
    penetrationHistory(idx) = contactInfo.penetration;
    normalForceHistory(idx) = contactInfo.normalForce;
    kineticTransHistory(idx) = 0.5 * body.mass * dot(diagnostics.vCom, diagnostics.vCom);
    kineticRotHistory(idx) = 0.5 * state.wWorld.' * diagnostics.worldInertiaCom * state.wWorld;
    potentialHistory(idx) = -body.mass * dot(gravity, diagnostics.rCom);
    linearMomentumHistory(:, idx) = diagnostics.linearMomentum;
    angularMomentumHistory(:, idx) = diagnostics.angularMomentum;

    if idx == stepCount
        continue;
    end

    state = hw21utils.advanceComState(state, body, forceWorld, torqueWorld, config.dt);
end

result = struct();
result.name = config.name;
result.body = body;
result.gravity = gravity;
result.time = time;
result.originHistory = rOriginHistory;
result.rComHistory = rComHistory;
result.vOriginHistory = vOriginHistory;
result.vComHistory = vComHistory;
result.vSpatialHistory = vSpatialHistory;
result.wWorldHistory = wWorldHistory;
result.rotationHistory = rotationHistory;
result.tipPositionHistory = tipPositionHistory;
result.tipVelocityHistory = tipVelocityHistory;
result.contactForceHistory = contactForceHistory;
result.totalForceHistory = totalForceHistory;
result.penetrationHistory = penetrationHistory;
result.normalForceHistory = normalForceHistory;
result.kineticTransHistory = kineticTransHistory;
result.kineticRotHistory = kineticRotHistory;
result.potentialHistory = potentialHistory;
result.linearMomentumHistory = linearMomentumHistory;
result.angularMomentumHistory = angularMomentumHistory;
result.trackedPointsBody = trackedPointsBody;
result.trackedPointNames = trackedPointNames;
result.trackedPointHistory = trackedPointHistory;
result.velocityDescription = 'spatial';
end
