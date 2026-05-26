function nextState = advanceComState(state, body, forceWorld, torqueWorld, dt)
% 用质心平动与转动状态推进刚体一步。

diagnostics = hw21utils.computeComStateDiagnostics(state, body);
torqueCom = torqueWorld - cross(diagnostics.rCom, forceWorld);
dwWorld = diagnostics.worldInertiaCom \ (torqueCom - cross(state.wWorld, diagnostics.worldInertiaCom * state.wWorld));
aCom = forceWorld / body.mass;

nextState = struct();
nextState.rCom = state.rCom + dt * state.vCom + 0.5 * dt^2 * aCom;
nextState.vCom = state.vCom + dt * aCom;
nextState.R = localIntegrateRotation(state.R, state.wWorld, dt);
nextState.wWorld = state.wWorld + dt * dwWorld;
end


function RNext = localIntegrateRotation(R, wWorld, dt)
% 按当前角速度推进一步姿态。

normW = norm(wWorld);

if normW < 1.0e-12
    RNext = R;
    return;
end

theta = normW * dt;
wUnit = wWorld / normW;
wHat = hw21utils.skew3(wUnit);
rot = eye(3) + sin(theta) * wHat + (1 - cos(theta)) * (wHat * wHat);
RNext = hw21utils.projectRotation(rot * R);
end
