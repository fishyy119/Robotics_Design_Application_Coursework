function diagnostics = computeComStateDiagnostics(state, body)
% 根据质心平动与转动状态恢复空间速度相关诊断量。

if isfield(body, 'comBody') && ~isempty(body.comBody)
    comBody = body.comBody(:);
else
    comBody = [0; 0; 0];
end

originPosition = state.rCom - state.R * comBody;
vSpatial = state.vCom - cross(state.wWorld, state.rCom);
vOrigin = state.vCom - cross(state.wWorld, state.R * comBody);
worldInertiaCom = state.R * body.inertiaBody * state.R.';
linearMomentum = body.mass * state.vCom;
angularMomentum = cross(state.rCom, linearMomentum) + worldInertiaCom * state.wWorld;

diagnostics = struct();
diagnostics.originPosition = originPosition;
diagnostics.rCom = state.rCom;
diagnostics.vOrigin = vOrigin;
diagnostics.vCom = state.vCom;
diagnostics.vSpatial = vSpatial;
diagnostics.linearMomentum = linearMomentum;
diagnostics.angularMomentum = angularMomentum;
diagnostics.worldInertiaCom = worldInertiaCom;
end
