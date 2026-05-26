function [forceWorld, torqueWorld, contactInfo] = evaluateTopContact(state, body, config)
% 计算质心平动与转动状态下陀螺的重力和接触受力。

gravityForce = body.mass * config.gravity(:);
originPosition = state.rCom - state.R * body.comBody;
tipPosition = originPosition + state.R * body.tipBody;
tipVelocity = state.vCom + cross(state.wWorld, tipPosition - state.rCom);

penetration = max(0.0, -tipPosition(3));
contactForce = zeros(3, 1);
normalForce = 0.0;

if penetration > 0
    contactForce(1:2) = -config.contactTangentialDamping * tipVelocity(1:2);
    normalForce = config.contactStiffness * penetration - config.contactNormalDamping * tipVelocity(3);
    normalForce = max(normalForce, 0.0);
    contactForce(3) = normalForce;
end

forceWorld = gravityForce + contactForce;
torqueWorld = cross(state.rCom, gravityForce) + cross(tipPosition, contactForce);

contactInfo = struct();
contactInfo.rCom = state.rCom;
contactInfo.tipPosition = tipPosition;
contactInfo.tipVelocity = tipVelocity;
contactInfo.contactForce = contactForce;
contactInfo.penetration = penetration;
contactInfo.normalForce = normalForce;
end
