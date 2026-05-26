function state = createComStateFromSpatialInit(config, body)
% 将空间速度输入转换为质心平动与转动状态。

if isfield(body, 'comBody') && ~isempty(body.comBody)
    comBody = body.comBody(:);
else
    comBody = [0; 0; 0];
end

state = struct();
state.rCom = config.origin0(:) + config.R0 * comBody;
state.vCom = config.vSpatial0(:) + cross(config.wWorld0(:), state.rCom);
state.R = config.R0;
state.wWorld = config.wWorld0(:);
end
