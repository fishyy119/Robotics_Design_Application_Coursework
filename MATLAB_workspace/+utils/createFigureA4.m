function fig = createFigureA4(params)
% 使用结构体参数创建论文尺寸的 figure。
%
% params 可包含以下字段：
%   params.Name        图窗名称，默认 ""
%   params.Width       图宽，单位 cm，默认 9
%   params.AspectRatio 高宽比，默认 2/3

% 设置默认值。
defaults.Name = "";
defaults.Width = 9;
defaults.AspectRatio = 2/3;

if nargin < 1
    params = struct();
end

% 补齐缺省字段。
if ~isfield(params, 'Name') || isempty(params.Name)
    params.Name = defaults.Name;
end
if ~isfield(params, 'Width') || isempty(params.Width)
    params.Width = defaults.Width;
end
if ~isfield(params, 'AspectRatio') || isempty(params.AspectRatio)
    params.AspectRatio = defaults.AspectRatio;
end

height_cm = params.Width * params.AspectRatio;

fig = figure( ...
    'Name', params.Name, ...
    'NumberTitle', 'off', ...
    'Units', 'centimeters', ...
    'Position', [5, 5, params.Width, height_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0, 0, params.Width, height_cm] ...
);

end
