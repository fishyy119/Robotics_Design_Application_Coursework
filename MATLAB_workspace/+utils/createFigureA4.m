function fig = createFigureA4(params)
% CREATEFIGUREA4 使用结构体参数生成论文尺寸 figure
%
% params 是一个结构体，可包含字段：
%   params.Name        图窗名称，默认 ""
%   params.Width       图宽 cm，默认 18
%   params.AspectRatio 高/宽比，默认 2/3

% 默认值
defaults.Name = "";
defaults.Width = 9;
defaults.AspectRatio = 2/3;

% 合并结构体，缺省字段使用默认值
if nargin < 1
    params = struct();
end

if ~isfield(params,'Name') || isempty(params.Name)
    params.Name = defaults.Name;
end
if ~isfield(params,'Width') || isempty(params.Width)
    params.Width = defaults.Width;
end
if ~isfield(params,'AspectRatio') || isempty(params.AspectRatio)
    params.AspectRatio = defaults.AspectRatio;
end

Height = params.Width * params.AspectRatio;

fig = figure( ...
    'Name', params.Name, ...
    'NumberTitle', 'off', ...
    'Units', 'centimeters', ...
    'Position', [5,5,params.Width,Height], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0,0,params.Width,Height] ...
);

end