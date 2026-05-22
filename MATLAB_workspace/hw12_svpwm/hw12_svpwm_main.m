clearvars;
close all;
clc;
utils.setDefaultGraphics;

%% 参数
dcBus = 48;
ud = 0;
uq = 24;
Ts = 1;
angleDeg = linspace(0, 1000, 4001);

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
outputDir = fullfile(projectRoot, 'LaTeX_workspace', 'figures', 'hw12');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% 数据计算
thetaElec = deg2rad(angleDeg);
uAlpha = ud .* cos(thetaElec) - uq .* sin(thetaElec);
uBeta = ud .* sin(thetaElec) + uq .* cos(thetaElec);

%% SVPWM 时序计算
theta = mod(atan2(uBeta, uAlpha), 2 * pi);
sector = floor(theta ./ (pi / 3)) + 1;
sector(sector == 7) = 6;

T1 = zeros(size(uAlpha));
T2 = zeros(size(uAlpha));

mask = sector == 1;
T1(mask) = Ts / dcBus * (1.5 * uAlpha(mask) - sqrt(3) / 2 * uBeta(mask));
T2(mask) = Ts / dcBus * (sqrt(3) * uBeta(mask));

mask = sector == 2;
T1(mask) = Ts / dcBus * (1.5 * uAlpha(mask) + sqrt(3) / 2 * uBeta(mask));
T2(mask) = Ts / dcBus * (-1.5 * uAlpha(mask) + sqrt(3) / 2 * uBeta(mask));

mask = sector == 3;
T1(mask) = Ts / dcBus * (sqrt(3) * uBeta(mask));
T2(mask) = Ts / dcBus * (-1.5 * uAlpha(mask) - sqrt(3) / 2 * uBeta(mask));

mask = sector == 4;
T1(mask) = Ts / dcBus * (-1.5 * uAlpha(mask) + sqrt(3) / 2 * uBeta(mask));
T2(mask) = Ts / dcBus * (-sqrt(3) * uBeta(mask));

mask = sector == 5;
T1(mask) = Ts / dcBus * (-1.5 * uAlpha(mask) - sqrt(3) / 2 * uBeta(mask));
T2(mask) = Ts / dcBus * (1.5 * uAlpha(mask) - sqrt(3) / 2 * uBeta(mask));

mask = sector == 6;
T1(mask) = Ts / dcBus * (-sqrt(3) * uBeta(mask));
T2(mask) = Ts / dcBus * (1.5 * uAlpha(mask) + sqrt(3) / 2 * uBeta(mask));

T0 = 0.5 * (Ts - T1 - T2);
T7 = T0;

if any(T0(:) < -1e-12)
    error('Reference voltage exceeds the linear SVPWM range.');
end

T0 = max(T0, 0);
T1 = max(T1, 0);
T2 = max(T2, 0);
T7 = max(T7, 0);

%% 反算电压分量
thetaStart = (sector - 1) * (pi / 3);
activeVoltage = 2 / 3 * dcBus;

uAlphaFromTimings = activeVoltage / Ts * ...
    (T1 .* cos(thetaStart) + T2 .* cos(thetaStart + pi / 3));
uBetaFromTimings = activeVoltage / Ts * ...
    (T1 .* sin(thetaStart) + T2 .* sin(thetaStart + pi / 3));

data = struct( ...
    'angleDeg', angleDeg, ...
    'ud', ud, ...
    'uq', uq, ...
    'dcBus', dcBus, ...
    'Ts', Ts, ...
    'uAlpha', uAlpha, ...
    'uBeta', uBeta, ...
    'uAlphaFromTimings', uAlphaFromTimings, ...
    'uBetaFromTimings', uBetaFromTimings, ...
    'T0', T0, ...
    'T1', T1, ...
    'T2', T2, ...
    'T7', T7, ...
    'tau0', T0 ./ Ts, ...
    'tau1', T1 ./ Ts, ...
    'tau2', T2 ./ Ts, ...
    'tau7', T7 ./ Ts, ...
    'sector', sector ...
);

%% 绘图
utils.setDefaultGraphics();

sectorFig = plotSvpwmSectorMap(dcBus);
exportgraphics(sectorFig, fullfile(outputDir, 'svpwm_sector_distribution.pdf'), 'ContentType', 'vector');

timingFig = plotSvpwmTimingCurves(data);
exportgraphics(timingFig, fullfile(outputDir, 'svpwm_timing_curves.pdf'), 'ContentType', 'vector');
