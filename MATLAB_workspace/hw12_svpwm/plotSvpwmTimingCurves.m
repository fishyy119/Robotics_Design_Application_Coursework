function fig = plotSvpwmTimingCurves(data)
% 根据预计算数据绘制 SVPWM 各段归一化时间曲线。

fig = utils.createFigureA4(struct('Name', 'SVPWM Timing Curves', 'Width', 16, 'AspectRatio', 0.92));
layout = tiledlayout(fig, 6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

timeDataList = {data.tau0, data.tau1, data.tau2, data.tau7};
timeLabelList = {'$T_0 / T_s$', '$T_1 / T_s$', '$T_2 / T_s$', '$T_7 / T_s$'};
timeColorList = {
    utils.tab10Color(1), ...
    utils.tab10Color(2), ...
    utils.tab10Color(3), ...
    utils.tab10Color(4)
};
timeLineStyleList = {'-', '-', '-', '-'};
timeUpper = max([data.tau0(:); data.tau1(:); data.tau2(:); data.tau7(:)]);

timeAxes = gobjects(4, 1);
for idx = 1:4
    timeAxes(idx) = nexttile(layout);
    hold(timeAxes(idx), 'on');
    grid(timeAxes(idx), 'on');
    box(timeAxes(idx), 'on');

    plot(timeAxes(idx), data.angleDeg, timeDataList{idx}, ...
        'Color', timeColorList{idx}, ...
        'LineStyle', timeLineStyleList{idx});

    ylabel(timeAxes(idx), timeLabelList{idx});
    ylim(timeAxes(idx), [0, 1.05 * timeUpper]);
    xlim(timeAxes(idx), [data.angleDeg(1), data.angleDeg(end)]);

    if idx == 1
        title(timeAxes(idx), sprintf(['$U_{dc} = %.0f\\,\\mathrm{V}$, ', ...
            '$u_d = %.0f\\,\\mathrm{V}$, $u_q = %.0f\\,\\mathrm{V}$'], data.dcBus, data.ud, data.uq));
    end

    timeAxes(idx).XTickLabel = [];
end

ax5 = nexttile(layout, [2, 1]);
hold(ax5, 'on');
grid(ax5, 'on');
box(ax5, 'on');

yyaxis(ax5, 'left');
hSector = stairs(ax5, data.angleDeg, data.sector, 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.1);
ylabel(ax5, 'Sector');
ylim(ax5, [0.5, 6.5]);
yticks(ax5, 1:6);

yyaxis(ax5, 'right');
hAlpha = plot(ax5, data.angleDeg, data.uAlphaFromTimings, 'Color', utils.tab10Color(5), 'DisplayName', '$u_{\alpha}$');
hBeta = plot(ax5, data.angleDeg, data.uBetaFromTimings, 'Color', utils.tab10Color(6), 'DisplayName', '$u_{\beta}$');
ylabel(ax5, 'Voltage (V)');

xlabel(ax5, 'Electrical angle (deg)');
lgd = legend(ax5, [hSector, hAlpha, hBeta], {'Sector', '$u_{\alpha}$', '$u_{\beta}$'}, ...
    'Location', 'eastoutside', 'Orientation', 'vertical');
lgd.Interpreter = 'latex';
xlim(ax5, [data.angleDeg(1), data.angleDeg(end)]);

linkaxes([timeAxes; ax5], 'x');
end
