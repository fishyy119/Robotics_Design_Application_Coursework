function fig = plotSvpwmSectorMap(dcBus)
% 绘制 SVPWM 六扇区与基本矢量分布图。

if nargin < 1 || isempty(dcBus)
    dcBus = 48;
end

activeRadius = 2 / 3 * dcBus;
vertexAngles = deg2rad(0:60:360);
sectorAngles = deg2rad(30:60:330);

xVertex = activeRadius * cos(vertexAngles);
yVertex = activeRadius * sin(vertexAngles);

fig = utils.createFigureA4(struct('Name', 'SVPWM Sector Distribution', 'Width', 13.5, 'AspectRatio', 0.78));
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');
box(ax, 'on');

sectorLabels = {'I', 'II', 'III', 'IV', 'V', 'VI'};

for k = 1:6
    sectorColor = utils.tab10Color(k);
    patch(ax, ...
        [0, xVertex(k), xVertex(k + 1)], ...
        [0, yVertex(k), yVertex(k + 1)], ...
        sectorColor, ...
        'FaceAlpha', 0.22, ...
        'EdgeColor', 'none');
end

plot(ax, xVertex, yVertex, 'Color', [0.15, 0.15, 0.15], 'LineWidth', 1.4);

for k = 1:6
    plot(ax, [0, activeRadius * cos(vertexAngles(k))], [0, activeRadius * sin(vertexAngles(k))], ...
        '--', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 0.8);
    quiver(ax, 0, 0, activeRadius * cos(vertexAngles(k)), activeRadius * sin(vertexAngles(k)), ...
        0, 'Color', [0, 0, 0], 'LineWidth', 1.0, 'MaxHeadSize', 0.14);
    text(1.08 * activeRadius * cos(vertexAngles(k)), 1.08 * activeRadius * sin(vertexAngles(k)), ...
        sprintf('$V_{%d}$', k), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    text(0.56 * activeRadius * cos(sectorAngles(k)), 0.56 * activeRadius * sin(sectorAngles(k)), ...
        sprintf('Sector %s', sectorLabels{k}), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold');
end

plot(ax, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4.5);
text(0.06 * activeRadius, 0.05 * activeRadius, 'O', ...
    'Interpreter', 'latex', ...
    'HorizontalAlignment', 'left');

xlabel(ax, '$u_{\alpha}$ (V)');
ylabel(ax, '$u_{\beta}$ (V)');
title(ax, sprintf('$U_{dc} = %.0f\\,\\mathrm{V}$', dcBus));

limit = 1.22 * activeRadius;
xlim(ax, [-limit, limit]);
ylim(ax, [-limit, limit]);
end
