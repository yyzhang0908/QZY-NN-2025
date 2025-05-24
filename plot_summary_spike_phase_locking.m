function plot_summary_spike_phase_locking(rzpmks, ttl, clr)
if nargin < 2
    ttl = 'spike-LFP phase locking';
end
if nargin < 3
    clr = 'k';
end

ns = size(rzpmks, 1);
[R, Z, P, M, K] = circular_statistics(rzpmks(:, 4));
V = rad2deg(sqrt(-2*log(R-1e-8))/sqrt(ns));
for ii = 1:ns
    if rzpmks(ii, 3)<0.05
        polarplot(rzpmks(ii, 4), rzpmks(ii, 5), 'o', ...
            'MarkerFaceColor', clr, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    else
        polarplot(rzpmks(ii, 4), rzpmks(ii, 5), 'o', ...
            'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    end
    hold on;
end
set(gca, 'RLim', [0, 6], 'ThetaZeroLocation', 'top', 'ThetaDir', ...
    'clockwise', 'ThetaAxisUnits', 'degrees', 'ThetaLim', [-180, 180]);
set(get(get(gca, 'RAxis'), 'Label'), 'String', '\kappa', ...
    'Position', [0.3, 3.05, 0], 'Rotation', 0);
title(sprintf('%s\nR=%.4f, Z=%.2f, \\itp\\rm=%.3f, \\mu=%.0f\\circ (\\pm%.0f\\circ), \\kappa=%.3f', ...
    ttl, R, Z, P, rad2deg(M), V, K));

hold off;
end