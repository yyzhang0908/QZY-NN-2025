function plot_summary_lfp_eeg_shift(dts, dphis, ttl, clr)
time_rg = [-0.2, 0.2];
phase_rg = [-pi, pi];
if nargin < 3
    ttl = 'LFP refer to EEG delta peak';
end
if nargin < 4
    clr = 'g';
end

ns = length(dts);
cmean = mean([dphis, dts], 1, 'omitnan');
csem = std([dphis, dts], 1, 1, 'omitnan')/sqrt(ns);

hold on;

for jj = 1:ns
    plot(dphis(jj), dts(jj), 'o', 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', clr, 'MarkerSize', 6);
end
plot([cmean(1)-csem(1), cmean(1)+csem(1)], [cmean(2), cmean(2)], '-k');
plot([cmean(1), cmean(1)], [cmean(2)-csem(2), cmean(2)+csem(2)], '-k');
plot(cmean(1), cmean(2), 'or', 'MarkerFaceColor', clr, ...
    'MarkerEdgeColor', 'w', 'MarkerSize', 10);

set(gca, 'TickDir', 'out', 'YLim', time_rg, 'XLim', phase_rg, 'XTick', ...
    (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
xlabel('Shift phase');
ylabel('Shift time (s)');
title(sprintf('%s\nPhase: %.2f (%.0f\\circ), Time: %.2fs', ...
    ttl, cmean(1), rad2deg(cmean(1)), cmean(2)));

hold off;
end