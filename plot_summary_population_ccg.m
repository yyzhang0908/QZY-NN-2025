function plot_summary_population_ccg(lag, cbsl, cburst, citr, cemer, ttl)
if nargin < 6
    ttl = 'Cross-correlogram of population';
end
ns = size(cbsl, 1);
hold on;
shadedErrorBar(lag, mean(cbsl, 1, 'omitnan'), std(cbsl, 1, 1, 'omitnan')/sqrt(ns), ...
    'lineProps', {'-', 'Color', 'k', 'LineWidth', 1.5});
shadedErrorBar(lag, mean(citr, 1, 'omitnan'), std(citr, 1, 1, 'omitnan')/sqrt(ns), ...
    'lineProps', {'-', 'Color', 'b', 'LineWidth', 1.5});
shadedErrorBar(lag, mean(cburst, 1, 'omitnan'), std(cburst, 1, 1, 'omitnan')/sqrt(ns), ...
    'lineProps', {'-', 'Color', 'r', 'LineWidth', 1.5});
shadedErrorBar(lag, mean(cemer, 1, 'omitnan'), std(cemer, 1, 1, 'omitnan')/sqrt(ns), ...
    'lineProps', {'-', 'Color', 'm', 'LineWidth', 1.5});
xline(0, ':k');
set(gca, 'TickDir', 'out', 'XLim', [min(lag), max(lag)]);
xlabel('Lag (ms)');
ylabel('Z-scored correlation');
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off', ...
    'Location', 'best');
title(sprintf('%s', ttl));
hold off;
end
