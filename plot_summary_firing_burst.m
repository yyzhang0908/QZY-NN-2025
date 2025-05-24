function plot_summary_firing_burst(idx1, idx0, burst_rg, fr_rg)
burst1 = idx1.spike_burst_perc;
burst0 = idx0.spike_burst_perc;
fr1 = idx1.spike_fr;
fr0 = idx0.spike_fr;
if nargin < 3
    burst_rg = [0, max(max(burst0), max(burst1))+0.1];
end
if nargin < 4
    fr_rg = [0, max(max(fr0), max(fr1))+0.1];
end

subplot('Position', [0.12, 0.54, 0.35, 0.35]);
plot(burst0, burst1, 'o', 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'w');
set(gca, 'TickDir', 'out', 'YLim', burst_rg, 'XLim', burst_rg);
xlabel('Inter-burst burst percentage (%)');
ylabel('Delta-burst burst percentage (%)');
title('Burst percentage');

subplot('Position', [0.60, 0.54, 0.35, 0.35]);
plot(fr1, burst1, 'o', 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'w');
set(gca, 'TickDir', 'out', 'YLim', burst_rg, 'XLim', fr_rg);
xlabel('Delta-burst firing rate (Hz)');
ylabel('Delta-burst burst percentage (%)');
title('Delta-burst firing rate - burst percentage');

subplot('Position', [0.12, 0.07, 0.35, 0.35]);
plot(burst0, fr0, 'o', 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'w');
set(gca, 'TickDir', 'out', 'YLim', fr_rg, 'XLim', burst_rg);
xlabel('Inter-burst burst percentage (%)');
ylabel('Inter-burst firing rate (Hz)');
title('Inter-burst firing rate - burst percentage');

subplot('Position', [0.60, 0.07, 0.35, 0.35]);
plot(fr1, fr0, 'o', 'MarkerFaceColor', [.5, .5, .5], 'MarkerEdgeColor', 'w');
set(gca, 'TickDir', 'out', 'YLim', fr_rg, 'XLim', fr_rg);
xlabel('Delta-burst firing rate (Hz)');
ylabel('Inter-burst firing rate (Hz)');
title('Firing rate');

end