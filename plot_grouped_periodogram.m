function plot_grouped_periodogram(pxxs, fxx, pvals, alpha, ...
    xlim, xticks, xlog, ylb, lgd, ttl, clr)
hold on;

ngrp = numel(pxxs);
for kk = 1:ngrp
    shadedErrorBar(fxx, mean(pxxs{kk}, 1, 'omitnan'), ...
        std(pxxs{kk}, 1, 1, 'omitnan') / sqrt(size(pxxs{kk}, 1)), ...
        'lineProps', {'-', 'Color', clr{kk}, 'LineWidth', 1});
end

yunit = range(get(gca, 'YLim')) * 0.02;
for kk = 1:ngrp
    cp = double(pvals(kk, :) < alpha);
    cp(cp == 0) = nan;
    plot(fxx, -kk * yunit * cp, '.', 'Color', clr{kk}, 'LineWidth', 1.5);
end

if xlog; set(gca, 'XScale', 'log'); end
set(gca, 'XLim', xlim, 'XTick', xticks, 'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel(ylb);
legend(lgd, 'Box', 'off', 'Location', 'best');
title(sprintf('%s', ttl));
hold off;
end