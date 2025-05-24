function plot_summary_phase_locking_dist(rzpmks, ttl, clr)
deg_win = (-180:20:180);
prob_scale = 0.4;
theta_shift = -90;
theta_grid = (0.25:0.25:0.75);
rad_tick_grid = (-180+45:45:180);
arr_len = 1.4;
error_rad = 0.85;
tick_len = 0.05;
if nargin < 2
    ttl = 'spike-LFP phase locking';
end
if nargin < 3
    clr = [.5, .5, .5];
end

ns = size(rzpmks, 1);
[R, Z, P, M, K] = circular_statistics(rzpmks(:, 4));
V = rad2deg(sqrt(-2*log(R-1e-8))/sqrt(ns));
deg_win_cnt = (deg_win(1:end-1)+deg_win(2:end))/2;
deg_prob = histcounts(rad2deg(rzpmks(:, 4)), deg_win, 'Normalization', 'probability');

hold on;

[xv, yv] = deg_rho_to_xy((-180:3:180), 1, theta_shift);
plot(xv, yv, '-', 'Color', 'k', 'LineWidth', 1.3);
for kk = 1:length(theta_grid)
    [xv, yv] = deg_rho_to_xy((-180:3:180), theta_grid(kk), theta_shift);
    plot(xv, yv, '-', 'Color', [.7, .7, .7], 'LineWidth', 1);
end
for kk = 1:length(rad_tick_grid)
    [xv, yv] = deg_rho_to_xy(rad_tick_grid(kk), [0, 1], theta_shift);
    plot(xv, yv, '-', 'Color', [.7, .7, .7], 'LineWidth', 1);
    [xv, yv] = deg_rho_to_xy(rad_tick_grid(kk), [1-tick_len, 1], theta_shift);
    plot(xv, yv, '-', 'Color', 'k', 'LineWidth', 1.3);
    [xv, yv] = deg_rho_to_xy(rad_tick_grid(kk), 1.15, theta_shift);
    text(xv, yv, num2str(rad_tick_grid(kk)), 'HorizontalAlignment', 'center');
end

[xv, yv] = deg_rho_to_xy(deg_win_cnt, deg_prob/prob_scale, theta_shift);
patch(xv, yv, clr, 'EdgeColor', clr);
[xv, yv] = deg_rho_to_xy([0, -2, 0, 2]+rad2deg(M), [0, -0.1, -0.07, -0.1]+arr_len, ...
    theta_shift);
patch(xv, yv, 'r', 'EdgeColor', 'r');
[xv, yv] = deg_rho_to_xy(rad2deg(M), [0, arr_len], theta_shift);
plot(xv, yv, '-', 'Color', 'r', 'LineWidth', 1.3);
[xv, yv] = deg_rho_to_xy((-V:min(V/5, 3):V)+rad2deg(M), error_rad, theta_shift);
plot(xv, yv, '-', 'Color', 'r', 'LineWidth', 1.3);
[xv, yv] = deg_rho_to_xy(repmat([-V, V]+rad2deg(M), 2, 1), ...
    repmat([-1; 1]*tick_len/2+error_rad, 1, 2), theta_shift);
plot(xv, yv, '-', 'Color', 'r', 'LineWidth', 1.3);

axis equal; axis off;
set(gca, 'XLim', [-1, 1]*arr_len, 'YLim', [-1, 1]*arr_len);
title(sprintf('%s\nR=%.4f, Z=%.2f, \\itp\\rm=%.3f, \\mu=%.0f\\circ (\\pm%.0f\\circ), \\kappa=%.3f', ...
    ttl, R, Z, P, rad2deg(M), V, K));

end


function [x, y] = deg_rho_to_xy(deg, rho, dshift)
[x, y] = pol2cart(deg2rad(-deg-dshift), rho);
end




