function outpath = analyze_related_ca2(ca2_sig, ca2_ts, ca2_fs, opt, db_out, workpath)

t = tic;
fprintf('\n\tcalculating ... ');
load(db_out, 'eeg_sig', 'eeg_fs', 'eeg_ts', 'ev_mrk', 'delta_sig', ...
    'delta_burst_peak_t', 'delta_burst_peak_t_anes', ...
    'delta_burst_start', 'delta_burst_stop');
if ~isfolder(workpath); mkdir(workpath); end

%% compute metrics
% period info
t_anes = ev_mrk.t_anes;
t_bsl = ev_mrk.t_bsl;
t_emer = ev_mrk.t_emer;


% delta (inter) burst information
delta_burst_center = (delta_burst_start + delta_burst_stop)/2;
delta_burst_start_anes = mask_time_range(delta_burst_center, t_anes, delta_burst_start);
delta_burst_stop_anes = mask_time_range(delta_burst_center, t_anes, delta_burst_stop);
delta_burst_anes = [delta_burst_start_anes, delta_burst_stop_anes];
inter_burst_anes = [delta_burst_stop_anes(1:end-1), delta_burst_start_anes(2:end)];

% calcium signal spectrum
freq_limit = minmax(cat(2, opt.ca2_band_freq{:}));
[ca2_pxx, ca2_fxx] = wavelet_amplitude_spectrum(ca2_sig, ca2_ts, ca2_fs, ...
    freq_limit, opt.wavelet_scale);
ca2_pxx_burst = mean(interp1(ca2_ts, ca2_pxx, mean(delta_burst_anes, 2), ...
    'nearest'), 1, 'omitnan');
ca2_pxx_itr = mean(interp1(ca2_ts, ca2_pxx, mean(inter_burst_anes, 2), ...
    'nearest'), 1, 'omitnan');
ca2_pxx_bsl = mean(ca2_pxx(mask_time_range(ca2_ts, t_bsl), :), 1, 'omitnan');
ca2_pxx_anes = mean(ca2_pxx(mask_time_range(ca2_ts, t_anes), :), 1, 'omitnan');
ca2_pxx_emer = mean(ca2_pxx(mask_time_range(ca2_ts, t_emer), :), 1, 'omitnan');

ca2_pxx_ch_burst = ca2_pxx_burst./ca2_pxx_bsl-1;
ca2_pxx_ch_itr = ca2_pxx_itr./ca2_pxx_bsl-1;
[~, max_index_burst] = max(ca2_pxx_ch_burst);
ca2_cent_freq_burst = ca2_fxx(max_index_burst);
[~, max_index_itr] = max(ca2_pxx_ch_itr);
ca2_cent_freq_itr = ca2_fxx(max_index_itr);



ca2_pxx_per = ca2_pxx./ trapz(ca2_fxx, ca2_pxx, 2);
ca2_pxx_burst_per = mean(interp1(ca2_ts, ca2_pxx_per, mean(delta_burst_anes, 2), ...
    'nearest'), 1, 'omitnan');
ca2_pxx_itr_per = mean(interp1(ca2_ts, ca2_pxx_per, mean(inter_burst_anes, 2), ...
    'nearest'), 1, 'omitnan');
ca2_pxx_bsl_per = mean(ca2_pxx_per(mask_time_range(ca2_ts, t_bsl), :), 1, 'omitnan');
ca2_pxx_anes_per = mean(ca2_pxx_per(mask_time_range(ca2_ts, t_anes), :), 1, 'omitnan');
ca2_pxx_emer_per = mean(ca2_pxx_per(mask_time_range(ca2_ts, t_emer), :), 1, 'omitnan');


% calcium signal bandpower
ca2_bp_raw = integrate_psd_bandpower(ca2_pxx, ca2_fxx, opt.ca2_band_freq);
ca2_bp_raw_per = integrate_psd_bandpower(ca2_pxx_per, ca2_fxx, opt.ca2_band_freq);

[~, ca2_bp_bin, ca2_bp_bin_center] = get_event_timecourse( ...
    ca2_bp_raw, ca2_ts, opt.band_power_bin);
[~, ca2_bp_bin_per, ~] = get_event_timecourse( ...
    ca2_bp_raw_per, ca2_ts, opt.band_power_bin);

ca2_bp_bsl = integrate_psd_bandpower(ca2_pxx_bsl, ca2_fxx, opt.ca2_band_freq);
ca2_bp_anes = integrate_psd_bandpower(ca2_pxx_anes, ca2_fxx, opt.ca2_band_freq);
ca2_bp_burst = integrate_psd_bandpower(ca2_pxx_burst, ca2_fxx, opt.ca2_band_freq);
ca2_bp_itr = integrate_psd_bandpower(ca2_pxx_itr, ca2_fxx, opt.ca2_band_freq);
ca2_bp_emer = integrate_psd_bandpower(ca2_pxx_emer, ca2_fxx, opt.ca2_band_freq);

ca2_bp_bsl_per = integrate_psd_bandpower(ca2_pxx_bsl_per, ca2_fxx, opt.ca2_band_freq);
ca2_bp_anes_per = integrate_psd_bandpower(ca2_pxx_anes_per, ca2_fxx, opt.ca2_band_freq);
ca2_bp_burst_per = integrate_psd_bandpower(ca2_pxx_burst_per, ca2_fxx, opt.ca2_band_freq);
ca2_bp_itr_per = integrate_psd_bandpower(ca2_pxx_itr_per, ca2_fxx, opt.ca2_band_freq);
ca2_bp_emer_per = integrate_psd_bandpower(ca2_pxx_emer_per, ca2_fxx, opt.ca2_band_freq);

% calcium signal peak
is_delta_band = strcmpi(opt.ca2_band_label, 'delta');
[~, delta_ca2] = get_band_envelope( ...
    ca2_sig, ca2_fs, opt.ca2_band_freq(is_delta_band));
delta_burst_ca2_peak_dt = get_peri_event_dt(delta_ca2, ca2_ts, ...
    delta_burst_peak_t);
delta_burst_ca2_peak_dt_med = median(delta_burst_ca2_peak_dt);
delta_burst_ca2_peak_dt_anes = get_peri_event_dt(delta_ca2, ca2_ts, ...
    delta_burst_peak_t_anes);
delta_burst_ca2_peak_dt_anes_med = median(delta_burst_ca2_peak_dt_anes);

% delta burst related calcium timecourse
pets_delta_sig = get_peri_event_signal( ...
    delta_sig, eeg_ts, delta_burst_peak_t_anes, opt.peri_event_window);
pets_ca2_delta_sig = get_peri_event_signal( ...
    delta_ca2, ca2_ts, delta_burst_peak_t_anes, opt.peri_event_window);

% calcium signal cross-correlation with delta
[delta_ca2_xcorr_anes, delta_ca2_lag_ts, delta_ca2_xcorr_anes_cc, ...
    delta_ca2_xcorr_anes_lag] = time_mask_crosscorr( ...
    delta_sig, interp1(ca2_ts, delta_ca2, eeg_ts), opt.delta_ccg_lag, ...
    eeg_fs, mask_time_range(eeg_ts, t_anes));
[delta_ca2_xcorr_bsl, ~, delta_ca2_xcorr_bsl_cc, ...
    delta_ca2_xcorr_bsl_lag] = time_mask_crosscorr( ...
    delta_sig, interp1(ca2_ts, delta_ca2, eeg_ts), opt.delta_ccg_lag, ...
    eeg_fs, mask_time_range(eeg_ts, t_bsl));
[delta_ca2_xcorr_emer, ~, delta_ca2_xcorr_emer_cc, ...
    delta_ca2_xcorr_emer_lag] = time_mask_crosscorr( ...
    delta_sig, interp1(ca2_ts, delta_ca2, eeg_ts), opt.delta_ccg_lag, ...
    eeg_fs, mask_time_range(eeg_ts, t_emer));
[delta_ca2_xcorr_burst, ~, delta_ca2_xcorr_burst_cc, ...
    delta_ca2_xcorr_burst_lag] = time_mask_crosscorr( ...
    delta_sig, interp1(ca2_ts, delta_ca2, eeg_ts), opt.delta_ccg_lag, ...
    eeg_fs, mask_time_range(eeg_ts, mask_time_range(delta_burst_center, ...
    t_anes, [delta_burst_start, delta_burst_stop])));
[delta_ca2_xcorr_itr, ~, delta_ca2_xcorr_itr_cc, ...
    delta_ca2_xcorr_itr_lag] = time_mask_crosscorr( ...
    delta_sig, interp1(ca2_ts, delta_ca2, eeg_ts), opt.delta_ccg_lag, ...
    eeg_fs, mask_time_range(eeg_ts, mask_time_range(delta_burst_center, ...
    t_anes, [delta_burst_stop(1:end-1), delta_burst_start(2:end)])));

fprintf('%gs.', toc(t));

%% save outputs
t = tic;
fprintf('\n\tsaving ... ');
outpath = fullfile(workpath, 'outputs.mat');
if ~isfile(outpath)
    save(outpath, 'eeg_sig', 'eeg_fs', 'eeg_ts', 'ev_mrk', 'delta_sig', ...
        'delta_burst_peak_t_anes', 'delta_burst_start', 'delta_burst_stop', ...
        '-v7.3');
end    
save(outpath, 'ca2_sig', 'ca2_ts', 'ca2_fs', 'opt', 'db_out', 'ca2_fxx', 'ca2_pxx',...
    'ca2_cent_freq_burst', 'ca2_cent_freq_itr', ...
    'ca2_pxx_bsl', 'ca2_pxx_anes', 'ca2_pxx_burst', 'ca2_pxx_itr', 'ca2_pxx_emer', ...
    'ca2_pxx_bsl_per', 'ca2_pxx_anes_per', 'ca2_pxx_burst_per', 'ca2_pxx_itr_per', 'ca2_pxx_emer_per', ...
    'ca2_bp_bin', 'ca2_bp_bin_per', 'ca2_bp_bin_center', ...
    'ca2_bp_bsl', 'ca2_bp_anes', 'ca2_bp_emer', 'ca2_bp_burst', 'ca2_bp_itr', ...
    'ca2_bp_bsl_per', 'ca2_bp_anes_per', 'ca2_bp_emer_per', 'ca2_bp_burst_per', 'ca2_bp_itr_per', ...
    'delta_ca2', 'delta_burst_ca2_peak_dt', ...
    'delta_burst_ca2_peak_dt_med', 'delta_burst_ca2_peak_dt_anes', ...
    'delta_burst_ca2_peak_dt_anes_med', 'pets_ca2_delta_sig', 'pets_delta_sig', ...
    'delta_ca2_lag_ts', 'delta_ca2_xcorr_anes', 'delta_ca2_xcorr_bsl', ...
    'delta_ca2_xcorr_burst', 'delta_ca2_xcorr_itr', 'delta_ca2_xcorr_emer', ...
    'delta_ca2_xcorr_anes_cc', 'delta_ca2_xcorr_anes_lag', ...
    'delta_ca2_xcorr_emer_cc', 'delta_ca2_xcorr_emer_lag', ...
    'delta_ca2_xcorr_bsl_cc', 'delta_ca2_xcorr_bsl_lag', ...
    'delta_ca2_xcorr_burst_cc', 'delta_ca2_xcorr_burst_lag', ...
    'delta_ca2_xcorr_itr_cc', 'delta_ca2_xcorr_itr_lag', ...
    '-append');
fprintf('%gs.', toc(t));

%% visualize results
t = tic;
fprintf('\n\tplotting ... ');

plot_ca2_spectra(opt, workpath, ca2_fxx, ...
    ca2_pxx_bsl, ca2_pxx_itr, ca2_pxx_burst, ca2_pxx_anes, ca2_pxx_emer, ...
    ca2_pxx_bsl_per, ca2_pxx_itr_per, ca2_pxx_burst_per, ca2_pxx_anes_per, ca2_pxx_emer_per);

plot_delta_burst_ca2(opt, workpath, delta_burst_peak_t_anes, ...
    pets_delta_sig, pets_ca2_delta_sig);

plot_delta_ca2_xcorr(opt, workpath, delta_ca2_lag_ts, ...
    delta_ca2_xcorr_bsl, delta_ca2_xcorr_anes, delta_ca2_xcorr_burst, ...
    delta_ca2_xcorr_itr, delta_ca2_xcorr_emer);

fprintf('%gs.', toc(t));

end

function dt = get_peri_event_dt(sig, ts, ev)
[~, pt] = findpeaks(sig, ts);
[~, imin] = min(abs(ev - pt'), [], 2);
dt = pt(imin) - ev;
end

function plot_ca2_spectra(opt, workpath, fxx, ...
    spec_bsl, spec_itr, spec_burst, spec_anes, spec_emer, ...
    spec_bsl_per, spec_itr_per, spec_burst_per, spec_anes_per, spec_emer_per)

spec_itr_ch = spec_itr./spec_bsl-1;
spec_burst_ch = spec_burst./spec_bsl-1;
spec_anes_ch = spec_anes./spec_bsl-1;
spec_emer_ch = spec_emer./spec_bsl-1;

spec_itr_ch_per =  spec_itr_per./spec_bsl_per-1;
spec_burst_ch_per =spec_burst_per./spec_bsl_per-1;
spec_anes_ch_per = spec_anes_per./spec_bsl_per-1;
spec_emer_ch_per = spec_emer_per./spec_bsl_per-1;


freq_bd = unique(cat(2, opt.ca2_band_freq{:}));
figure; set(gcf, 'Position', [100, 100, 790, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot(2, 2, 1); hold on;
plot(fxx, spec_bsl, 'Color', 'k');
plot(fxx, spec_itr, 'Color', 'b');
plot(fxx, spec_burst, 'Color', 'r');
plot(fxx, spec_emer, 'Color', 'm');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Amplitude (a.u.)');
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 2); hold on;
plot(fxx, spec_itr_ch, 'Color', 'b');
plot(fxx, spec_burst_ch, 'Color', 'r');
plot(fxx, spec_emer_ch, 'Color', 'm');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Amplitude change ratio (a.u)');
legend({'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 3); hold on;
plot(fxx, spec_bsl_per, 'Color', 'k');
plot(fxx, spec_itr_per, 'Color', 'b');
plot(fxx, spec_burst_per, 'Color', 'r');
plot(fxx, spec_emer_per, 'Color', 'm');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative amplitude (%)');
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 4); hold on;
plot(fxx, spec_itr_ch_per, 'Color', 'b');
plot(fxx, spec_burst_ch_per, 'Color', 'r');
plot(fxx, spec_emer_ch_per, 'Color', 'm');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative amplitude change ratio (a.u)');
legend({'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'calcium_signal_spectra_1'), ...
    opt.fig_save_format);
close;


figure; set(gcf, 'Position', [100, 100, 790, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot(2, 2, 1); hold on;
plot(fxx, spec_bsl, 'Color', 'k');
plot(fxx, spec_anes, 'Color', 'b');
plot(fxx, spec_emer, 'Color', 'r');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Amplitude (a.u.)');
legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 2); hold on;
plot(fxx, spec_anes_ch, 'Color', 'b');
plot(fxx, spec_emer_ch, 'Color', 'r');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Amplitude change ratio (a.u)');
legend({'Anesthesia', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 3); hold on;
plot(fxx, spec_bsl_per, 'Color', 'k');
plot(fxx, spec_anes_per, 'Color', 'b');
plot(fxx, spec_emer_per, 'Color', 'r');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative amplitude (%)');
legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 4); hold on;
plot(fxx, spec_anes_ch_per, 'Color', 'b');
plot(fxx, spec_emer_ch_per, 'Color', 'r');
set(gca, 'XLim', minmax(freq_bd), 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative amplitude change ratio (a.u)');
legend({'Anesthesia', 'Emergence'}, 'Box', 'off');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'calcium_signal_spectra_2'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_ca2(opt, workpath, pkt, pets_delta_sig, pets_ca2_delta_sig)
amp_lim = [-350, 350];

clr_lim = prctile(pets_ca2_delta_sig(:), [0.1, 99.9]);
nev = length(pkt);
figure; set(gcf, 'Position', [100, 100, 400, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot('Position', [0.15, 0.54, 0.70, 0.35]); hold on;
yyaxis(gca, 'left');
shadedErrorBar(opt.peri_event_window, ...
    mean(pets_delta_sig, 'omitnan'), std(pets_delta_sig, 'omitnan')/sqrt(nev), ...
    'LineProps', {'-', 'Color', [0.3, 0.5, 1], 'LineWidth', 1});
ylabel('Delta signal (\muV)');
set(gca, 'YLim', amp_lim, 'YColor', [0.3, 0.5, 1]);
yyaxis(gca, 'right');
shadedErrorBar(opt.peri_event_window, ...
    mean(pets_ca2_delta_sig, 'omitnan'), std(pets_ca2_delta_sig, 'omitnan')/sqrt(nev), ...
    'LineProps', {'-', 'Color', [0, 0.7, 0.1], 'LineWidth', 1});
plot([0, 0], get(gca, 'YLim'), '--k');
set(gca, 'YColor', [0, 0.7, 0.1]);
xlabel('Time to delta peak (s)');
ylabel('\DeltaF/F');

subplot('Position', [0.15, 0.07, 0.70, 0.35]); hold on;
imagesc(opt.peri_event_window, (1:nev), pets_ca2_delta_sig, clr_lim); 
colormap(jet);
c = colorbar('northoutside');
c.Label.String = '\DeltaF/F (%)';
set(c, 'TickDirection', 'out', 'Box', 'off');
set(gca, 'XLim', minmax(opt.peri_event_window), 'YLim', [1, nev]);
xlabel(sprintf('Time to %s peak of delta burst (s)', opt.delta_burst_align));
ylabel('Number of delta peaks (#)');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'calcium_timecourse'), ...
    opt.fig_save_format);
close;
end

function plot_delta_ca2_xcorr(opt, workpath, lag, cc_bsl, cc_anes, ...
    cc_burst, cc_itr, cc_emer)
lag_tx = lag*1000;
figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape'); hold on;

subplot('Position', [0.60, 0.14, 0.35, 0.70]); hold on;
plot(lag_tx, cc_bsl, 'k', 'LineWidth', 1);
plot(lag_tx, cc_itr, 'b', 'LineWidth', 1);
plot(lag_tx, cc_burst, 'r', 'LineWidth', 1);
plot(lag_tx, cc_emer, 'm', 'LineWidth', 1);

set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000);
xlabel('Lag (ms)');
ylabel('Cross-correlation coefficient');
title('Calcium signal ref. to delta')
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
    'Box', 'off', 'Location', 'best');
ylim = get(gca, 'YLim');

subplot('Position', [0.12, 0.14, 0.35, 0.70]); hold on;
plot(lag_tx, cc_bsl, 'k', 'LineWidth', 1);
plot(lag_tx, cc_anes, 'b', 'LineWidth', 1);
plot(lag_tx, cc_emer, 'r', 'LineWidth', 1);

set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ylim);
xlabel('Lag (ms)');
ylabel('Cross-correlation coefficient');
title('Calcium signal ref. to delta')
legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off', 'Location', 'best');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'calcium_delta_xcorr'), ...
    opt.fig_save_format);
close;
end

