function outpath = analyze_delta_burst(eeg_sig, eeg_ts, eeg_fs, ev_mrk, opt, workpath)

t = tic;
fprintf('\n\tcalculating ... ');
if ~isfolder(workpath); mkdir(workpath); end

%% compute metrics
% period info
t_med = ev_mrk.t_med;
t_anes = ev_mrk.t_anes;
t_bsl = ev_mrk.t_bsl;
t_emer = ev_mrk.t_emer;

% notch filter
% eeg_sig = get_notched_sig(eeg_sig, eeg_fs);

% envolope
is_delta_band = strcmpi(opt.band_label, 'delta');
is_theta_band = strcmpi(opt.band_label, 'theta');
is_alpha_band = strcmpi(opt.band_label, 'alpha');
is_beta_band = strcmpi(opt.band_label, 'beta');
is_gamma_band = strcmpi(opt.band_label, 'gamma');

[band_amp_delta, delta_sig] = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(is_delta_band));
[~, gamma_sig] = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(is_gamma_band));
[~, theta_sig] = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(is_theta_band));
[~, alpha_sig] = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(is_alpha_band));
[~, beta_sig] = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(is_beta_band));
band_amp_other = get_band_envelope( ...
    eeg_sig, eeg_fs, opt.band_freq(~is_delta_band));
band_amp_delta_sm = movmean(band_amp_delta, opt.delta_smooth_window + 1/eeg_fs, ...
    'SamplePoints', eeg_ts);
delta_amp_peak_anes = findpeaks(mask_time_range(eeg_ts, t_anes, band_amp_delta_sm));

% delta burst threshold & state label
[delta_burst_thr, is_delta_burst, delta_amp_peak_gm_anes] = detect_delta_burst( ...
    delta_sig, band_amp_delta_sm, eeg_ts, t_bsl, t_anes, opt.delta_burst_method);

% peak & endpoints of delta burst
[delta_burst_start, delta_burst_stop] = clip_delta_burst( ...
    delta_sig, eeg_ts, is_delta_burst);
delta_burst_peak_t = align_delta_peak(delta_sig, eeg_ts, ...
    delta_burst_start, delta_burst_stop, opt.delta_burst_align, delta_burst_thr);
delta_burst_peak_t_anes = mask_time_range(delta_burst_peak_t, ...
    t_anes, delta_burst_peak_t);



% burst period index
delta_burst_center = (delta_burst_start + delta_burst_stop)/2;
delta_burst_duration = delta_burst_stop - delta_burst_start + 1/eeg_fs;
delta_burst_amp = get_delta_burst_amp(abs(delta_sig), eeg_ts, ...
    delta_burst_start, delta_burst_stop);
delta_burst_itr = (delta_burst_start(2:end) + delta_burst_stop(1:end-1))/2;
delta_peak_num = get_delta_cycle(delta_sig, eeg_ts, ...
    mask_time_range(delta_burst_center, t_anes, delta_burst_start), ...
    mask_time_range(delta_burst_center, t_anes, delta_burst_stop));

% burst timecourse
[delta_burst_bin_count, delta_burst_bin_duration, delta_burst_bin_center] = ...
    get_event_timecourse(delta_burst_duration, delta_burst_center, opt.delta_burst_bin);
[~, delta_burst_bin_amp, ~] = ...
    get_event_timecourse(delta_burst_amp, delta_burst_center, opt.delta_burst_bin);
delta_burst_duration_gm_anes = fitdist(mask_time_range( ...
    delta_burst_center, t_anes, delta_burst_duration), 'Lognormal');
delta_burst_interval_gm_anes = fitdist(diff(mask_time_range( ...
    delta_burst_center, t_anes, delta_burst_center)), 'Lognormal');

% spectrum
freq_limit = minmax(cat(2, opt.band_freq{:}));
[spectrum_sxx, spectrum_fxx, spectrum_txx] = multi_taper_spectrum( ...
    eeg_sig, opt.fft_window, opt.fft_step, eeg_fs, freq_limit, eeg_ts);

spectrum_pxx = step_window_spectrum(eeg_sig, opt.fft_window, opt.fft_step, ...
    eeg_fs, round(opt.apply_win_len * eeg_fs), freq_limit, eeg_ts);
spectrum_burst = median(interp1(spectrum_txx, spectrum_pxx', ...
    mask_time_range(delta_burst_center, t_anes, delta_burst_center), ...
    'nearest'), 1, 'omitnan');
spectrum_itr = median(interp1(spectrum_txx, spectrum_pxx', ...
    mask_time_range(delta_burst_itr, t_anes, delta_burst_itr), ...
    'nearest'), 1, 'omitnan');
spectrum_bsl = median( ...
    spectrum_pxx(:, mask_time_range(spectrum_txx, t_bsl)), 2, 'omitnan')';
spectrum_anes = median( ...
    spectrum_pxx(:, mask_time_range(spectrum_txx, t_anes)), 2, 'omitnan')';
spectrum_emer = median( ...
    spectrum_pxx(:, mask_time_range(spectrum_txx, t_emer)), 2, 'omitnan')';

spectrum_pxx_per = spectrum_pxx./trapz(spectrum_fxx, spectrum_pxx, 1);
spectrum_burst_per = median(interp1(spectrum_txx, spectrum_pxx_per', ...
    mask_time_range(delta_burst_center, t_anes, delta_burst_center), ...
    'nearest'), 1, 'omitnan');
spectrum_itr_per = median(interp1(spectrum_txx, spectrum_pxx_per', ...
    mask_time_range(delta_burst_itr, t_anes, delta_burst_itr), ...
    'nearest'), 1, 'omitnan');
spectrum_bsl_per = median( ...
    spectrum_pxx_per(:, mask_time_range(spectrum_txx, t_bsl)), 2, 'omitnan')';
spectrum_anes_per = median( ...
    spectrum_pxx_per(:, mask_time_range(spectrum_txx, t_anes)), 2, 'omitnan')';
spectrum_emer_per = median( ...
    spectrum_pxx_per(:, mask_time_range(spectrum_txx, t_emer)), 2, 'omitnan')';
[~, max_index] = max(spectrum_burst_per);
cent_freq_burst = spectrum_fxx(max_index);

% spectrum_bsl_per = spectrum_bsl/trapz(spectrum_fxx, spectrum_bsl);
% spectrum_itr_per = spectrum_itr/trapz(spectrum_fxx,spectrum_itr);
% spectrum_burst_per = spectrum_burst/trapz(spectrum_fxx, spectrum_burst);
% spectrum_anes_per = spectrum_anes/trapz(spectrum_fxx, spectrum_anes);
% spectrum_emer_per= spectrum_emer/trapz(spectrum_fxx, spectrum_emer);


% bandpower
bandpower_raw = 10*log10(integrate_psd_bandpower(spectrum_pxx, spectrum_fxx, opt.band_freq));
bandpower_raw_per = integrate_psd_bandpower(spectrum_pxx_per, spectrum_fxx, opt.band_freq);

[~, bandpower_bin, bandpower_bin_center] = get_event_timecourse( ...
    bandpower_raw, spectrum_txx', opt.band_power_bin);
[~, bandpower_bin_per, ~] = get_event_timecourse( ...
    bandpower_raw_per, spectrum_txx', opt.band_power_bin);


bandpower_bsl = 10*log10(integrate_psd_bandpower(spectrum_bsl, spectrum_fxx, opt.band_freq));
bandpower_anes = 10*log10(integrate_psd_bandpower(spectrum_anes, spectrum_fxx, opt.band_freq));
bandpower_burst = 10*log10(integrate_psd_bandpower(spectrum_burst, spectrum_fxx, opt.band_freq));
bandpower_itr = 10*log10(integrate_psd_bandpower(spectrum_itr, spectrum_fxx, opt.band_freq));
bandpower_emer = 10*log10(integrate_psd_bandpower(spectrum_emer, spectrum_fxx, opt.band_freq));


bandpower_bsl_per = integrate_psd_bandpower(spectrum_bsl_per, spectrum_fxx, opt.band_freq);
bandpower_anes_per = integrate_psd_bandpower(spectrum_anes_per, spectrum_fxx, opt.band_freq);
bandpower_burst_per = integrate_psd_bandpower(spectrum_burst_per, spectrum_fxx, opt.band_freq);
bandpower_itr_per = integrate_psd_bandpower(spectrum_itr_per, spectrum_fxx, opt.band_freq);
bandpower_emer_per = integrate_psd_bandpower(spectrum_emer_per, spectrum_fxx, opt.band_freq);

% band amplitude timecourse
pets_delta_amp = get_peri_event_signal( ...
    band_amp_delta, eeg_ts, delta_burst_peak_t_anes, opt.peri_event_window);
pets_other_amp = get_peri_event_signal( ...
    band_amp_other, eeg_ts, delta_burst_peak_t_anes, opt.peri_event_window);

% band amplitude cross-correlation
[band_amp_xcorr_anes, band_amp_lag_ts, band_amp_xcorr_anes_cc, ...
    band_amp_xcorr_anes_lag] = time_mask_crosscorr( ...
    band_amp_delta, band_amp_other, opt.delta_ccg_lag, eeg_fs, ...
    mask_time_range(eeg_ts, t_anes));
[band_amp_xcorr_bsl, ~, band_amp_xcorr_bsl_cc, ...
    band_amp_xcorr_bsl_lag] = time_mask_crosscorr( ...
    band_amp_delta, band_amp_other, opt.delta_ccg_lag, eeg_fs, ...
    mask_time_range(eeg_ts, t_bsl));
[band_amp_xcorr_emer, ~, band_amp_xcorr_emer_cc, ...
    band_amp_xcorr_emer_lag] = time_mask_crosscorr( ...
    band_amp_delta, band_amp_other, opt.delta_ccg_lag, eeg_fs, ...
    mask_time_range(eeg_ts, t_emer));
[band_amp_xcorr_burst, ~, band_amp_xcorr_burst_cc, ...
    band_amp_xcorr_burst_lag] = time_mask_crosscorr( ...
    band_amp_delta, band_amp_other, opt.delta_ccg_lag, eeg_fs, ...
    mask_time_range(eeg_ts, mask_time_range(delta_burst_center, ...
    t_anes, [delta_burst_start, delta_burst_stop])));
[band_amp_xcorr_itr, ~, band_amp_xcorr_itr_cc, ...
    band_amp_xcorr_itr_lag] = time_mask_crosscorr( ...
    band_amp_delta, band_amp_other, opt.delta_ccg_lag, eeg_fs, ...
    mask_time_range(eeg_ts, mask_time_range(delta_burst_itr, ...
    t_anes, [delta_burst_stop(1:end-1), delta_burst_start(2:end)])));

fprintf('%gs.', toc(t));

%% save outputs
t = tic;
fprintf('\n\tsaving ... ');
outpath = fullfile(workpath, 'outputs.mat');
save(outpath, 'eeg_sig', 'eeg_fs', 'eeg_ts', 'ev_mrk', 'opt', 'delta_sig', 'gamma_sig', ...
    'theta_sig', 'alpha_sig', 'beta_sig', ...
    'band_amp_delta_sm', 'delta_amp_peak_anes', 'delta_amp_peak_gm_anes', ...
    'delta_burst_thr', 'delta_burst_start', 'delta_burst_stop', ...
    'delta_burst_amp', 'delta_peak_num', 'cent_freq_burst', ...
    'delta_burst_bin_center', 'delta_burst_bin_count', ...
    'delta_burst_bin_duration', 'delta_burst_bin_amp', ...
    'delta_burst_duration_gm_anes', 'delta_burst_interval_gm_anes', ...
    'delta_burst_peak_t', 'delta_burst_peak_t_anes', ...
    'spectrum_sxx', 'spectrum_fxx', 'spectrum_txx','spectrum_pxx', ...
    'spectrum_bsl', 'spectrum_anes', 'spectrum_itr', 'spectrum_burst', 'spectrum_emer', ...
    'spectrum_bsl_per', 'spectrum_anes_per', 'spectrum_itr_per', 'spectrum_burst_per', 'spectrum_emer_per', ...
    'bandpower_bsl', 'bandpower_anes', 'bandpower_itr', 'bandpower_burst', 'bandpower_emer', ...
    'bandpower_bsl_per', 'bandpower_anes_per', 'bandpower_itr_per', 'bandpower_burst_per', 'bandpower_emer_per', ...
    'bandpower_bin', 'bandpower_bin_per', 'bandpower_bin_center', 'pets_delta_amp', 'pets_other_amp', ...
    'band_amp_lag_ts', 'band_amp_xcorr_bsl', 'band_amp_xcorr_anes', ...
    'band_amp_xcorr_burst', 'band_amp_xcorr_itr', 'band_amp_xcorr_emer', ...
    'band_amp_xcorr_anes_cc', 'band_amp_xcorr_anes_lag', ...
    'band_amp_xcorr_emer_cc', 'band_amp_xcorr_emer_lag', ...
    'band_amp_xcorr_bsl_cc', 'band_amp_xcorr_bsl_lag', ...
    'band_amp_xcorr_burst_cc', 'band_amp_xcorr_burst_lag', ...
    'band_amp_xcorr_itr_cc', 'band_amp_xcorr_itr_lag', ...
    '-v7.3');
fprintf('%gs.', toc(t));

%% visualize results
t = tic;
fprintf('\n\tplotting ... ');

plot_delta_burst_timecourse(opt, workpath, delta_sig, band_amp_delta_sm, ...
    eeg_ts, t_med, t_anes, delta_burst_thr, delta_burst_peak_t, ...
    delta_burst_start, delta_burst_stop, delta_burst_bin_center, ...
    delta_burst_bin_count, delta_burst_bin_duration, delta_burst_bin_amp, ...
    spectrum_txx, spectrum_fxx, spectrum_sxx);

plot_delta_burst_threshold(opt, workpath, band_amp_delta_sm, eeg_ts, ...
    t_bsl, t_anes, delta_burst_thr, delta_amp_peak_anes, delta_amp_peak_gm_anes);

plot_delta_burst_trace(opt, workpath, delta_sig, eeg_ts, ...
    delta_burst_peak_t_anes);

plot_delta_burst_metrics(opt, workpath, t_anes, ...
    delta_burst_start, delta_burst_stop, ...
    delta_burst_duration_gm_anes, delta_burst_interval_gm_anes);

plot_delta_burst_spectra(opt, workpath, spectrum_fxx, ...
    spectrum_bsl, spectrum_itr, spectrum_burst, spectrum_anes, spectrum_emer, ...
    spectrum_bsl_per, spectrum_itr_per, spectrum_burst_per, spectrum_anes_per, spectrum_emer_per);

plot_delta_burst_band_amp(opt, workpath, delta_burst_peak_t_anes, ...
    pets_delta_amp, pets_other_amp);

plot_delta_burst_xcorr(opt, workpath, band_amp_lag_ts, ...
    band_amp_xcorr_bsl, band_amp_xcorr_anes, band_amp_xcorr_burst, ...
    band_amp_xcorr_itr, band_amp_xcorr_emer);

fprintf('%gs.', toc(t));

end

function [thr, states, gm] = detect_delta_burst(sig, amp, ts, tbsl, tanes, method)
rng(1024);
gm = fitgmdist(findpeaks(mask_time_range(ts, tanes, amp)), 2, ...
    'RegularizationValue', 0.1);
if strcmp(method, 'anes_2gm')
    kfold = 2;
    thr =min(gm.mu) * kfold;
elseif strcmp(method, 'bsl_rms')
    kfold = 2;
    thr = sqrt(mean(mask_time_range(ts, tbsl, sig).^2)) * kfold;
else
    error('Not implemented.');
end
states = amp > thr;
% starts = find(diff([0; states])==1);
% stops = find(diff([states; 0])==-1);
% for tt = 1:length(starts)
%     if ~any(findpeaks(sig(starts(tt):stops(tt))))
%         states(starts(tt):stops(tt)) = false;
%     end
% end
end

function pk_amp = get_delta_burst_amp(sig, ts, tstart, tstop)
[~, starts] = ismember(tstart, ts);
[~, stops] = ismember(tstop, ts);
nb = length(starts);
pk_amp = zeros(nb, 1);
for kk = 1:nb
    pk_amp(kk) = max(sig(starts(kk):stops(kk)));    
end
end

function pk_num = get_delta_cycle(sig, ts, tstart, tstop)
[~, starts] = ismember(tstart, ts);
[~, stops] = ismember(tstop, ts);
nb = length(starts);
pk_num = zeros(nb, 1);
for kk = 1:nb
    [~, locs] = findpeaks(sig(starts(kk):stops(kk))); 
    pk_num(kk) = length(locs);
end
end

function [pk_t, pk_val] = align_delta_peak(sig, ts, tstart, tstop, method, thr)
ttrans = reshape([tstart, tstop]', [], 1);
[pk_val_o, pk_t_o] = findpeaks(sig, ts);
[~, ~, idx_bin] = histcounts(pk_t_o, ttrans);
is_pk = false(size(pk_t_o));
switch method
    case 'first'
        [~, idx_unique] = unique(idx_bin);
        idx_unique(mod(idx_bin(idx_unique), 2)==0) = [];
        is_pk(idx_unique) = true;
    case 'highest'
        unique_bin = unique(idx_bin);
        for kk = 1:length(unique_bin)
            if mod(unique_bin(kk), 2)==0
                continue;
            end
            tmp_idx = find(idx_bin==unique_bin(kk));
            [~, imax] = max(pk_val_o(idx_bin==unique_bin(kk)));
            is_pk(tmp_idx(imax)) = true;
        end
    case 'first-highest'
        is_surp = pk_val_o > thr;
        unique_bin = unique(idx_bin);
        for kk = 1:length(unique_bin)
            if mod(unique_bin(kk), 2)==0
                continue;
            end
            tmp_idx = find(idx_bin==unique_bin(kk));
            imax = find(is_surp(tmp_idx), 1);
            if isempty(imax)
                [~, imax] = max(pk_val_o(idx_bin==unique_bin(kk)));
            end
            is_pk(tmp_idx(imax)) = true;
        end
    case 'first-skip'
        is_surp = pk_val_o > thr;
        unique_bin = unique(idx_bin);
        for kk = 1:length(unique_bin)
            if mod(unique_bin(kk), 2)==0
                continue;
            end
            tmp_idx = find(idx_bin==unique_bin(kk));
            imax = find(is_surp(tmp_idx), 1);
            if isempty(imax)
                imax = 1;
                pk_t_o(tmp_idx(imax)) = nan;
                pk_val_o(tmp_idx(imax)) = nan;
            end
            is_pk(tmp_idx(imax)) = true;
        end
end
pk_t = pk_t_o(is_pk);
pk_val = pk_val_o(is_pk);
end

function [tstart, tstop] = clip_delta_burst(sig, ts, isstate)
tstart = ts(diff([0; isstate])==1);
tstop = ts(diff([isstate; 0])==-1);

[~, tmin] = findpeaks(-sig, ts);
ispass = tstart < tmin(1) | tstop > tmin(end);
nt = length(tstart);
for tt = 1:nt
    if ispass(tt); continue; end
    tmp = tstart(tt) - tmin;
    tstart(tt) = tstart(tt) - min(tmp(tmp >= 0));
    tmp = tmin - tstop(tt);
    tstop(tt) = tstop(tt) + min(tmp(tmp >= 0));
end
ispass = find(tstart(2:end) - tstop(1:end-1) <= 0);
tstop(ispass) = [];
tstart(ispass+1) = [];

[~, tmax] = findpeaks(sig, ts);
nt = length(tstart);
ispass = false(nt, 1);
for tt = 1:nt
    if ~any(tmax - tstart(tt) >= 0 & tstop(tt) - tmax >= 0)
        ispass(tt) = true;
    end
end
tstart = tstart(~ispass);
tstop = tstop(~ispass);
end

function plot_delta_burst_timecourse(opt, workpath, sig, amp, ts, tmed, tanes, ...
    thr, pkt, tstart, tstop, bcnt, bcounts, bdur, bamp, txx, fxx, pxx)
tw = 300;
tx = (ceil((min(ts) - tmed)/tw):floor((max(ts) - tmed)/tw));
ttick = (ceil((min(ts) - tmed)/tw):floor((max(ts) - tmed)/tw))*tw+tmed;
tlabel = tx*tw/60;
figure; set(gcf, 'Position', [0, 100, 1120, 790], 'Visible', false, ...
    'PaperOrientation', 'landscape');

subplot('Position', [0.10, 0.78, 0.85, 0.15]); hold on;
temp = griddedInterpolant(ts, sig);
plot(ts, sig);
plot(ts, amp);
plot(minmax(ts'), [1, 1]*thr, '-k', 'LineWidth', 1);
plot(pkt, temp(pkt), '.r');
plot(tstart, temp(tstart), '.m');
plot(tstop, temp(tstop), '.m');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.005], ...
    'XLim', tanes, 'XTick', ttick, 'XTickLabel', tlabel);
xlabel('Time in anesthesia (min)');
ylabel('Delta amplitude (\muV)');

subplot('Position', [0.10, 0.61, 0.85, 0.12]); hold on;
plot(bcnt, bdur, '-m');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.005], 'YLim', [-0.1, 4.1], ...
    'XLim', minmax(txx(:)'), 'XTick', ttick, 'XTickLabel', tlabel);
ylabel('Duration (s)');

subplot('Position', [0.10, 0.46, 0.85, 0.12]); hold on;
plot(bcnt, bcounts/opt.delta_burst_bin, '-k');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.005], 'YLim', [-0.1, 1.2], ...
    'XLim', minmax(txx(:)'), 'XTick', ttick, 'XTickLabel', tlabel);
ylabel('Bouts (s^{-1})');

subplot('Position', [0.10, 0.31, 0.85, 0.12]); hold on;
plot(bcnt, bamp, '-c');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.005], 'YLim', [110, 400], ...
    'XLim', minmax(txx(:)'), 'XTick', ttick, 'XTickLabel', tlabel);
ylabel('Amplitude (\muV)');

subplot('Position', [0.10, 0.07, 0.85, 0.20]); hold on;
imagesc(txx, fxx, 10*log10(pxx), [0, 40]); colormap(jet);
xlabel('Time from injection (min)');
ylabel('Frequency (Hz)');
set(gca, 'TickDir', 'out', 'TickLength', [0.005, 0.005], 'YLim', minmax(fxx(:)'), ...
    'XLim', minmax(txx(:)'), 'XTick', ttick, 'XTickLabel', tlabel);

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_timecourse'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_threshold(opt, workpath, amp, ts, tbsl, tanes, thr, pk, pk_gm)
amp_bin = 10;
amp_rg = [0, 400];

amp_edges = (amp_rg(1):amp_bin:amp_rg(2))';
amp_edges_high = (amp_rg(1):amp_bin/10:amp_rg(2))';
figure; set(gcf, 'Position', [100, 100, 400, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot('Position', [0.20, 0.54, 0.70, 0.35]); hold on;
histogram(mask_time_range(ts, tbsl, amp), amp_edges, ...
    'Normalization', 'probability', 'FaceColor', [0.5, 0.8, 1]);
histogram(mask_time_range(ts, tanes, amp), amp_edges, ...
    'Normalization', 'probability', 'FaceColor', [1, 0.3, 0.3]); 
plot([1, 1]*thr, get(gca, 'YLim'), '--k', 'LineWidth', 1.5);
xlabel('Amplitude of delta envelope (\muV)');
ylabel('Probability');
title(sprintf('Threshold: %.2f \\muV', thr));
legend({'Pre-injection', 'Anesthesia', opt.delta_burst_method}, ...
    'Box', 'off', 'Interpreter', 'none');

subplot('Position', [0.20, 0.07, 0.70, 0.35]); hold on;
histogram(pk, amp_edges, 'Normalization', 'probability', ...
    'FaceColor', [0.7, 0.7, 0.7]); 
plot(amp_edges_high, pdf(pk_gm, amp_edges_high)*amp_bin, ...
    'LineWidth', 1.5);
plot(amp_edges_high, pdf(gmdistribution(pk_gm.mu(1), pk_gm.Sigma(1)), ...
    amp_edges_high)*pk_gm.ComponentProportion(1)*amp_bin, ...
    'LineWidth', 1.5);
plot(amp_edges_high, pdf(gmdistribution(pk_gm.mu(2), pk_gm.Sigma(2)), ...
    amp_edges_high)*pk_gm.ComponentProportion(2)*amp_bin, ...
    'LineWidth', 1.5);
plot([1, 1]*min(pk_gm.mu)*2, get(gca, 'YLim'), '--k', 'LineWidth', 1.5);
xlabel('Peak amplitude of anesthesia delta envelope (\muV)');
ylabel('Probability');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_threshold'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_trace(opt, workpath, sig, ts, pkt)
pk_win = (-1:0.001:2);

nev = length(pkt);
traces = interp1(ts, sig, pkt+pk_win);
figure; set(gcf, 'Position', [100, 100, 400, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot('Position', [0.20, 0.54, 0.70, 0.35]); hold on;
plot(pk_win, traces', 'Color', [.7, .7, .7], 'LineWidth', 0.3);
plot(pk_win, mean(traces, 1), '-r', 'LineWidth', 1.3);
xlabel(sprintf('Time to %s peak of delta burst (s)', opt.delta_burst_align));
ylabel('Amplitude (\muV)');
title(sprintf('N = %d', length(pkt)));

subplot('Position', [0.20, 0.07, 0.70, 0.35]); hold on;
imagesc(pk_win, (1:nev), traces, [-350, 350]); colormap(jet);
c = colorbar('northoutside');
c.Label.String = 'Amplitude (\muV)';
set(c, 'TickDirection', 'out', 'Box', 'off');
set(gca, 'XLim', minmax(pk_win), 'YLim', [1, nev]);
xlabel(sprintf('Time to %s peak of delta burst (s)', opt.delta_burst_align));
ylabel('Number of delta peaks (#)');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_trace'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_metrics(opt, workpath, tanes, tstart, tstop, dur_gm, itv_gm)
dur_bin = 0.1;
dur_rg = [0, 4];
itv_bin = 0.4;
itv_rg = [0, 20];

cnt = (tstart+tstop)/2;
dur_edges = (dur_rg(1):dur_bin:dur_rg(2));
itv_edges = (itv_rg(1):itv_bin:itv_rg(2));
dur_edges_high = (dur_rg(1):dur_bin/10:dur_rg(2));
itv_edges_high = (itv_rg(1):itv_bin/10:itv_rg(2));

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');

subplot('Position', [0.12, 0.14, 0.35, 0.70]); hold on;
temp = mask_time_range(cnt, tanes, tstop-tstart);
histogram(temp, dur_edges, ...
    'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], ...
    'FaceAlpha', 0.4);
plot(dur_edges_high, pdf(dur_gm, dur_edges_high)*dur_bin, ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
xlabel('Delta burst duration (s)');
ylabel('Probability');
title(sprintf('log-normal \\mu=%.2f, median=%.2f', ...
    exp(dur_gm.mu), median(temp)));

subplot('Position', [0.60, 0.14, 0.35, 0.70]); hold on;
temp = diff(mask_time_range(cnt, tanes, cnt));
histogram(temp, itv_edges, ...
    'Normalization', 'probability', 'FaceColor', [0.8500, 0.3250, 0.0980], ...
    'FaceAlpha', 0.4);
plot(itv_edges_high, ...
    pdf(itv_gm, itv_edges_high)*itv_bin, ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
xlabel('Delta burst interval (s)');
ylabel('Probability');
title(sprintf('log-normal \\mu=%.2f, median=%.2f', ...
    exp(itv_gm.mu), median(temp)));

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_metrics'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_spectra(opt, workpath, fxx, ...
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

freq_bd = unique(cat(2, opt.band_freq{:}));
figure; set(gcf, 'Position', [100, 100, 790, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot(2, 2, 1);
semilogx(fxx, 10*log10(spec_bsl), 'Color', 'k');
hold on;
semilogx(fxx, 10*log10(spec_itr), 'Color', 'b');
hold on;
semilogx(fxx, 10*log10(spec_burst), 'Color', 'r');
hold on;
semilogx(fxx, 10*log10(spec_emer), 'Color', 'm');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 2);
semilogx(fxx, spec_itr_ch, 'Color', 'b');
hold on;
semilogx(fxx, spec_burst_ch, 'Color', 'r');
hold on;
semilogx(fxx, spec_emer_ch, 'Color', 'm');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Power change ratio (a.u.)');
legend({'Inter-burst', 'Delta burst','Emergence'}, 'Box', 'off');

subplot(2, 2, 3);
semilogx(fxx, spec_bsl_per, 'Color', 'k');
hold on;
semilogx(fxx, spec_itr_per, 'Color', 'b');
hold on;
semilogx(fxx, spec_burst_per, 'Color', 'r');
hold on;
semilogx(fxx, spec_emer_per, 'Color', 'm');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative PSD (%)');
legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 4);
semilogx(fxx, spec_itr_ch_per, 'Color', 'b');
hold on;
semilogx(fxx, spec_burst_ch_per, 'Color', 'r');
hold on;
semilogx(fxx, spec_emer_ch_per, 'Color', 'm');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative power change ratio (a.u.)');
legend({'Inter-burst', 'Delta burst','Emergence'}, 'Box', 'off');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_spectra_1'), ...
    opt.fig_save_format);
close;


figure; set(gcf, 'Position', [100, 100, 790, 800], 'Visible', false, ...
    'PaperOrientation', 'portrait');

subplot(2, 2, 1);
semilogx(fxx, 10*log10(spec_bsl), 'Color', 'k');
hold on;
semilogx(fxx, 10*log10(spec_anes), 'Color', 'b');
hold on;
semilogx(fxx, 10*log10(spec_emer), 'Color', 'r');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
legend({'Pre-injection', 'Anesthesia','Emergence'}, 'Box', 'off');

subplot(2, 2, 2);
semilogx(fxx, spec_anes_ch, 'Color', 'b');
hold on;
semilogx(fxx, spec_emer_ch, 'Color', 'r');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Power change ratio');
legend({'Anesthesia', 'Emergence'}, 'Box', 'off');

subplot(2, 2, 3);
semilogx(fxx, spec_bsl_per, 'Color', 'k');
hold on;
semilogx(fxx, spec_anes_per, 'Color', 'b');
hold on;
semilogx(fxx, spec_emer_per, 'Color', 'r');
set(gca, 'XLim', [0, 100],'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative PSD');
legend({'Pre-injection', 'Anesthesia','Emergence'}, 'Box', 'off');

subplot(2, 2, 4);
semilogx(fxx, spec_anes_ch_per, 'Color', 'b');
hold on;
semilogx(fxx, spec_emer_ch_per, 'Color', 'r');
set(gca, 'XLim', [0, 100], 'XTick', freq_bd, ...
    'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
xlabel('Frequency (Hz)');
ylabel('Relative power change ratio (a.u.)');
legend({'Anesthesia', 'Emergence'}, 'Box', 'off');

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_spectra_2'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_xcorr(opt, workpath, lag, cc_bsl, cc_anes, ...
    cc_burst, cc_itr, cc_emer)
lag_tx = lag*1000;
band_label_other = opt.band_label(~strcmpi(opt.band_label, 'delta'));
num_band = length(band_label_other);

figure; set(gcf, 'Position', [100, 100, 900, 300], 'Visible', false, ...
    'PaperOrientation', 'landscape');

subplot(1, 3, 1); hold on;
plot(lag_tx, cc_anes, 'LineWidth', 1);
set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000);
xlabel('Lag (ms)');
ylabel('Cross-correlation coefficient');
title('Anesthesia');
ylim = get(gca, 'YLim');

subplot(1, 3, 2); hold on;
plot(lag_tx, cc_bsl, 'LineWidth', 1);
set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ylim);
xlabel('Lag (ms)');
ylabel('Cross-correlation coefficient');
title('Pre-injection');

subplot(1, 3, 3); hold on;
plot(lag_tx, cc_emer, 'LineWidth', 1);
set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ylim);
xlabel('Lag (ms)');
ylabel('Cross-correlation coefficient');
title('Emergence');

legend(opt.band_label(~strcmpi(opt.band_label, 'delta')), ...
    'Box', 'off', 'Location', 'north', 'NumColumns', 2);
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_xcorr'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 1100, 660], 'Visible', false, ...
    'PaperOrientation', 'landscape');

for bb = 1:num_band
    subplot(2, 3, bb); hold on;
    plot(lag_tx, cc_bsl(:, bb), 'k', 'LineWidth', 1);
    plot(lag_tx, cc_itr(:, bb), 'b', 'LineWidth', 1);
    plot(lag_tx, cc_burst(:, bb), 'r', 'LineWidth', 1);
    plot(lag_tx, cc_emer(:, bb), 'm', 'LineWidth', 1);
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, ...
        'YLim', ylim * 1.1);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title(band_label_other{bb}, 'Interpreter', 'none');
    if bb == 1
        legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
            'Box', 'off', 'Location', 'best');
    end
end

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_xcorr_by_band'), ...
    opt.fig_save_format);
close;
end

function plot_delta_burst_band_amp(opt, workpath, pkt, pets_delta_amp, pets_other_amp)
clrs = get_default_colors;
band_label_other = opt.band_label(~strcmpi(opt.band_label, 'delta'));

num_band = length(band_label_other);
nev = length(pkt);
figure; set(gcf, 'Position', [100, 100, 1100, 660], 'Visible', false, ...
    'PaperOrientation', 'landscape');

for bb = 1:num_band
    subplot(2, 3, bb);
    yyaxis(gca, 'left');
    shadedErrorBar(opt.peri_event_window, ...
        mean(pets_delta_amp, 'omitnan'), std(pets_delta_amp, 'omitnan')/sqrt(nev), ...
        'LineProps', {'-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1});
    ylabel('Amplitude of delta (\muV)');
    set(gca, 'YColor', [0.5, 0.5, 0.5]);
    yyaxis(gca, 'right');
    shadedErrorBar(opt.peri_event_window, mean(pets_other_amp(:, :, bb), 'omitnan'), ...
        std(pets_other_amp(:, :, bb), 'omitnan')/sqrt(nev), ...
        'LineProps', {'-', 'Color', clrs(bb, :), 'LineWidth', 1});
    xline(0, '--k');
    set(gca, 'YColor', clrs(bb, :));
    xlabel('Time to delta peak (s)');
    ylabel(sprintf('Amplitude of %s (\\muV)', band_label_other{bb}));
end

sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'delta_burst_band_amp'), ...
    opt.fig_save_format);
close;
end
