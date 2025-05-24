function [outpath, stats] = summarize_multi_outputs(outpaths, workpath)

t = tic;
fprintf('\nsummarizing ... ');
if ~isfolder(workpath); mkdir(workpath); end
do_ca2 = false;
do_spk = false;

%% gather recording outputs
num_out = length(outpaths);
stats = struct;
for nn = 1:num_out
    stats(nn).OutPath = outpaths{nn};
    var_list = fieldnames(matfile(fullfile(outpaths{nn}, 'outputs.mat')));
    if ~ismember('delta_burst_start', var_list)
        continue;
    end

    load(fullfile(outpaths{nn}, 'outputs.mat'), 'opt', 'ev_mrk', ...
        'delta_burst_start', 'delta_burst_stop', 'delta_peak_num', 'cent_freq_burst',...
        'delta_burst_duration_gm_anes',  'delta_burst_interval_gm_anes',  ...
        'spectrum_fxx', 'spectrum_bsl', 'spectrum_anes', 'spectrum_itr', 'spectrum_burst', 'spectrum_emer', ...
        'spectrum_bsl_per', 'spectrum_anes_per', 'spectrum_itr_per', 'spectrum_burst_per', 'spectrum_emer_per', ...
        'bandpower_bin', 'bandpower_bin_center', 'bandpower_bin_per', ...
        'bandpower_bsl', 'bandpower_anes', 'bandpower_emer', 'bandpower_burst', 'bandpower_itr', ...
        'bandpower_bsl_per', 'bandpower_anes_per', 'bandpower_emer_per', 'bandpower_burst_per', 'bandpower_itr_per', ...
        'pets_delta_amp', 'pets_other_amp', ...
        'band_amp_lag_ts', 'band_amp_xcorr_bsl', 'band_amp_xcorr_anes', ...
        'band_amp_xcorr_burst', 'band_amp_xcorr_itr', 'band_amp_xcorr_emer', ...
        'band_amp_xcorr_anes_cc', 'band_amp_xcorr_anes_lag', ...
        'band_amp_xcorr_bsl_cc', 'band_amp_xcorr_bsl_lag', ...
        'band_amp_xcorr_burst_cc', 'band_amp_xcorr_burst_lag', ...
        'band_amp_xcorr_itr_cc', 'band_amp_xcorr_itr_lag', ...
        'band_amp_xcorr_emer_cc', 'band_amp_xcorr_emer_lag');

    delta_burst_center = (delta_burst_start+delta_burst_stop)/2;
    delta_burst_duration = delta_burst_stop-delta_burst_start;
    stats(nn).Duration = [exp(delta_burst_duration_gm_anes.mu), ...
        median(mask_time_range(delta_burst_center, ...
        ev_mrk.t_anes, delta_burst_duration))];
    stats(nn).Interval = [exp(delta_burst_interval_gm_anes.mu), ...
        median(diff(mask_time_range(delta_burst_center, ...
        ev_mrk.t_anes, delta_burst_center)))];
    stats(nn).DeltaPeakNum = delta_peak_num;
  

    stats(nn).SpectrumBSL = spectrum_bsl;
    stats(nn).SpectrumInter = spectrum_itr;
    stats(nn).SpectrumBurst = spectrum_burst;
    stats(nn).SpectrumAnes = spectrum_anes;
    stats(nn).SpectrumEmer = spectrum_emer;
    stats(nn).SpectrumInterCh = spectrum_itr./spectrum_bsl-1;
    stats(nn).SpectrumBurstCh = spectrum_burst./spectrum_bsl-1;
    stats(nn).SpectrumAnesCh = spectrum_anes./spectrum_bsl-1;
    stats(nn).SpectrumEmerCh = spectrum_emer./spectrum_bsl-1;
    

    stats(nn).SpectrumBSLPer = spectrum_bsl_per;
    stats(nn).SpectrumInterPer = spectrum_itr_per;
    stats(nn).SpectrumBurstPer = spectrum_burst_per;
    stats(nn).SpectrumAnesPer = spectrum_anes_per;
    stats(nn).SpectrumEmerPer = spectrum_emer_per;

    stats(nn).CentFreqBurst = cent_freq_burst;
    
    stats(nn).SpectrumInterChPer = spectrum_itr_per./spectrum_bsl_per-1;
    stats(nn).SpectrumBurstChPer = spectrum_burst_per./spectrum_bsl_per-1;
    stats(nn).SpectrumAnesChPer = spectrum_anes_per./spectrum_bsl_per-1;
    stats(nn).SpectrumEmerChPer = spectrum_emer_per./spectrum_bsl_per-1;


    stats(nn).BandpowerBin = bandpower_bin;
    stats(nn).BandpowerBinPer = bandpower_bin_per;
    stats(nn).BandpowerBinCnt = bandpower_bin_center;

    stats(nn).BandpowerBSL = bandpower_bsl;
    stats(nn).BandpowerAnes = bandpower_anes;
    stats(nn).BandpowerInter = bandpower_itr;
    stats(nn).BandpowerBurst = bandpower_burst;
    stats(nn).BandpowerEmer = bandpower_emer;

    stats(nn).BandpowerBSLPer = bandpower_bsl_per;
    stats(nn).BandpowerAnesPer = bandpower_anes_per;
    stats(nn).BandpowerInterPer = bandpower_itr_per;
    stats(nn).BandpowerBurstPer = bandpower_burst_per;
    stats(nn).BandpowerEmerPer = bandpower_emer_per;

    stats(nn).PSTHDeltaAmp = mean(pets_delta_amp, 1);
    stats(nn).PSTHOtherAmp = mean(pets_other_amp, 1);

    stats(nn).CrossCorrBSL = permute(band_amp_xcorr_bsl, [3, 1, 2]);
    stats(nn).CrossCorrAnes = permute(band_amp_xcorr_anes, [3, 1, 2]);
    stats(nn).CrossCorrInter = permute(band_amp_xcorr_itr, [3, 1, 2]);
    stats(nn).CrossCorrBurst = permute(band_amp_xcorr_burst, [3, 1, 2]);
    stats(nn).CrossCorrEmer = permute(band_amp_xcorr_emer, [3, 1, 2]);

    stats(nn).CrossCorrBSLPk = permute([band_amp_xcorr_bsl_cc; ...
        band_amp_xcorr_bsl_lag], [3, 1, 2]);
    stats(nn).CrossCorrAnesPk = permute([band_amp_xcorr_anes_cc; ...
        band_amp_xcorr_anes_lag], [3, 1, 2]);
    stats(nn).CrossCorrInterPk = permute([band_amp_xcorr_itr_cc; ...
        band_amp_xcorr_itr_lag], [3, 1, 2]);
    stats(nn).CrossCorrBurstPk = permute([band_amp_xcorr_burst_cc; ...
        band_amp_xcorr_burst_lag], [3, 1, 2]);
    stats(nn).CrossCorrEmerPk = permute([band_amp_xcorr_emer_cc; ...
        band_amp_xcorr_emer_lag], [3, 1, 2]);

    if any(strcmp(var_list, 'ca2_sig'))
        do_ca2 = true;
        load(fullfile(outpaths{nn}, 'outputs.mat'), 'ca2_fxx', ...
             'ca2_cent_freq_burst', 'ca2_cent_freq_itr', ...
            'ca2_pxx_bsl', 'ca2_pxx_anes', 'ca2_pxx_burst', 'ca2_pxx_itr', 'ca2_pxx_emer', ...
            'ca2_pxx_bsl_per', 'ca2_pxx_anes_per', 'ca2_pxx_burst_per', 'ca2_pxx_itr_per', 'ca2_pxx_emer_per', ...
            'ca2_bp_bin', 'ca2_bp_bin_center', 'ca2_bp_bin_per',...
            'ca2_bp_bsl', 'ca2_bp_anes', 'ca2_bp_burst', 'ca2_bp_itr', 'ca2_bp_emer', ...
            'ca2_bp_bsl_per', 'ca2_bp_anes_per', 'ca2_bp_burst_per', 'ca2_bp_itr_per', 'ca2_bp_emer_per', ...
            'pets_delta_sig', 'pets_ca2_delta_sig', 'delta_ca2_lag_ts', ...
            'delta_ca2_xcorr_anes', 'delta_ca2_xcorr_bsl', 'delta_ca2_xcorr_emer', ...
            'delta_ca2_xcorr_burst', 'delta_ca2_xcorr_itr',...
            'delta_ca2_xcorr_anes_cc', 'delta_ca2_xcorr_anes_lag', ...
            'delta_ca2_xcorr_bsl_cc', 'delta_ca2_xcorr_bsl_lag', ...
            'delta_ca2_xcorr_burst_cc', 'delta_ca2_xcorr_burst_lag', ...
            'delta_ca2_xcorr_itr_cc', 'delta_ca2_xcorr_itr_lag', ...
            'delta_ca2_xcorr_emer_cc', 'delta_ca2_xcorr_emer_lag', ...
            'delta_burst_ca2_peak_dt_med', 'delta_burst_ca2_peak_dt_anes_med');

        stats(nn).DtMed = delta_burst_ca2_peak_dt_med;
        stats(nn).DtMedAnes = delta_burst_ca2_peak_dt_anes_med;

        stats(nn).Ca2SpectrumBSL = ca2_pxx_bsl;
        stats(nn).Ca2SpectrumInter = ca2_pxx_itr;
        stats(nn).Ca2SpectrumBurst = ca2_pxx_burst;
        stats(nn).Ca2SpectrumAnes = ca2_pxx_anes;
        stats(nn).Ca2SpectrumEmer = ca2_pxx_emer;
        stats(nn).Ca2SpectrumInterCh = ca2_pxx_itr./ca2_pxx_bsl-1;
        stats(nn).Ca2SpectrumBurstCh = ca2_pxx_burst./ca2_pxx_bsl-1;
        stats(nn).Ca2SpectrumAnesCh = ca2_pxx_anes./ca2_pxx_bsl-1;
        stats(nn).Ca2SpectrumEmerCh = ca2_pxx_emer./ca2_pxx_bsl-1;
        stats(nn).Ca2CentFreqBurst = ca2_cent_freq_burst;
        stats(nn).Ca2CentFreqInter = ca2_cent_freq_itr;

        stats(nn).Ca2SpectrumBSLPer = ca2_pxx_bsl_per;
        stats(nn).Ca2SpectrumInterPer = ca2_pxx_itr_per;
        stats(nn).Ca2SpectrumBurstPer = ca2_pxx_burst_per;
        stats(nn).Ca2SpectrumAnesPer = ca2_pxx_anes_per;
        stats(nn).Ca2SpectrumEmerPer = ca2_pxx_emer_per;

        stats(nn).Ca2SpectrumInterChPer = ca2_pxx_itr_per./ca2_pxx_bsl_per-1;
        stats(nn).Ca2SpectrumBurstChPer = ca2_pxx_burst_per./ca2_pxx_bsl_per-1;
        stats(nn).Ca2SpectrumAnesChPer = ca2_pxx_anes_per./ca2_pxx_bsl_per-1;
        stats(nn).Ca2SpectrumEmerChPer = ca2_pxx_emer_per./ca2_pxx_bsl_per-1;

        stats(nn).Ca2BandpowerBin = ca2_bp_bin;
        stats(nn).Ca2BandpowerBinPer = ca2_bp_bin_per;
        stats(nn).Ca2BandpowerBinCnt = ca2_bp_bin_center;

        stats(nn).Ca2BandpowerBSL = ca2_bp_bsl;
        stats(nn).Ca2BandpowerAnes = ca2_bp_anes;
        stats(nn).Ca2BandpowerInter = ca2_bp_itr;
        stats(nn).Ca2BandpowerBurst = ca2_bp_burst;
        stats(nn).Ca2BandpowerEmer = ca2_bp_emer;

        stats(nn).Ca2BandpowerBSLPer = ca2_bp_bsl_per;
        stats(nn).Ca2BandpowerAnesPer = ca2_bp_anes_per;
        stats(nn).Ca2BandpowerInterPer = ca2_bp_itr_per;
        stats(nn).Ca2BandpowerBurstPer = ca2_bp_burst_per;
        stats(nn).Ca2BandpowerEmerPer = ca2_bp_emer_per;

        stats(nn).PSTHDelta = mean(pets_delta_sig, 1);
        stats(nn).PSTHCa2 = mean(pets_ca2_delta_sig, 1);

        stats(nn).CrossCorrCa2Anes = delta_ca2_xcorr_anes';
        stats(nn).CrossCorrCa2BSL = delta_ca2_xcorr_bsl';
        stats(nn).CrossCorrCa2Inter = delta_ca2_xcorr_itr';
        stats(nn).CrossCorrCa2Burst = delta_ca2_xcorr_burst';
        stats(nn).CrossCorrCa2Emer = delta_ca2_xcorr_emer';

        stats(nn).CrossCorrCa2AnesPk = ...
            [delta_ca2_xcorr_anes_cc, delta_ca2_xcorr_anes_lag];
        stats(nn).CrossCorrCa2BSLPk = ...
            [delta_ca2_xcorr_bsl_cc, delta_ca2_xcorr_bsl_lag];
        stats(nn).CrossCorrCa2BurstPk = ...
            [delta_ca2_xcorr_burst_cc, delta_ca2_xcorr_burst_lag];
        stats(nn).CrossCorrCa2InterPk = ...
            [delta_ca2_xcorr_itr_cc, delta_ca2_xcorr_itr_lag];
        stats(nn).CrossCorrCa2EmerPk = ...
            [delta_ca2_xcorr_emer_cc, delta_ca2_xcorr_emer_lag];
    end

    if any(strcmp(var_list, 'spike_ts'))
        do_spk = true;
        load(fullfile(outpaths{nn}, 'outputs.mat'), 'spike_ts', 'spike_unit', ...
            'spike_burst_vec', 'lfp_eeg_delta_dt', 'lfp_eeg_delta_dphi', ...
            'spike_lfp_rzpmk_bsl', 'spike_eeg_rzpmk_bsl', ...
            'spike_lfp_rzpmk_burst', 'spike_eeg_rzpmk_burst', ...
            'spike_lfp_rzpmk_itr', 'spike_eeg_rzpmk_itr', ...
            'spike_lfp_rzpmk_emer', 'spike_eeg_rzpmk_emer', ...
            'spike_led_psth', 'burst_ratio_per_unit', 'burst_ratio_per_delta_burst', ...
            'spike_burst_idx_bsl', 'spike_burst_idx_burst', ...
            'spike_burst_idx_itr', 'spike_burst_idx_emer', ...
            'sync_jbsi_bsl', 'sync_jbsi_itr', 'sync_jbsi_burst', 'sync_jbsi_emer', ...
            'ccg_t_lag', ...
            'ccg_zscore_bsl', 'ccg_zscore_burst', 'ccg_zscore_itr', 'ccg_zscore_emer', ...
            'ccg_perc_bsl', 'ccg_perc_burst', 'ccg_perc_itr', 'ccg_perc_emer');

        stats(nn).LFPdT = cellfun(@(x)(mean(x)), lfp_eeg_delta_dt);
        stats(nn).LFPdPhi = cellfun(@(x)(mean(x)), lfp_eeg_delta_dphi);

        stats(nn).SpkLFPRaylBSL = cat(1, spike_lfp_rzpmk_bsl{:});
        stats(nn).SpkEEGRaylBSL = cat(1, spike_eeg_rzpmk_bsl{:});
        stats(nn).SpkLFPRaylBurst = cat(1, spike_lfp_rzpmk_burst{:});
        stats(nn).SpkEEGRaylBurst = cat(1, spike_eeg_rzpmk_burst{:});
        stats(nn).SpkLFPRaylInter = cat(1, spike_lfp_rzpmk_itr{:});
        stats(nn).SpkEEGRaylInter = cat(1, spike_eeg_rzpmk_itr{:});
        stats(nn).SpkLFPRaylEmer = cat(1, spike_lfp_rzpmk_emer{:});
        stats(nn).SpkEEGRaylEmer = cat(1, spike_eeg_rzpmk_emer{:});


        unit_spike_ts_led = arrayfun(@(x)(spike_ts(spike_unit == x) ...
            - ev_mrk.t_led(1)), unique(spike_unit), 'UniformOutput', false);
        stats(nn).SpkLEDTs = cellfun(@(x)(mask_time_range(x, ...
            minmax(opt.peri_led_window), x)), unit_spike_ts_led, ...
            'UniformOutput', false);
        stats(nn).SpkLEDBurst = cellfun(@(x, y)(mask_time_range(x(y), ...
            minmax(opt.peri_led_window), x(y))), unit_spike_ts_led, ...
            spike_burst_vec, 'UniformOutput', false);
        stats(nn).SpkLEDPSTH = permute(cat(3, spike_led_psth{:}), [3, 2, 1]);
        stats(nn).SpkLEDPSTHNorm = stats(nn).SpkLEDPSTH ./ ...
            cat(1, spike_burst_idx_bsl{:}).spike_fr;

        stats(nn).SpkBurstIdxBSL = cat(1, spike_burst_idx_bsl{:});
        stats(nn).SpkBurstIdxBurst = cat(1, spike_burst_idx_burst{:});
        stats(nn).SpkBurstIdxInter = cat(1, spike_burst_idx_itr{:});
        stats(nn).SpkBurstIdxEmer = cat(1, spike_burst_idx_emer{:});
        stats(nn).BurstRatioPerUnit = burst_ratio_per_unit(:);
        stats(nn).BurstRatioPerEv = burst_ratio_per_delta_burst(:);

        stats(nn).SyncJBSIBSL = sync_jbsi_bsl;
        stats(nn).SyncJBSIBurst = sync_jbsi_burst;
        stats(nn).SyncJBSIInter = sync_jbsi_itr;
        stats(nn).SyncJBSIEmer = sync_jbsi_emer;

        stats(nn).CCGZscoreBSL = mean(ccg_zscore_bsl, 1, 'omitnan');
        stats(nn).CCGZscoreBurst = mean(ccg_zscore_burst, 1, 'omitnan');
        stats(nn).CCGZscoreInter = mean(ccg_zscore_itr, 1, 'omitnan');
        stats(nn).CCGZscoreEmer = mean(ccg_zscore_emer, 1, 'omitnan');

        stats(nn).CCGZLag0BSL = ccg_zscore_bsl( ...
            eye(sqrt(size(ccg_zscore_bsl, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGZLag0Burst = ccg_zscore_burst( ...
            eye(sqrt(size(ccg_zscore_burst, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGZLag0Inter = ccg_zscore_itr( ...
            eye(sqrt(size(ccg_zscore_itr, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGZLag0Emer = ccg_zscore_emer( ...
            eye(sqrt(size(ccg_zscore_emer, 1))) ~= 1, ccg_t_lag == 0);

        stats(nn).CCGPercLag0BSL = ccg_perc_bsl( ...
            eye(sqrt(size(ccg_perc_bsl, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGPercLag0Burst = ccg_perc_burst( ...
            eye(sqrt(size(ccg_perc_burst, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGPercLag0Inter = ccg_perc_itr( ...
            eye(sqrt(size(ccg_perc_itr, 1))) ~= 1, ccg_t_lag == 0);
        stats(nn).CCGPercLag0Emer = ccg_perc_emer( ...
            eye(sqrt(size(ccg_perc_emer, 1))) ~= 1, ccg_t_lag == 0);
    end
end

% do significance test
[spectrum_ch_burst_pval, spectrum_ch_burst_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumBurstCh', true), 10000, 0.01);
[spectrum_ch_inter_pval, spectrum_ch_inter_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumInterCh', true), 10000, 0.01);
[spectrum_ch_anes_pval, spectrum_ch_anes_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumAnesCh', true), 10000, 0.01);
[spectrum_ch_emer_pval, spectrum_ch_emer_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumEmerCh', true), 10000, 0.01);
[spectrum_ratio_ch_burst_pval, spectrum_ratio_ch_burst_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumBurstChPer', true), 10000, 0.01);
[spectrum_ratio_ch_inter_pval, spectrum_ratio_ch_inter_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumInterChPer', true), 10000, 0.01);
[spectrum_ratio_ch_anes_pval, spectrum_ratio_ch_anes_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumAnesChPer', true), 10000, 0.01);
[spectrum_ratio_ch_emer_pval, spectrum_ratio_ch_emer_ci] = bootstrap_1d( ...
    get_mean_sem_matrix(stats, 'SpectrumEmerChPer', true), 10000, 0.01);

if do_ca2
    [ca2_spectrum_ch_burst_pval, ca2_spectrum_ch_burst_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumBurstCh', true), 10000, 0.01);
    [ca2_spectrum_ch_inter_pval, ca2_spectrum_ch_inter_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumInterCh', true), 10000, 0.01);
    [ca2_spectrum_ch_anes_pval, ca2_spectrum_ch_anes_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumAnesCh', true), 10000, 0.01);
    [ca2_spectrum_ch_emer_pval, ca2_spectrum_ch_emer_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumEmerCh', true), 10000, 0.01);
    [ca2_spectrum_ratio_ch_burst_pval, ca2_spectrum_ratio_ch_burst_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumBurstChPer', true), 10000, 0.01);
    [ca2_spectrum_ratio_ch_inter_pval, ca2_spectrum_ratio_ch_inter_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumInterChPer', true), 10000, 0.01);
    [ca2_spectrum_ratio_ch_anes_pval, ca2_spectrum_ratio_ch_anes_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumAnesChPer', true), 10000, 0.01);
    [ca2_spectrum_ratio_ch_emer_pval, ca2_spectrum_ratio_ch_emer_ci] = bootstrap_1d( ...
        get_mean_sem_matrix(stats, 'Ca2SpectrumEmerChPer', true), 10000, 0.01);
end

fprintf('%gs.', toc(t));

%% save results
t = tic;
fprintf('\n\tsaving ... ');
outpath = fullfile(workpath, 'outputs.mat');
save(outpath, 'outpaths', 'stats', 'opt', 'spectrum_fxx', 'band_amp_lag_ts', ...
    'spectrum_ch_burst_pval', 'spectrum_ch_burst_ci', ...
    'spectrum_ch_inter_pval', 'spectrum_ch_inter_ci', ...
    'spectrum_ch_anes_pval', 'spectrum_ch_anes_ci', ...
    'spectrum_ch_emer_pval', 'spectrum_ch_emer_ci', ...
    'spectrum_ratio_ch_burst_pval', 'spectrum_ratio_ch_burst_ci', ...
    'spectrum_ratio_ch_inter_pval', 'spectrum_ratio_ch_inter_ci', ...
    'spectrum_ratio_ch_anes_pval', 'spectrum_ratio_ch_anes_ci', ...
    'spectrum_ratio_ch_emer_pval', 'spectrum_ratio_ch_emer_ci', ...
    '-v7.3');

if do_ca2
    save(outpath, 'ca2_fxx', 'delta_ca2_lag_ts', ...
        'ca2_spectrum_ch_burst_pval', 'ca2_spectrum_ch_burst_ci', ...
        'ca2_spectrum_ch_inter_pval', 'ca2_spectrum_ch_inter_ci', ...
        'ca2_spectrum_ch_anes_pval', 'ca2_spectrum_ch_anes_ci', ...
        'ca2_spectrum_ch_emer_pval', 'ca2_spectrum_ch_emer_ci', ...
        'ca2_spectrum_ratio_ch_burst_pval', 'ca2_spectrum_ratio_ch_burst_ci', ...
        'ca2_spectrum_ratio_ch_inter_pval', 'ca2_spectrum_ratio_ch_inter_ci', ...
        'ca2_spectrum_ratio_ch_anes_pval', 'ca2_spectrum_ratio_ch_anes_ci', ...
        'ca2_spectrum_ratio_ch_emer_pval', 'ca2_spectrum_ratio_ch_emer_ci', ...
        '-append');
end

if do_spk
    save(outpath, 'ccg_t_lag', '-append');
end

fprintf('%gs.', toc(t));

t = tic;
fprintf('\n\tplotting ... ');

%% delta burst duration
[cmean, csem, cmat] = get_mean_sem_matrix(stats, 'Duration');
figure; set(gcf, 'Position', [100, 100, 500, 450], 'Visible', false, ...
    'PaperOrientation', 'landscape'); hold on;
plot(cmat(:, 1), cmat(:, 2), 'o', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'r', 'MarkerSize', 6);
plot([cmean(1)-csem(1), cmean(1)+csem(1)], [cmean(2), cmean(2)], '-k');
plot([cmean(1), cmean(1)], [cmean(2)-csem(2), cmean(2)+csem(2)], '-k');
plot(get(gca, 'XLim'), get(gca, 'XLim'), '--k');
plot(cmean(1), cmean(2), 'or', 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'w', 'MarkerSize', 10);
xlabel('Log-normal \mu (s)');
ylabel('Median (s)');
title('Burst duration');
save_multi_formats(gcf, fullfile(workpath, 'stats_burst_duration'), ...
    opt.fig_save_format);
close;

%% delta burst interval
[cmean, csem, cmat] = get_mean_sem_matrix(stats, 'Interval');
figure; set(gcf, 'Position', [100, 100, 500, 450], 'Visible', false, ...
    'PaperOrientation', 'landscape'); hold on;
plot(cmat(:, 1), cmat(:, 2), 'o', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'b', 'MarkerSize', 6);
plot([cmean(1)-csem(1), cmean(1)+csem(1)], [cmean(2), cmean(2)], '-k');
plot([cmean(1), cmean(1)], [cmean(2)-csem(2), cmean(2)+csem(2)], '-k');
plot(get(gca, 'XLim'), get(gca, 'XLim'), '--k');
plot(cmean(1), cmean(2), 'or', 'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'w', 'MarkerSize', 10);
xlabel('Log-normal \mu (s)');
ylabel('Median (s)');
title('Burst interval');
save_multi_formats(gcf, fullfile(workpath, 'stats_burst_interval'), ...
    opt.fig_save_format);
close;

%% periodogram
figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape');
plot_grouped_periodogram(cellfun(@(x)(10*log10(get_mean_sem_matrix(stats, x, true))), ...
    {'SpectrumBSL', 'SpectrumInter', 'SpectrumBurst', 'SpectrumEmer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, nan(size(spectrum_fxx)), spectrum_ch_inter_pval(1, :), ...
    spectrum_ch_burst_pval(1, :), spectrum_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'PSD (dB/Hz)', ...
    {'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
    'EEG periodogram', {'k', 'b', 'r', 'm'});
save_multi_formats(gcf, fullfile(workpath, 'stats_periodogram_burst'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(10*log10(get_mean_sem_matrix(stats, x, true))), ...
    {'SpectrumBSL', 'SpectrumAnes', 'SpectrumEmer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, nan(size(spectrum_fxx)), ...
    spectrum_ch_anes_pval(1, :), spectrum_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'PSD (dB/Hz)', ...
    {'Pre-injection', 'Anesthesia', 'Emergence'}, ...
    'EEG periodogram', {'k', 'b', 'r'});
save_multi_formats(gcf, fullfile(workpath, 'stats_periodogram_anes'), ...
    opt.fig_save_format);
close;

%% relative periodogram
figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumBSLPer', 'SpectrumInterPer', 'SpectrumBurstPer', 'SpectrumEmerPer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, nan(size(spectrum_fxx)), spectrum_ratio_ch_inter_pval(1, :), ...
    spectrum_ratio_ch_burst_pval(1, :), spectrum_ratio_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Relative PSD (%)', ...
    {'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
    'EEG periodogram', {'k', 'b', 'r', 'm'});
save_multi_formats(gcf, fullfile(workpath, 'stats_relative_periodogram_burst'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumBSLPer', 'SpectrumAnesPer', 'SpectrumEmerPer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, nan(size(spectrum_fxx)), ...
    spectrum_ratio_ch_anes_pval(1, :), spectrum_ratio_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Relative PSD (%)', ...
    {'Pre-injection', 'Anesthesia', 'Emergence'}, ...
    'EEG periodogram', {'k', 'b', 'r'});
save_multi_formats(gcf, fullfile(workpath, 'stats_relative_periodogram_anes'), ...
    opt.fig_save_format);
close;

%% periodogram of power change
figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumInterCh', 'SpectrumBurstCh', 'SpectrumEmerCh'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, spectrum_ch_inter_pval(1, :), ...
    spectrum_ch_burst_pval(1, :), spectrum_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Power change ratio', ...
    {'Inter-burst', 'Delta burst', 'Emergence'}, 'EEG periodogram', {'b', 'r', 'm'});
save_multi_formats(gcf, fullfile(workpath, 'stats_periodogram_change_burst'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumAnesCh', 'SpectrumEmerCh'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, spectrum_ch_anes_pval(1, :), spectrum_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Power change ratio', ...
    {'Anesthesia', 'Emergence'}, 'EEG periodogram', {'b', 'r'});
save_multi_formats(gcf, fullfile(workpath, 'stats_periodogram_change_anes'), ...
    opt.fig_save_format);
close;

%%  Relative periodogram of power change
figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumInterChPer', 'SpectrumBurstChPer', 'SpectrumEmerChPer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, spectrum_ratio_ch_inter_pval(1, :), ...
    spectrum_ratio_ch_burst_pval(1, :), spectrum_ratio_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Relative Power change ratio', ...
    {'Inter-burst', 'Delta burst', 'Emergence'}, 'EEG periodogram', {'b', 'r', 'm'});
save_multi_formats(gcf, fullfile(workpath, 'stats_relative_periodogram_change_burst'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); 
plot_grouped_periodogram(cellfun(@(x)(get_mean_sem_matrix(stats, x, true)), ...
    {'SpectrumAnesChPer', 'SpectrumEmerChPer'}, 'UniformOutput', false), ...
    spectrum_fxx, cat(1, spectrum_ratio_ch_anes_pval(1, :), spectrum_ratio_ch_emer_pval(1, :)), ...
    0.05, [0, 100], unique(cat(2, opt.band_freq{:})), true, 'Relative Power change ratio', ...
    {'Anesthesia', 'Emergence'}, 'EEG periodogram', {'b', 'r'});
save_multi_formats(gcf, fullfile(workpath, 'stats_relative_periodogram_change_anes'), ...
    opt.fig_save_format);
close;

%% band amplitude timecourse
clrs = get_default_colors;
band_label_other = opt.band_label(~strcmpi(opt.band_label, 'delta'));
num_band = length(band_label_other);
[cmean0, csem0] = get_mean_sem_matrix(stats, 'PSTHDeltaAmp');
[cmean, csem] = get_mean_sem_matrix(stats, 'PSTHOtherAmp');
for bb = 1:num_band
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;

    yyaxis(gca, 'left');
    shadedErrorBar(opt.peri_event_window, cmean0, csem0, ...
        'LineProps', {'-', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1});
    ylabel('Amplitude of delta (\muV)');
    set(gca, 'YColor', [0.5, 0.5, 0.5]);

    yyaxis(gca, 'right');
    shadedErrorBar(opt.peri_event_window, cmean(:, :, bb), csem(:, :, bb), ...
        'LineProps', {'-', 'Color', clrs(bb, :), 'LineWidth', 1});
    xline(0, '--k');
    set(gca, 'YColor', clrs(bb, :));
    xlabel('Time to delta peak (s)');
    ylabel(sprintf('Amplitude of %s (\\muV)', band_label_other{bb}));

    figname = ['stats_band_amplitude_timecourse_', band_label_other{bb}];
    save_multi_formats(gcf, fullfile(workpath, figname), opt.fig_save_format);
    close;
end

%% cross-correlation
clrs = get_default_colors;
band_label_other = opt.band_label(~strcmpi(opt.band_label, 'delta'));
num_band = length(band_label_other);
[cmean0, csem0, cmat0] = get_mean_sem_matrix(stats, 'CrossCorrBSL');
[cmean1, csem1, cmat1] = get_mean_sem_matrix(stats, 'CrossCorrAnes');
[cmean2, csem2, cmat2] = get_mean_sem_matrix(stats, 'CrossCorrInter');
[cmean3, csem3, cmat3] = get_mean_sem_matrix(stats, 'CrossCorrBurst');
[cmean4, csem4, cmat4] = get_mean_sem_matrix(stats, 'CrossCorrEmer');


for bb = 1:num_band
    figure; set(gcf, 'Position', [100, 50, 1000, 400], 'Visible', false, ...
        'PaperOrientation', 'landscape');

    subplot(1, 2, 1); hold on;
    shadedErrorBar(band_amp_lag_ts*1000, cmean0(:, :, bb), csem0(:, :, bb), ...
        'lineProps', {'--', 'Color', clrs(bb, :)});
    shadedErrorBar(band_amp_lag_ts*1000, cmean1(:, :, bb), csem1(:, :, bb), ...
        'lineProps', {'-', 'Color', clrs(bb, :), 'LineWidth', 2});
    shadedErrorBar(band_amp_lag_ts*1000, cmean4(:, :, bb), csem4(:, :, bb), ...
        'lineProps', {'--', 'Color', clrs(bb, :)});

    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off', 'Location', 'best');
    title(sprintf('%s refer to delta', band_label_other{bb}));

    subplot(1, 2, 2); hold on;
    cmat = zeros(num_out, 2, 2);
    for nn = 1:num_out
        [~, imax] = max(abs(cmat1(nn, :, bb)), [], 2);
        cmat(nn, 1, 1) = cmat1(nn, imax, bb);
        cmat(nn, 1, 2) = -band_amp_lag_ts(imax)*1000;
        [~, imax] = max(abs(cmat0(nn, :, bb)), [], 2);
        cmat(nn, 2, 1) = cmat0(nn, imax, bb);
        cmat(nn, 2, 2) = -band_amp_lag_ts(imax)*1000;
        plot(cmat(nn, :, 2), cmat(nn, :, 1), '-k', 'Color', [.5, .5, .5]);
        plot(cmat(nn, 2, 2), cmat(nn, 2, 1), 'o', 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', [.5, .5, .5], 'LineWidth', 2, 'MarkerSize', 8);
        plot(cmat(nn, 1, 2), cmat(nn, 1, 1), 'o', 'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', clrs(bb, :), 'LineWidth', 2, 'MarkerSize', 8);
    end
    cmean = median(cmat);
    csem = prctile(cmat, [25, 75], 1);
    cp1 = signrank(cmat(:, 1, 1));
    cp2 = signrank(cmat(:, 1, 2));
    plot(csem(:, 2, 2), [1, 1]*cmean(1, 2, 1), '-k', 'LineWidth', 1.5);
    plot([1, 1]*cmean(1, 2, 2), csem(:, 2, 1), '-k', 'LineWidth', 1.5);
    plot(csem(:, 1, 2), [1, 1]*cmean(1, 1, 1), '-k', 'LineWidth', 1.5);
    plot([1, 1]*cmean(1, 1, 2), csem(:, 1, 1), '-k', 'LineWidth', 1.5);
    plot(cmean(1, 2, 2), cmean(1, 2, 1), 'o', 'MarkerFaceColor', [.5, .5, .5], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 12);
    plot(cmean(1, 1, 2), cmean(1, 1, 1), 'o', 'MarkerFaceColor', clrs(bb, :), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 12);
    set(gca, 'TickDir', 'out');
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title(sprintf('lag median: %.2f, \\itp\\rm=%.3f\n\\itr\\rm median: %.2f, \\itp\\rm=%.3f', ...
        cmean(1, 1, 2), cp2, cmean(1, 1, 1), cp1));

    figname = ['stats_crosscorr_anes_', band_label_other{bb}];
    save_multi_formats(gcf, fullfile(workpath, figname), opt.fig_save_format);
    close;


    figure; set(gcf, 'Position', [100, 50, 500, 400], 'Visible', false, ...
        'PaperOrientation', 'landscape'); hold on;

    shadedErrorBar(band_amp_lag_ts*1000, cmean0(:, :, bb), csem0(:, :, bb), ...
        'lineProps', {'--', 'Color', clrs(bb, :)});
    shadedErrorBar(band_amp_lag_ts*1000, cmean2(:, :, bb), csem2(:, :, bb), ...
        'lineProps', {':', 'Color', clrs(bb, :), 'LineWidth', 1});
    shadedErrorBar(band_amp_lag_ts*1000, cmean3(:, :, bb), csem3(:, :, bb), ...
        'lineProps', {'-', 'Color', clrs(bb, :), 'LineWidth', 2});
    shadedErrorBar(band_amp_lag_ts*1000, cmean4(:, :, bb), csem4(:, :, bb), ...
        'lineProps', {'-', 'Color', clrs(bb, :), 'LineWidth', 2});

    cp = ttest(cmat3(:, :, bb)-cmat2(:, :, bb)); cp(cp==0) = nan;
    plot(band_amp_lag_ts*1000, (1-0.02*3)*max(get(gca, 'YLim'))*cp, ...
        '.m', 'LineWidth', 1.5);
    cp = ttest(cmat2(:, :, bb)-cmat0(:, :, bb)); cp(cp==0) = nan;
    plot(band_amp_lag_ts*1000, (1-0.02*2)*max(get(gca, 'YLim'))*cp, ...
        '.b', 'LineWidth', 1.5);
    cp = ttest(cmat3(:, :, bb)-cmat0(:, :, bb)); cp(cp==0) = nan;
    plot(band_amp_lag_ts*1000, (1-0.02*1)*max(get(gca, 'YLim'))*cp, ...
        '.r', 'LineWidth', 1.5);
    cp = ttest(cmat4(:, :, bb)-cmat0(:, :, bb)); cp(cp==0) = nan;
    plot(band_amp_lag_ts*1000, (1-0.02*1)*max(get(gca, 'YLim'))*cp, ...
        '.k', 'LineWidth', 1.5);

    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
        'Box', 'off', 'Location', 'best');
    title(sprintf('%s refer to delta', band_label_other{bb}));

    figname = ['stats_crosscorr_burst_', band_label_other{bb}];
    save_multi_formats(gcf, fullfile(workpath, figname), opt.fig_save_format);
    close;
end

%% calcium periodogram
if do_ca2
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat0] = get_mean_sem_matrix(stats, 'Ca2SpectrumBSL');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'k'});
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumInter');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumBurst');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});
    [cmean, csem, cmat3] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'm'});

    cp = ttest(cmat1./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 2*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);
    cp = ttest(cmat3./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 3*0.02*max(get(gca, 'YLim'))*cp, '.m', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (a.u.)');
    legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, 'stats_calcium_periodogram_burst'), ...
        opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat0] = get_mean_sem_matrix(stats, 'Ca2SpectrumBSL');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'k'});
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumAnes');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});

    cp = ttest(cmat1./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (a.u.)');
    legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, 'stats_calcium_periodogram_anes'), ...
        opt.fig_save_format);
    close;
end

%%  relative calcium periodogram
if do_ca2
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat0] = get_mean_sem_matrix(stats, 'Ca2SpectrumBSLPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'k'});
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumInterPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumBurstPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});
    [cmean, csem, cmat3] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'm'});

    cp = ttest(cmat1./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 2*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);
    cp = ttest(cmat3./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 3*0.02*max(get(gca, 'YLim'))*cp, '.m', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Relative Amplitude (%)');
    legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, 'stats_relative_calcium_periodogram_burst'), ...
        opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat0] = get_mean_sem_matrix(stats, 'Ca2SpectrumBSLPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'k'});
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumAnesPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});

    cp = ttest(cmat1./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2./cmat0-1); cp(cp==0) = nan;
    plot(ca2_fxx, 1*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Relative Amplitude (%)');
    legend({'Pre-injection', 'Anesthesia', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, 'stats_relative_calcium_periodogram_anes'), ...
        opt.fig_save_format);
    close;
end

%% calcium periodogram of power change
if do_ca2
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumInterCh');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumBurstCh');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});
    [cmean, csem, cmat3] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerCh');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'm'});


    cp = ttest(cmat2-cmat1); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat3-cmat2); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.m', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude change ratio');
    legend({'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_calcium_periodogram_change_burst'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumAnesCh');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerCh');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});

    cp = ttest(cmat1); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude change ratio');
    legend({'Anesthesia', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_calcium_periodogram_change_anes'), opt.fig_save_format);
    close;
end


%% relative calcium periodogram of power change
if do_ca2
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumInterChPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumBurstChPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});
    [cmean, csem, cmat3] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerChPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'm'});


    cp = ttest(cmat2-cmat1); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat3-cmat2); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.m', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Relative Amplitude change ratio');
    legend({'Inter-burst', 'Delta burst', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_relative_calcium_periodogram_change_burst'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;
    [cmean, csem, cmat1] = get_mean_sem_matrix(stats, 'Ca2SpectrumAnesChPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'b'});
    [cmean, csem, cmat2] = get_mean_sem_matrix(stats, 'Ca2SpectrumEmerChPer');
    shadedErrorBar(ca2_fxx, cmean, csem, 'lineProps', {'Color', 'r'});

    cp = ttest(cmat1); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.b', 'LineWidth', 1.5);
    cp = ttest(cmat2); cp(cp==0) = nan;
    plot(ca2_fxx, -3*0.02*max(get(gca, 'YLim'))*cp, '.r', 'LineWidth', 1.5);

    set(gca, 'XLim', [0.01, 5], 'XTick', unique(cat(2, opt.ca2_band_freq{:})), ...
        'XGrid', 'on', 'XMinorGrid', 'off', 'Box', 'off');
    xlabel('Frequency (Hz)');
    ylabel('Relative Amplitude change ratio');
    legend({'Anesthesia', 'Emergence'}, 'Box', 'off');
    title('Wavelet transform of calcium signal');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_relative_calcium_periodogram_change_anes'), opt.fig_save_format);
    close;
end

%% calcium timecourse
if do_ca2
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;

    yyaxis(gca, 'left');
    [cmean, csem] = get_mean_sem_matrix(stats, 'PSTHDelta');
    shadedErrorBar(opt.peri_event_window, cmean, csem, ...
        'LineProps', {'-', 'Color', [0.3, 0.5, 1], 'LineWidth', 1});
    ylabel('Delta signal (\muV)');
    set(gca, 'YLim', [-350, 350], 'YColor', [0.3, 0.5, 0.1]);

    yyaxis(gca, 'right');
    [cmean, csem] = get_mean_sem_matrix(stats, 'PSTHCa2');
    shadedErrorBar(opt.peri_event_window, cmean, csem, ...
        'LineProps', {'-', 'Color', [0, 0.7, 0.1], 'LineWidth', 1});
    plot([0, 0], get(gca, 'YLim'), '--k');
    set(gca, 'YColor', [0, 0.7, 0.1]);
    xlabel('Time to delta peak (s)');
    ylabel('\DeltaF/F');

    save_multi_formats(gcf, fullfile(workpath, 'stats_calcium_timecourse'), ...
        opt.fig_save_format);
    close;
end

%% delta calcium cross-correlation
if do_ca2
    clrs = get_default_colors;
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'landscape'); hold on;

    [cmean0, csem0, cmat0] = get_mean_sem_matrix(stats, 'CrossCorrCa2BSL');
    [cmean1, csem1, cmat1] = get_mean_sem_matrix(stats, 'CrossCorrCa2Anes');
    [cmean2, csem2, cmat2] = get_mean_sem_matrix(stats, 'CrossCorrCa2Inter');
    [cmean3, csem3, cmat3] = get_mean_sem_matrix(stats, 'CrossCorrCa2Burst');
    [cmean4, csem4, cmat4] = get_mean_sem_matrix(stats, 'CrossCorrCa2Emer');

    subplot(2, 3, 1); hold on;
    plot(delta_ca2_lag_ts*1000, cmat0, 'Color', [0.8, 0.8, 0.8]);
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean0, csem0, ...
        'LineProps', {'-', 'Color', clrs(1, :), 'LineWidth', 1.5});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('Pre-injection');

    subplot(2, 3, 2); hold on;
    plot(delta_ca2_lag_ts*1000, cmat1, 'Color', [0.8, 0.8, 0.8]);
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean1, csem1, ...
        'LineProps', {'-', 'Color', clrs(2, :), 'LineWidth', 1.5});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('anesthesia');

    subplot(2, 3, 3); hold on;
    plot(delta_ca2_lag_ts*1000, cmat2, 'Color', [0.8, 0.8, 0.8]);
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean2, csem2, ...
        'LineProps', {'-', 'Color', clrs(3, :), 'LineWidth', 1.5});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('inter-burst');

    subplot(2, 3, 4); hold on;
    plot(delta_ca2_lag_ts*1000, cmat3, 'Color', [0.8, 0.8, 0.8]);
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean3, csem3, ...
        'LineProps', {'-', 'Color', clrs(4, :), 'LineWidth', 1.5});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('delta-burst');

    subplot(2, 3, 5); hold on;
    plot(delta_ca2_lag_ts*1000, cmat4, 'Color', [0.8, 0.8, 0.8]);
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean4, csem4, ...
        'LineProps', {'-', 'Color', clrs(5, :), 'LineWidth', 1.5});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('Emergence');
    save_multi_formats(gcf, fullfile(workpath, 'stats_ca2_xcorr_1'), ...
        opt.fig_save_format);
    close;


    figure; set(gcf, 'Position', [100, 50, 1000, 400], 'Visible', false, ...
        'PaperOrientation', 'landscape');

    subplot(1, 2, 1); hold on;
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean0, csem0, ...
        'LineProps', {'--', 'Color', clrs(1, :), 'LineWidth', 1});
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean1, csem1, ...
        'LineProps', {'-', 'Color', clrs(2, :), 'LineWidth', 1});
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean4, csem4, ...
        'LineProps', {'--', 'Color', clrs(5, :), 'LineWidth', 1});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('Calcium signal ref. to delta');
    legend({'Pre-injection', 'anesthesia', 'emergence'}, 'Box', 'off', 'Location', 'best');

    subplot(1, 2, 2); hold on;
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean0, csem0, ...
        'LineProps', {'--', 'Color', clrs(1, :), 'LineWidth', 1});
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean2, csem2, ...
        'LineProps', {':', 'Color', clrs(3, :), 'LineWidth', 1});
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean3, csem3, ...
        'LineProps', {'-', 'Color', clrs(4, :), 'LineWidth', 1});
    shadedErrorBar(delta_ca2_lag_ts*1000, cmean4, csem4, ...
        'LineProps', {'--', 'Color', clrs(5, :), 'LineWidth', 1});
    set(gca, 'TickDir', 'out', 'XLim', [-1, 1]*opt.delta_ccg_lag*1000, 'YLim', ...
        [min(min(get(gca, 'YLim')), -0.1), max(get(gca, 'YLim'))]);
    xlabel('Lag (ms)');
    ylabel('Cross-correlation coefficient');
    title('Calcium signal ref. to delta');
    legend({'Pre-injection', 'Inter-burst', 'Delta burst', 'Emergence'}, ...
        'Box', 'off', 'Location', 'best');
    save_multi_formats(gcf, fullfile(workpath, 'stats_ca2_xcorr_2'), ...
        opt.fig_save_format);
    close;
end

%% LFP-EEG time/phase shift
if do_spk
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_lfp_eeg_shift( ...
        cat(1, stats(:).LFPdT), cat(1, stats(:).LFPdPhi));
    title('LFP-EEG shift');
    save_multi_formats(gcf, fullfile(workpath, 'stats_lfp_eeg_shift'), ...
        opt.fig_save_format);
    close;
end

%% spike-LFP phase locking
if do_spk
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkLFPRaylBSL), 'Spike-LFP phase locking (Bsaeline)', 'c');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_lfp_phase_locking (Bsaeline)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkLFPRaylBSL), 'Spike-LFP phase locking (Bsaeline)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_lfp_phase_locking_dist (Bsaeline)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkLFPRaylBurst), 'Spike-LFP phase locking (Delta burst)', 'c');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_lfp_phase_locking (Delta burst)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkLFPRaylBurst), 'Spike-LFP phase locking (Delta burst)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_lfp_phase_locking_dist (Delta burst)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkLFPRaylInter), 'Spike-LFP phase locking (Inter burst)', 'c');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_lfp_phase_locking (Inter burst)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkLFPRaylInter), 'Spike-LFP phase locking (Inter burst)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_lfp_phase_locking_dist (Inter burst)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkLFPRaylEmer), 'Spike-LFP phase locking (Emergence)', 'c');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_lfp_phase_locking (Emergence)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkLFPRaylEmer), 'Spike-LFP phase locking (Emergence)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_lfp_phase_locking_dist (Emergence)'), opt.fig_save_format);
    close;
end

%% spike-EEG phase locking
if do_spk
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkEEGRaylBSL), 'Spike-EEG phase locking (Bsaeline)', 'm');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_eeg_phase_locking (Bsaeline)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkEEGRaylBSL), 'Spike-EEG phase locking (Bsaeline)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_eeg_phase_locking_dist (Bsaeline)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkEEGRaylBurst), 'Spike-EEG phase locking (Delta burst)', 'm');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_eeg_phase_locking (Delta burst)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkEEGRaylBurst), 'Spike-EEG phase locking (Delta burst)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_eeg_phase_locking_dist (Delta burst)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkEEGRaylInter), 'Spike-EEG phase locking (Inter burst)', 'm');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_eeg_phase_locking (Inter burst)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkEEGRaylInter), 'Spike-EEG phase locking (Inter burst)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_eeg_phase_locking_dist (Inter burst)'), opt.fig_save_format);
    close;


    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_spike_phase_locking( ...
        cat(1, stats(:).SpkEEGRaylEmer), 'Spike-EEG phase locking (Emergence)', 'm');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_spike_eeg_phase_locking (Emergence)'), opt.fig_save_format);
    close;

    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_phase_locking_dist( ...
        cat(1, stats(:).SpkEEGRaylEmer), 'Spike-EEG phase locking (Emergence)');
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_eeg_phase_locking_dist (Emergence)'), opt.fig_save_format);
    close;
end

%% firing rate & burst percentage
if do_spk
    figure; set(gcf, 'Position', [100, 100, 790, 790], 'Visible', false, ...
        'PaperOrientation', 'portrait');
    plot_summary_firing_burst(cat(1, stats(:).SpkBurstIdxBurst), ...
        cat(1, stats(:).SpkBurstIdxInter));
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_firing_rate_burst_percentage'), opt.fig_save_format);
    close;
end

%% raster plot peri-LED pulse
if do_spk
    bar_yv = [-0.3; 0.3];
    sm_sigma = 5;
    sm_kernel = normpdf((-3*sm_sigma:1:3*sm_sigma), 0, sm_sigma);
    figure; set(gcf, 'Position', [100, 100, 600, 800], 'Visible', false, ...
        'PaperOrientation', 'portrait');

    cmat1 = cat(1, stats(:).SpkLEDTs);
    cmat2 = cat(1, stats(:).SpkLEDBurst);
    num_unit = size(cmat, 1);
    subplot('Position', [0.15, 0.32, 0.70, 0.63]); hold on;
    plot([0, 0], [0, num_unit + 1], '--k');
    for ii = 1:num_unit
        if ~any(cmat1{ii}); continue; end
        plot(repmat(cmat1{ii}', 2, 1), cmat1{ii}'*0+ii+bar_yv, '-', ...
            'Color', [.5, .5, .5], 'LineWidth', 0.3);
        if ~any(cmat2{ii}); continue; end
        plot(repmat(cmat2{ii}', 2, 1), cmat2{ii}'*0+ii+bar_yv, '-b');
    end
    set(gca, 'Visible', 'off', 'YLim', [0, num_unit + 1], 'XLim', ...
        minmax(opt.peri_led_window));

    cmat = arrayfun(@(x)(x.SpkLEDPSTHNorm(:, :, 1)), stats, 'UniformOutput', false);
    cmat = conv2(cat(1, cmat{:}), sm_kernel, 'same');
    subplot('Position', [0.15, 0.10, 0.70, 0.20]); hold on;
    xline(0, '--k');
    shadedErrorBar((opt.peri_led_window(1:end-1) + opt.peri_led_window(2:end)) / 2, ...
        mean(cmat, 1, 'omitnan'), std(cmat, 0, 1, 'omitnan') / sqrt(num_unit), ...
        'lineProps', {'-k', 'LineWidth', 1});
    set(gca, 'TickDir', 'out', 'XLim', minmax(opt.peri_led_window), 'YLim', ...
        min(get(gca, 'YLim'), [0, inf]));
    xlabel('Time to LED onset (s)');
    ylabel('Averaged firing rate (Hz)');

    save_multi_formats(gcf, fullfile(workpath, 'stats_raster_plot_led_pulse'), ...
        opt.fig_save_format);
    close;
end

%% population cross-correlogram of unit spikes
if do_spk
    % workpath = 'G:\in vivo\sorting\普通tetrode\普通tetrode统计分析\test-cluster(1.5-4)';
    figure; set(gcf, 'Visible', false, 'PaperOrientation', 'portrait');
    plot_summary_population_ccg(ccg_t_lag, cat(1, stats(:).CCGZscoreBSL), ...
        cat(1, stats(:).CCGZscoreBurst), cat(1, stats(:).CCGZscoreInter), ...
        cat(1, stats(:).CCGZscoreEmer));
    save_multi_formats(gcf, fullfile(workpath, ...
        'stats_population_cross-correlogram'), opt.fig_save_format);
    close;
end

fprintf('%gs.', toc(t));
end

function [cmean, csem, cmat] = get_mean_sem_matrix(s, fd, flag)
if nargin < 3; flag = false; end
cmat = cat(1, s(:).(fd));
if flag; cmean = cmat; return; end
cmean = mean(cmat, 1, 'omitnan');
csem = std(cmat, 0, 1, 'omitnan')/sqrt(size(cmat, 1));
end

function [pval, ci] = bootstrap_1d(cmat, nrnd, alpha)
if nargin < 2 || isempty(nrnd); nrnd = 10000; end
if nargin < 3 || isempty(alpha); alpha = 0.05; end

[N, np] = size(cmat);
pval = nan(1, np);
ci = nan(2, np);

for pp = 1:np
    cx = cmat(:, pp);
    cp  = mean(cx(randi(N, N, nrnd)), 1, 'omitnan');
    pval(pp) = 1 - max(mean(cp < 0), mean(cp > 0));
    ci(:, pp) = prctile(cp, [alpha, 1-alpha] * 100);
end

padj = mafdr(pval, 'BHFDR', true);
pval = [pval; padj];
end
