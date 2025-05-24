function outpath = analyze_related_lfp_spike(lfp_sig, lfp_ts, lfp_fs, ...
    spike_ts, spike_cluster, cluster_ch, ch_grp, ex_ch, ex_cluster, opt, ...
    db_out, workpath)

t = tic;
fprintf('\n\tcalculating basic indexes ... ');
load(db_out, 'eeg_ts', 'ev_mrk', 'delta_sig', 'delta_burst_peak_t_anes', ...
    'delta_burst_start', 'delta_burst_stop');
if ~isfolder(workpath); mkdir(workpath); end

%% compute metrics
% period info
t_anes = ev_mrk.t_anes;
t_bsl = ev_mrk.t_bsl;
t_led = ev_mrk.t_led;
t_emer = ev_mrk.t_emer;

% delta (inter) burst information
delta_burst_center = (delta_burst_start + delta_burst_stop)/2;
delta_burst_start_anes = mask_time_range(delta_burst_center, t_anes, delta_burst_start);
delta_burst_stop_anes = mask_time_range(delta_burst_center, t_anes, delta_burst_stop);
delta_burst_anes = [delta_burst_start_anes, delta_burst_stop_anes];
inter_burst_anes = [delta_burst_stop_anes(1:end-1), delta_burst_start_anes(2:end)];

% screen & merge clusters / channels into units
[spike_ts, spike_cluster, cluster_ch, ch_grp] = exclude_channel_cluster( ...
    spike_ts, spike_cluster, cluster_ch, ch_grp, ex_ch, ex_cluster);
[spike_unit, unit_cluster, unit_ch, unit_grp] = merge_spike_to_unit( ...
    spike_cluster, cluster_ch, ch_grp, opt.spike_merge_by);

% create outputs related to included channels
num_channel = size(lfp_sig, 2);
lfp_delta = zeros(size(lfp_sig), 'like', lfp_sig);
[pets_lfp_sig, pets_led_lfp, lfp_eeg_delta_dt, lfp_eeg_delta_dphi] = ...
    deal(cell(num_channel, 1));

for jj = 1:num_channel
    if ismember(jj, ex_ch); continue; end
    
    % filter lfp signal
    [~, lfp_delta(:, jj)] = get_band_envelope(get_notched_sig( ...
        lfp_sig(:, jj), lfp_fs), lfp_fs, ...
        opt.band_freq(strcmpi(opt.band_label, 'delta')));
    
    % delta lfp peri-delta-peak-event signal 
    pets_lfp_sig{jj} = get_peri_event_signal(lfp_delta(:, jj), ...
        lfp_ts, delta_burst_peak_t_anes, opt.peri_event_window);
    
    % delta lfp peri-led-event signal 
    pets_led_lfp{jj} = get_peri_event_signal(lfp_delta(:, jj), ...
        lfp_ts, t_led, opt.peri_led_window);

    % lfp phase/time shift refer to eeg
    [lfp_eeg_delta_dt{jj}, lfp_eeg_delta_dphi{jj}] = get_peak_shift_t_phi( ...
        lfp_delta(:, jj), lfp_ts, delta_sig, eeg_ts, ...
        delta_burst_peak_t_anes);
end

% create outputs related to units
num_unit = length(unit_ch);
[spike_eeg_phase_bsl,spike_eeg_rzpmk_bsl, spike_lfp_phase_bsl, spike_lfp_rzpmk_bsl, ...
    spike_eeg_phase_burst,spike_eeg_rzpmk_burst, spike_lfp_phase_burst, spike_lfp_rzpmk_burst, ...
    spike_eeg_phase_itr,spike_eeg_rzpmk_itr, spike_lfp_phase_itr, spike_lfp_rzpmk_itr, ...
    spike_eeg_phase_emer,spike_eeg_rzpmk_emer, spike_lfp_phase_emer, spike_lfp_rzpmk_emer, ...
    spike_burst_vec, spike_burst_time, burst_per_delta_burst, ...
    spike_burst_idx_bsl, spike_burst_idx_burst, spike_burst_idx_itr, spike_burst_idx_emer, ...
    spike_delta_psth, spike_led_psth] = deal(cell(num_unit, 1));

for nn = 1:num_unit
    uspk = spike_ts(spike_unit==nn);
    
    % histogram of peri-delta-peak-event spikes
    spike_delta_psth{nn} = get_peri_event_hist(uspk, delta_burst_peak_t_anes, ...
        opt.peri_event_window);

    % histogram of peri-led-event spikes
    spike_led_psth{nn} = get_peri_event_hist(uspk, t_led, opt.peri_led_window);
    
    % detect burst sequences
    [spike_burst_vec{nn}, idx_beg, idx_end] = detect_spike_burst(uspk, ...
        opt.beginning_isi_ms, opt.end_isi_ms, opt.minimum_sequence, opt.hyperpolarization_ms);
    spike_burst_time{nn} = cat(2, uspk(idx_beg), uspk(idx_end));
    
    % compute indexes of spikes and burst spiking
    spike_burst_idx_bsl{nn} = compute_spiking_index(uspk, ...
        spike_burst_vec{nn}, idx_beg, idx_end, t_bsl);
    spike_burst_idx_burst{nn} = compute_spiking_index(uspk, ...
        spike_burst_vec{nn}, idx_beg, idx_end, delta_burst_anes);
    spike_burst_idx_itr{nn} = compute_spiking_index(uspk, ...
        spike_burst_vec{nn}, idx_beg, idx_end, inter_burst_anes);
    spike_burst_idx_emer{nn} = compute_spiking_index(uspk, ...
        spike_burst_vec{nn}, idx_beg, idx_end, t_emer);
    

    % compute number of bursts within each delta-burst event
    burst_per_delta_burst{nn} = arrayfun(@(x, y)(sum(mask_time_range( ...
        uspk(idx_beg), [x, y]))), delta_burst_start_anes, delta_burst_stop_anes);
       
    % spike-EEG phase locking within burst period
    [spike_eeg_phase_bsl{nn}, spike_eeg_rzpmk_bsl{nn}] = spike_phase_locking( ...
        uspk, delta_sig, eeg_ts, t_bsl);
    [spike_eeg_phase_burst{nn}, spike_eeg_rzpmk_burst{nn}] = spike_phase_locking( ...
        uspk, delta_sig, eeg_ts, delta_burst_anes);
    [spike_eeg_phase_itr{nn}, spike_eeg_rzpmk_itr{nn}] = spike_phase_locking( ...
        uspk, delta_sig, eeg_ts, inter_burst_anes);
    [spike_eeg_phase_emer{nn}, spike_eeg_rzpmk_emer{nn}] = spike_phase_locking( ...
        uspk, delta_sig, eeg_ts, t_emer);
    
    
    % spike-LFP phase locking within burst period
    [spike_lfp_phase_bsl{nn}, spike_lfp_rzpmk_bsl{nn}] = spike_phase_locking( ...
        uspk, mean(lfp_delta(:, unit_ch{nn}), 2), lfp_ts, t_bsl);
    [spike_lfp_phase_burst{nn}, spike_lfp_rzpmk_burst{nn}] = spike_phase_locking( ...
        uspk, mean(lfp_delta(:, unit_ch{nn}), 2), lfp_ts, delta_burst_anes);
    [spike_lfp_phase_itr{nn}, spike_lfp_rzpmk_itr{nn}] = spike_phase_locking( ...
        uspk, mean(lfp_delta(:, unit_ch{nn}), 2), lfp_ts, inter_burst_anes);
    [spike_lfp_phase_emer{nn}, spike_lfp_rzpmk_emer{nn}] = spike_phase_locking( ...
        uspk, mean(lfp_delta(:, unit_ch{nn}), 2), lfp_ts, t_emer);
end

% ratio of burst spiking occurrence in delta-burst for each unit
burst_ratio_per_unit = mean(cat(2, burst_per_delta_burst{:}) > 0, 1);

% ratio of burst spiking occurrence by unit for each delta-burst event
burst_ratio_per_delta_burst = mean(cat(2, burst_per_delta_burst{:}) > 0, 2);

% population JBSI synchrony within baseline / burst / inter-burst
unit_spike_ts = arrayfun(@(x)(spike_ts(spike_unit==x)), ...
    unique(spike_unit), 'UniformOutput', false);
sync_jbsi_bsl = get_population_jbsi(unit_spike_ts, opt.sync_window, t_bsl);
sync_jbsi_burst = get_population_jbsi(unit_spike_ts, opt.sync_window, delta_burst_anes);
sync_jbsi_itr = get_population_jbsi(unit_spike_ts, opt.sync_window, inter_burst_anes);
sync_jbsi_emer = get_population_jbsi(unit_spike_ts, opt.sync_window, t_emer);

% population CCG within baseline / burst / inter-burst
[ccg_prob_bsl, ccg_zscore_bsl, ccg_perc_bsl, ccg_t_lag] = ...
    get_population_ccg(unit_spike_ts, opt.ccg_max_lag, ...
    opt.ccg_tbin, opt.jitter_window, t_bsl);
[ccg_prob_burst, ccg_zscore_burst, ccg_perc_burst] = ...
    get_population_ccg(unit_spike_ts, opt.ccg_max_lag, ...
    opt.ccg_tbin, opt.jitter_window, delta_burst_anes);
[ccg_prob_itr, ccg_zscore_itr, ccg_perc_itr] = ...
    get_population_ccg(unit_spike_ts, opt.ccg_max_lag, ...
    opt.ccg_tbin, opt.jitter_window, inter_burst_anes);
[ccg_prob_emer, ccg_zscore_emer, ccg_perc_emer] = ...
    get_population_ccg(unit_spike_ts, opt.ccg_max_lag, ...
    opt.ccg_tbin, opt.jitter_window, t_emer);

fprintf('%gs.', toc(t));

%% save outputs
t = tic;
fprintf('\n\tsaving ... ');
outpath = fullfile(workpath, 'outputs.mat');
if ~isfile(outpath)
    save(outpath, 'eeg_ts', 'ev_mrk', 'delta_sig', 'delta_burst_peak_t_anes', ...
    'delta_burst_start', 'delta_burst_stop', '-v7.3');
end    
save(outpath, 'lfp_sig', 'lfp_ts', 'lfp_fs', 'spike_ts', 'spike_cluster', ...
    'cluster_ch', 'ch_grp', 'ex_ch', 'ex_cluster', 'opt', 'db_out', ...
    'spike_unit', 'unit_cluster', 'unit_ch', 'unit_grp', ...
    'pets_lfp_sig', 'pets_led_lfp', 'lfp_eeg_delta_dt', 'lfp_eeg_delta_dphi', ...
    'spike_eeg_phase_bsl', 'spike_eeg_rzpmk_bsl', 'spike_lfp_phase_bsl', 'spike_lfp_rzpmk_bsl', ...
    'spike_eeg_phase_burst', 'spike_eeg_rzpmk_burst', 'spike_lfp_phase_burst', 'spike_lfp_rzpmk_burst', ...
    'spike_eeg_phase_itr', 'spike_eeg_rzpmk_itr', 'spike_lfp_phase_itr', 'spike_lfp_rzpmk_itr', ...
    'spike_eeg_phase_emer', 'spike_eeg_rzpmk_emer', 'spike_lfp_phase_emer', 'spike_lfp_rzpmk_emer', ...
    'spike_delta_psth', 'spike_led_psth', 'spike_burst_vec', 'spike_burst_time', ...
    'spike_burst_idx_bsl', 'spike_burst_idx_burst', 'spike_burst_idx_itr', 'spike_burst_idx_emer', ...
    'burst_per_delta_burst', 'burst_ratio_per_unit', 'burst_ratio_per_delta_burst', ...
    'sync_jbsi_bsl', 'sync_jbsi_itr', 'sync_jbsi_burst', 'sync_jbsi_emer', 'ccg_t_lag', ...
    'ccg_prob_bsl', 'ccg_zscore_bsl', 'ccg_perc_bsl', ...
    'ccg_prob_burst', 'ccg_zscore_burst', 'ccg_perc_burst', ...
    'ccg_prob_itr', 'ccg_zscore_itr', 'ccg_perc_itr', ...
    'ccg_prob_emer', 'ccg_zscore_emer', 'ccg_perc_emer', ...
    '-append');
fprintf('%gs.', toc(t));

%% visualize results
t = tic;
fprintf('\n\tplotting lfp coherence ... ');
workpath_chs = fullfile(workpath, 'lfp_eeg_coherence');
if ~isfolder(workpath_chs); mkdir(workpath_chs); end
plot_lfp_eeg_shift(opt, workpath_chs, ...
    lfp_eeg_delta_dt, lfp_eeg_delta_dphi, ch_grp);
fprintf('%gs.', toc(t));

t = tic;
fprintf('\n\tplotting unit indexes ... ');
workpath_units = fullfile(workpath, 'unit_spike_metrics');
if ~isfolder(workpath_units); mkdir(workpath_units); end
plot_unit_spike_metrics(opt, workpath_units, spike_ts, spike_unit, ...
    unit_cluster, unit_ch, unit_grp, spike_burst_vec, delta_burst_peak_t_anes, ...
    spike_delta_psth, pets_lfp_sig, spike_burst_idx_burst, spike_burst_idx_itr, ...
    spike_eeg_phase_burst,spike_eeg_rzpmk_burst, spike_lfp_phase_burst, spike_lfp_rzpmk_burst, ...
    spike_eeg_phase_itr,spike_eeg_rzpmk_itr, spike_lfp_phase_itr, spike_lfp_rzpmk_itr, ...
    spike_eeg_phase_bsl,spike_eeg_rzpmk_bsl, spike_lfp_phase_bsl, spike_lfp_rzpmk_bsl, ...
    spike_eeg_phase_emer,spike_eeg_rzpmk_emer, spike_lfp_phase_emer, spike_lfp_rzpmk_emer);

fprintf('%gs.', toc(t));

if opt.plot_details
    t = tic;
    fprintf('\n\tplotting burst raster ... ');
    workpath_events = fullfile(workpath, 'burst_lfp_spike');
    if ~isfolder(workpath_events); mkdir(workpath_events); end
    plot_burst_lfp_spike(opt, workpath_events, delta_sig, eeg_ts, spike_ts, ...
        spike_unit, unit_ch, spike_burst_vec, delta_burst_peak_t_anes, ...
        pets_lfp_sig);
    fprintf('%gs.', toc(t));
end

t = tic;
fprintf('\n\tplotting summary ... ');
plot_led_lfp_spike(opt, workpath, delta_sig, eeg_ts, spike_ts, spike_unit, unit_ch, ...
    spike_burst_vec, t_led, pets_led_lfp, spike_led_psth, spike_burst_idx_bsl);
plot_channel_unit_summary(opt, workpath, lfp_eeg_delta_dt, lfp_eeg_delta_dphi, ...
    spike_lfp_rzpmk_bsl, spike_eeg_rzpmk_bsl, ...
    spike_lfp_rzpmk_burst, spike_eeg_rzpmk_burst, ...
    spike_lfp_rzpmk_itr, spike_eeg_rzpmk_itr, ...
    spike_lfp_rzpmk_emer, spike_eeg_rzpmk_emer, ...
    spike_burst_idx_burst, spike_burst_idx_itr, ...
    ccg_t_lag, ccg_zscore_bsl, ccg_zscore_burst, ccg_zscore_itr, ccg_zscore_emer);
fprintf('%gs.', toc(t));

end

function [spks, clusters, chs, grps] = exclude_channel_cluster( ...
    spks, clusters, chs, grps, ex_ch, ex_cluster)
xch = ismember(cat(1, grps{:, 1}), ex_ch);
xcluster = ismember(chs(:, 2), ex_ch) | ismember(chs(:, 1), ex_cluster);
xspk = ismember(clusters, chs(xcluster, 1));
spks = spks(~xspk);
clusters = clusters(~xspk);
chs = chs(~xcluster, :);
grps = grps(~xch, :);
end

function [spku, ucluster, uch, ugrp] = merge_spike_to_unit(clusters, chs, grps, mtype)
spku = zeros(size(clusters));
switch mtype
    case 'cluster'
        ucluster = num2cell(unique(chs(:, 1)));
        num_unit = length(ucluster);
        [uch, ugrp] = deal(cell(num_unit, 1));
        for nn = 1:num_unit
            uch{nn} = chs(chs(:, 1)==ucluster{nn}, 2);
            ugrp{nn} = grps{cat(1, grps{:, 1})==uch{nn}, 2};
            spku(clusters==ucluster{nn}) = nn;
        end
    case 'channel'
        uch = num2cell(unique(chs(:, 2)));
        num_unit = length(uch);
        [ucluster, ugrp] = deal(cell(num_unit, 1));
        for nn = 1:num_unit
            ugrp{nn} = grps{cat(1, grps{:, 1})==uch{nn}, 2};
            ucluster{nn} = sort(chs(chs(:, 2)==uch{nn}, 1));
            spku(ismember(clusters, ucluster{nn})) = nn;
        end
    case 'group'
        ugrp = unique(grps(:, 2));
        num_unit = length(ugrp);
        [ucluster, uch] = deal(cell(num_unit, 1));
        for nn = 1:num_unit
            uch{nn} = cat(1, grps{strcmp(grps(:, 2), ugrp{nn}), 1});
            ucluster{nn} = sort(chs(ismember(chs(:, 2), uch{nn}), 1));
            spku(ismember(clusters, ucluster{nn})) = nn;
        end
end
end

function [dt, dphi] = get_peak_shift_t_phi(shift_sig, shift_ts, ref_sig, ref_ts, tev)
    [~, pkt] = findpeaks(shift_sig, shift_ts);
    [~, imax] = min(abs(tev - pkt'), [], 2);
    dt = pkt(imax)-tev;
    dphi = interp1(ref_ts, angle(hilbert(ref_sig)), pkt(imax));  
end

function [pe_hist, hist_cnt] = get_peri_event_hist(spk, ev, twin)
nev = length(ev);
nt = length(twin);
pe_hist = zeros(nev, nt-1);
for kk = 1:nev
    if ~isnan(ev(kk))
        pe_hist(kk, :) = histcounts(spk, twin + ev(kk));
    end
end
pe_hist = pe_hist/median(diff(twin));
hist_cnt = (twin(1:end-1) + twin(2:end))/2;
end

function idx = compute_spiking_index(spk, isburst, ibeg, iend, trg)
idx = table;
% firing rate
idx.spike_fr = sum(mask_time_range(spk, trg))/sum(diff(trg, 1, 2));
% burst spikes per second
idx.burst_fr = sum(mask_time_range(spk(isburst), trg))/sum(diff(trg, 1, 2));
% number of bursts per second
idx.burst_rate = sum(mask_time_range(spk(ibeg), trg))/sum(diff(trg, 1, 2));
% number of spikes per burst
idx.spike_per_burst = mean(diff(mask_time_range(spk(ibeg), ...
    trg, [ibeg, iend]), 1, 2)) + 1;
% burst spiking ratio
idx.spike_burst_perc = mean(mask_time_range(spk, trg, isburst))*100;
% inter-spike interval
idx.spike_isi = median(diff(mask_time_range(spk, trg, spk)));
% firing rate variation
idx.spike_fr_sd = std(histcounts(spk, (ceil(min(trg)):1:floor(max(trg)))));
end

function str = join_array_string(arr)
str = strjoin(arrayfun(@(x)(num2str(x)), arr, 'UniformOutput', false), ',');
end

function plot_lfp_eeg_shift(opt, workpath, dt, dphi, grps)
time_edges = (-0.2:0.01:0.2);
phase_edges = (-180:10:180);

for jj = 1:size(grps, 1)
    ch = grps{jj, 1};
    str_ch = num2str(ch, '%.2d');
    str_grp = grps{jj, 2};
    figure; set(gcf, 'Position', [100, 100, 790, 800], 'Visible', false, ...
        'PaperOrientation', 'portrait');
    
    subplot('Position', [0.35, 0.54, 0.36, 0.35]); hold on;
    histogram(dt{ch}, time_edges, 'Normalization', 'probability');
    set(gca, 'TickDir', 'out', 'XLim', minmax(time_edges));
    xlabel('Time to EEG delta peak (s)');
    ylabel('Probability');
    title(sprintf('Phase shift %.2fs', median(dt{ch})));
    
    subplot('Position', [0.35, 0.07, 0.36, 0.35]); hold on;
    histogram(rad2deg(dphi{ch}), phase_edges, 'Normalization', 'probability');
    set(gca, 'TickDir', 'out', 'XLim', minmax(phase_edges), 'XTick', ...
        (-180:90:180), 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('Phase to EEG delta peak');
    ylabel('Probability');
    title(sprintf('Phase shift: %.2f (%.0f\\circ)', median(dphi{ch}), ...
        rad2deg(median(dphi{ch}))));
    
    sgtitle(sprintf('%s\nChannel #%s, Group %s', workpath, ...
        str_ch, str_grp), 'Interpreter', 'none');
    save_multi_formats(gcf, fullfile(workpath, ...
        strcat('lfp_eeg_shift_ch#', str_ch)), opt.fig_save_format);
    close;
end
end

function plot_unit_spike_metrics(opt, workpath, spks, spku, uclusters, uchs, ugrps, ...
    is_burst, pkt, psths, pets_lfp_sig, spk_idx_burst, spk_idx_itr, ...
    eeg_phi_burst, eeg_rzpmk_burst, lfp_phi_burst, lfp_rzpmk_burst, ...
    eeg_phi_itr, eeg_rzpmk_itr, lfp_phi_itr, lfp_rzpmk_itr, ...
    eeg_phi_bsl, eeg_rzpmk_bsl, lfp_phi_bsl, lfp_rzpmk_bsl, ...    
    eeg_phi_emer, eeg_rzpmk_emer, lfp_phi_emer, lfp_rzpmk_emer)
    
ds_bin = 5;
phase_edges = (-180:10:180);
pk_rg = [-1, 1];

pk_win = opt.peri_event_window;
pk_win_low = pk_win(1:ds_bin:end);
pk_win_low_cnt = (pk_win_low(1:end-1)+pk_win_low(2:end))/2;
nu = length(unique(spku));
nev = length(pkt);
for ii = 1:nu
    str_cluster = join_array_string(uclusters{ii}(1:min(length(uclusters{ii}), 15)));
    str_ch = join_array_string(uchs{ii});
    str_grp = ugrps{ii};
    uspk = spks(spku==ii);    
    figure; set(gcf, 'Position', [100, 100, 400, 800], 'Visible', false, ...
        'PaperOrientation', 'portrait');
    
    subplot('Position', [0.10, 0.65, 0.80, 0.15]); hold on;
    temp = mean(cat(3, pets_lfp_sig{uchs{ii}}), 3);
    plot(pk_win, temp', 'Color', [.7, .7, .7], 'LineWidth', 0.3);
    plot(pk_win, mean(temp, 1), '-r', 'LineWidth', 1.3);
    plot([0, 0], get(gca, 'YLim'), '--k');
    axis off;
    set(gca, 'XLim', pk_rg);
    str_burst = sprintf('Delta-burst FR:%6.1f Hz,%4.1f%% burst', ...
        spk_idx_burst{ii}.spike_fr, spk_idx_burst{ii}.spike_burst_perc);
    str_inter = sprintf('Inter-burst FR:%6.1f Hz,%4.1f%% burst', ...
        spk_idx_itr{ii}.spike_fr, spk_idx_itr{ii}.spike_burst_perc);
    title(sprintf('%s\n%s', str_burst, str_inter));

    subplot('Position', [0.10, 0.35, 0.80, 0.30]); hold on;
    for tt = 1:nev
        in0 = uspk>pk_win(1)+pkt(tt) & uspk<pk_win(end)+pkt(tt);
        inb = in0 & is_burst{ii};
        plot(uspk(in0)-pkt(tt), uspk(in0)*0+tt, '.', ...
            'Color', [.7, .7, .7], 'LineWidth', 0.3);
        plot(uspk(inb)-pkt(tt), uspk(inb)*0+tt, '.b');
    end
    plot([0, 0], [0, nev+1], '--k');
    set(gca, 'YLim', [0, nev+1], 'XLim', pk_rg, 'Visible', 'off');
    
    subplot('Position', [0.10, 0.10, 0.80, 0.20]); hold on;
    temp = mean(squeeze(mean(reshape(psths{ii}, [], ds_bin, ...
        size(psths{ii}, 2)/ds_bin), 2)), 1);
    bar(pk_win_low_cnt, temp);
    plot([0, 0], get(gca, 'YLim'), '--k');
    set(gca, 'XLim', pk_rg);
    xlabel('Time to peak of filtered EEG signal (s)');
    ylabel('Averaged firing rate (Hz)');   

    sgtitle(sprintf('%s\nCluster #%s, Channel #%s, Group %s', workpath, ...
        str_cluster, str_ch, str_grp), 'Interpreter', 'none');
    save_multi_formats(gcf, fullfile(workpath, ...
        strcat('unit_spike_metrics_cluster#-1', str_cluster)), opt.fig_save_format);
    close;


    figure; set(gcf, 'Position', [100, 100, 1000, 500], 'Visible', false, ...
        'PaperOrientation', 'portrait');
    
    subplot(2, 4, 1); hold on;
    histogram(eeg_phi_burst{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'm', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(eeg_rzpmk_burst{ii}(4), (1+cos(eeg_rzpmk_burst{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('EEG phase');
    ylabel('Probability of spiks');
    title(sprintf('(Burst) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        eeg_rzpmk_burst{ii}(2), eeg_rzpmk_burst{ii}(3), eeg_rzpmk_burst{ii}(4), ...
        rad2deg(eeg_rzpmk_burst{ii}(4)), eeg_rzpmk_burst{ii}(5)));
    
    subplot(2, 4, 5); hold on;
    histogram(lfp_phi_burst{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'c', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(lfp_rzpmk_burst{ii}(4), (1+cos(lfp_rzpmk_burst{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('LFP phase');
    ylabel('Probability of spiks');
    title(sprintf('(Burst) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        lfp_rzpmk_burst{ii}(2), lfp_rzpmk_burst{ii}(3), lfp_rzpmk_burst{ii}(4), ...
        rad2deg(lfp_rzpmk_burst{ii}(4)), lfp_rzpmk_burst{ii}(5)));

    subplot(2, 4, 2); hold on;
    histogram(eeg_phi_itr{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'm', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(eeg_rzpmk_itr{ii}(4), (1+cos(eeg_rzpmk_itr{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('EEG phase');
    ylabel('Probability of spiks');
    title(sprintf('(itr) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        eeg_rzpmk_itr{ii}(2), eeg_rzpmk_itr{ii}(3), eeg_rzpmk_itr{ii}(4), ...
        rad2deg(eeg_rzpmk_itr{ii}(4)), eeg_rzpmk_itr{ii}(5)));
    
    subplot(2, 4, 6); hold on;
    histogram(lfp_phi_itr{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'c', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(lfp_rzpmk_itr{ii}(4), (1+cos(lfp_rzpmk_itr{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('LFP phase');
    ylabel('Probability of spiks');
    title(sprintf('(itr) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        lfp_rzpmk_itr{ii}(2), lfp_rzpmk_itr{ii}(3), lfp_rzpmk_itr{ii}(4), ...
        rad2deg(lfp_rzpmk_itr{ii}(4)), lfp_rzpmk_itr{ii}(5)));

    subplot(2, 4, 3); hold on;
    histogram(eeg_phi_bsl{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'm', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(eeg_rzpmk_bsl{ii}(4), (1+cos(eeg_rzpmk_bsl{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('EEG phase');
    ylabel('Probability of spiks');
    title(sprintf('(bsl) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        eeg_rzpmk_bsl{ii}(2), eeg_rzpmk_bsl{ii}(3), eeg_rzpmk_bsl{ii}(4), ...
        rad2deg(eeg_rzpmk_bsl{ii}(4)), eeg_rzpmk_bsl{ii}(5)));
    
    subplot(2, 4, 7); hold on;
    histogram(lfp_phi_bsl{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'c', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(lfp_rzpmk_bsl{ii}(4), (1+cos(lfp_rzpmk_bsl{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('LFP phase');
    ylabel('Probability of spiks');
    title(sprintf('(bsl) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        lfp_rzpmk_bsl{ii}(2), lfp_rzpmk_bsl{ii}(3), lfp_rzpmk_bsl{ii}(4), ...
        rad2deg(lfp_rzpmk_bsl{ii}(4)), lfp_rzpmk_bsl{ii}(5)));

    subplot(2, 4, 4); hold on;
    histogram(eeg_phi_emer{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'm', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(eeg_rzpmk_emer{ii}(4), (1+cos(eeg_rzpmk_emer{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('EEG phase');
    ylabel('Probability of spiks');
    title(sprintf('(emer) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        eeg_rzpmk_emer{ii}(2), eeg_rzpmk_emer{ii}(3), eeg_rzpmk_emer{ii}(4), ...
        rad2deg(eeg_rzpmk_emer{ii}(4)), eeg_rzpmk_emer{ii}(5)));
    
     subplot(2, 4, 8); hold on;
    histogram(lfp_phi_emer{ii}, deg2rad(phase_edges), 'Normalization', 'probability', ...
        'EdgeColor', 'w', 'FaceColor', 'c', 'FaceAlpha', 0.5);
    ymax = max(get(gca, 'YLim'))*1.05;
    plot((-pi:0.1:pi), (1+cos((-pi:0.1:pi)))*ymax/2, 'k', 'LineWidth', 1.5);
    plot([0, 0], [0, ymax], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
    plot(lfp_rzpmk_emer{ii}(4), (1+cos(lfp_rzpmk_emer{ii}(4)))*ymax/2, 'ok', ...
        'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
    set(gca, 'TickDir', 'out', 'XLim', [-pi, pi], 'XTick', ...
        (-1:0.5:1)*pi, 'XTickLabel', {'-\pi', '-1/2\pi', '0', '1/2\pi', '\pi'});
    xlabel('LFP phase');
    ylabel('Probability of spiks');
    title(sprintf('(emer) Z=%.2f, \\itp\\rm=%.3f, \\mu=%.2f (%.0f\\circ), \\kappa=%.3f', ...
        lfp_rzpmk_emer{ii}(2), lfp_rzpmk_emer{ii}(3), lfp_rzpmk_emer{ii}(4), ...
        rad2deg(lfp_rzpmk_emer{ii}(4)), lfp_rzpmk_emer{ii}(5)));
        
    sgtitle(sprintf('%s\nCluster #%s, Channel #%s, Group %s', workpath, ...
        str_cluster, str_ch, str_grp), 'Interpreter', 'none');
    save_multi_formats(gcf, fullfile(workpath, ...
        strcat('unit_spike_metrics_cluster#-2', str_cluster)), opt.fig_save_format);
    close;
end
end

function plot_burst_lfp_spike(opt, workpath, sig, ts, spks, spku, uchs, ...
    is_burst, pkt, pets_lfp_sig)
pk_rg = [-1, 2];
bar_yv = [-0.3; 0.3];

pk_win = opt.peri_event_window;
nu = length(unique(spku));
nev = length(pkt);
traces = interp1(ts, sig, pkt+pk_win);
lfps = cat(3, pets_lfp_sig{unique(cat(1, uchs{:}))}); 
for tt = 1:nev
    str_ev = num2str(tt, '%.3d');
    figure; set(gcf, 'Position', [100, 100, 400, 600], 'Visible', false, ...
        'PaperOrientation', 'portrait');
    
    subplot('Position', [0.15, 0.65, 0.70, 0.2]); hold on;
    temp = permute(lfps(tt, :, :), [2, 3, 1]);
    plot(pk_win, temp, 'Color', [.7, .7, .7], 'LineWidth', 0.3);
    plot(pk_win, traces(tt, :), '-r', 'LineWidth', 1.3);
    plot([0, 0], get(gca, 'YLim'), '--k');
    axis off;
    set(gca, 'XLim', pk_rg);
    title(sprintf('%s\nDelta-burst #%s, Peak time: %.3fs', workpath, ...
        str_ev, pkt(tt)), 'Interpreter', 'none', 'FontSize', 10);

    subplot('Position', [0.15, 0.15, 0.70, 0.50]); hold on;
    plot([0, 0], [0, nu+1], '--k');
    for ii = 1:nu
        uspk = spks(spku==ii); 
        in0 = uspk>pk_win(1)+pkt(tt) & uspk<pk_win(end)+pkt(tt);
        inb = in0 & is_burst{ii};
        plot(repmat(uspk(in0)'-pkt(tt), 2, 1), uspk(in0)'*0+ii+bar_yv, '-', ...
            'Color', [.5, .5, .5], 'LineWidth', 0.3);
        plot(repmat(uspk(inb)'-pkt(tt), 2, 1), uspk(inb)'*0+ii+bar_yv, '-b');
    end
    set(gca, 'YLim', [0, nu+1], 'XLim', pk_rg, 'TickDir', 'out');
    set(get(gca, 'YAxis'), 'Visible', 'off');
    xlabel('Time to peak of filtered EEG signal (s)');

    save_multi_formats(gcf, fullfile(workpath, ...
        strcat('raster_plot_delta_burst_#', str_ev)), ...
        opt.fig_save_format);
    close;
end
end

function plot_led_lfp_spike(opt, workpath, sig, ts, spks, spku, uchs, ...
    is_burst, evt, pets_lfp_sig, psths, spk_idx_bsl)
sm_sigma = 5;
bar_yv = [-0.3; 0.3];

led_win = opt.peri_led_window;
led_win_cnt = (led_win(1:end-1)+led_win(2:end))/2;
nu = numel(unique(spku));
nev = numel(evt);
sm_kernel = normpdf((-3*sm_sigma:1:3*sm_sigma), 0, sm_sigma);
psths = cellfun(@(x)(conv2(x, sm_kernel, 'same')), psths, 'UniformOutput', false);
psths = permute(cat(3, psths{:}), [2, 3, 1])./cellfun(@(x)(x.spike_fr), spk_idx_bsl)';
lfps = permute(cat(3, pets_lfp_sig{unique(cat(1, uchs{:}))}), [2, 3, 1]);

for tt = 1:nev
    str_ev = num2str(tt, '%.3d');
    figure; set(gcf, 'Position', [100, 100, 600, 800], 'Visible', false, ...
        'PaperOrientation', 'portrait');

    subplot('Position', [0.15, 0.65, 0.70, 0.20]); hold on;
    traces = interp1(ts, sig, evt(tt)+led_win);
    plot(led_win, lfps(:, :, tt), 'Color', [.7, .7, .7], 'LineWidth', 0.3);
    plot(led_win, traces, '-r', 'LineWidth', 0.8);
    xline(0, '--k');
    axis off;
    set(gca, 'XLim', minmax(led_win));

    subplot('Position', [0.15, 0.32, 0.70, 0.33]); hold on;
    plot([0, 0], [0, nu+1], '--k');
    for ii = 1:nu
        uspk = spks(spku==ii);
        in0 = uspk>led_win(1)+evt(tt) & uspk<led_win(end)+evt(tt);
        inb = in0 & is_burst{ii};
        plot(repmat(uspk(in0)'-evt(tt), 2, 1), uspk(in0)'*0+ii+bar_yv, '-', ...
            'Color', [.5, .5, .5], 'LineWidth', 0.3);
        plot(repmat(uspk(inb)'-evt(tt), 2, 1), uspk(inb)'*0+ii+bar_yv, '-b');
    end
    set(gca, 'YLim', [0, nu+1], 'XLim', minmax(led_win), 'Visible', 'off');

    subplot('Position', [0.15, 0.10, 0.70, 0.20]); hold on;
    xline(0, '--k');
    temp = psths(:, :, tt);
    shadedErrorBar(led_win_cnt, mean(temp, 2), std(temp, 1, 2)/sqrt(nu), ...
        'lineProps', {'-k', 'LineWidth', 1});
    set(gca, 'TickDir', 'out', 'XLim', minmax(led_win), 'YLim', ...
        min(get(gca, 'YLim'), [0, inf]));
    xlabel('Time to LED onset (s)');
    ylabel('Averaged normalized firing rate (a.u.)');

    sgtitle(sprintf('%s\nLED pulse #%s, Time: %.3fs', workpath, ...
        str_ev, evt(tt)), 'Interpreter', 'none');
    save_multi_formats(gcf, fullfile(workpath, ...
        strcat('raster_plot_led_pulse_#', str_ev)), ...
        opt.fig_save_format);
    close;
end
end

function plot_channel_unit_summary(opt, workpath, dts, dphis, ...
    lfp_rzpmk_bsl, eeg_rzpmk_bsl, lfp_rzpmk_burst, eeg_rzpmk_burst, ...
    lfp_rzpmk_itr, eeg_rzpmk_itr, lfp_rzpmk_emer, eeg_rzpmk_emer, ...
    spk_idx_burst, spk_idx_itr, lag, ...
    cc_bsl, cc_burst, cc_itr, cc_emer)

figure; set(gcf, 'Position', [100, 100, 400, 400], 'Visible', false, ...
    'PaperOrientation', 'portrait');
subplot('Position', [0.20, 0.12, 0.70, 0.70]);
plot_summary_lfp_eeg_shift( ...
    cellfun(@(x)(mean(x)), dts), cellfun(@(x)(mean(x)), dphis));
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'lfp_eeg_shift'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, lfp_rzpmk_bsl{:}), 'spike-LFP phase locking (Bsaeline)', 'c');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, eeg_rzpmk_bsl{:}), 'spike-EEG phase locking (Baseline)', 'm');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'spike_phase_locking (Bsaeline)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, lfp_rzpmk_burst{:}), 'spike-LFP phase locking (Delta burst)', 'c');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, eeg_rzpmk_burst{:}), 'spike-EEG phase locking (Delta burst)', 'm');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'spike_phase_locking (Delta burst)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, lfp_rzpmk_itr{:}), 'spike-LFP phase locking (Inter burst)', 'c');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, eeg_rzpmk_itr{:}), 'spike-EEG phase locking (Inter burst)', 'm');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'spike_phase_locking (Inter burst)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, lfp_rzpmk_emer{:}), 'spike-LFP phase locking (Emergence)', 'c');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_spike_phase_locking( ...
    cat(1, eeg_rzpmk_emer{:}), 'spike-EEG phase locking (Emergence)', 'm');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'spike_phase_locking (Emergence)'), ...
    opt.fig_save_format);
close;



figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, lfp_rzpmk_bsl{:}), 'spike-LFP phase locking (Baseline)');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, eeg_rzpmk_bsl{:}), 'spike-EEG phase locking (Baseline)');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'phase_locking_dist (Baseline)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, lfp_rzpmk_burst{:}), 'spike-LFP phase locking (Delta burst)');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, eeg_rzpmk_burst{:}), 'spike-EEG phase locking (Delta burst)');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'phase_locking_dist (Delta burst)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, lfp_rzpmk_itr{:}), 'spike-LFP phase locking (Inter burst)');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, eeg_rzpmk_itr{:}), 'spike-EEG phase locking (Inter burst)');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'phase_locking_dist (Inter burst)'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 800, 400], 'Visible', false, ...
    'PaperOrientation', 'landscape');
subplot('Position', [0.12, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, lfp_rzpmk_emer{:}), 'spike-LFP phase locking (Emergence)');
subplot('Position', [0.60, 0.10, 0.35, 0.65]);
plot_summary_phase_locking_dist( ...
    cat(1, eeg_rzpmk_emer{:}), 'spike-EEG phase locking (Emergence)');
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'phase_locking_dist (Emergence)'), ...
    opt.fig_save_format);
close;


figure; set(gcf, 'Position', [100, 100, 790, 790], 'Visible', false, ...
    'PaperOrientation', 'portrait');
plot_summary_firing_burst(cat(1, spk_idx_burst{:}), cat(1, spk_idx_itr{:}));
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'firing_rate_burst_percentage'), ...
    opt.fig_save_format);
close;

figure; set(gcf, 'Position', [100, 100, 400, 400], 'Visible', false, ...
    'PaperOrientation', 'portrait');
subplot('Position', [0.20, 0.12, 0.70, 0.70]);
plot_summary_population_ccg(lag, cc_bsl, cc_burst, cc_itr, cc_emer);
sgtitle(workpath, 'Interpreter', 'none');
save_multi_formats(gcf, fullfile(workpath, 'population_cross_correlogram'), ...
    opt.fig_save_format);
close;

end
