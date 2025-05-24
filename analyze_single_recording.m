function ops = analyze_single_recording(ds, opt, wp)

ch_grp = ds.ChGrp;
eeg_ch = cat(1, ch_grp{strcmpi(ch_grp(:, 2), 'eeg'), 1});
eeg_fs = ds.LFPFs;
eeg_ts = ds.LFPTs;
ev_mrk = ds.EvMrk;

if isfield(ds, 'Ca2Path') && ~isempty(ds.Ca2Path)
    do_ca2 = true; 
    ca2_sig = ds.Ca2Sig;
    ca2_fs = ds.Ca2Fs;
    ca2_ts = ds.Ca2Ts;
else
    do_ca2 = false; 
end

if isfield(ds, 'SpkTs') && ~isempty(ds.SpkTs)
    do_spk = true; 
    lfp_sig = ds.LFPSig;
    lfp_fs = ds.LFPFs;
    lfp_ts = ds.LFPTs;
    spike_ts = ds.SpkTs;
    spike_cluster = ds.SpkCluster;
    cluster_ch = ds.ClusterCh;
    ex_ch = ds.ExCh;
    ex_cluster = ds.ExCluster;
else
    do_spk = false; 
end

num_eeg = length(eeg_ch);
ops = cell(num_eeg, 1);
for ii = 1:num_eeg
    ti = tic;
    eeg_name = sprintf('EEGCh#%.2d', eeg_ch(ii));
    op = fullfile(wp, eeg_name);
    if ~isfolder(op); mkdir(op); end
    fprintf('\n<%d> %s', ii, op);
    
    fprintf('\ndelta burst ... ');
    eeg_sig = ds.LFPSig(:, eeg_ch(ii));
    try
        db_out = analyze_delta_burst(eeg_sig, eeg_ts, eeg_fs, ev_mrk, opt, op);
    catch except
        warning('\n%s', getReport(except));
    end
        
    if do_spk
        fprintf('\nLFP & spike ... ');
        try
            spk_out = analyze_related_lfp_spike(lfp_sig, lfp_ts, lfp_fs, ...
                spike_ts, spike_cluster, cluster_ch, ch_grp, ...
                ex_ch, ex_cluster, opt, db_out, op);
        catch except
            warning('\n%s', getReport(except));
        end
    end

    if do_ca2
        fprintf('\ncalcium signal ... ');
        try
            ca2_out = analyze_related_ca2(ca2_sig, ca2_ts, ca2_fs, opt, db_out, op);
        catch except
            warning('\n%s', getReport(except));
        end
    end
    
    ops{ii} = op;
    fprintf('\ndone in %gs.\n', toc(ti));
end
end

