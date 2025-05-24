%% set parameters
% define paths
script_path = '';      % folder where scripts stored
working_path = '';      % folder where outputs saved

% set configs
opt = struct;
% format of figure to be saved as, subset of {'pdf', 'jpg', 'tiff', 'fig'}
opt.fig_save_format = {'pdf', 'jpg', 'fig'};       
opt.fft_window = 2;     % window length (smoothness) used to calculate spectrogram, unit: sec
opt.fft_step = 0.05;     % step (resolution) for spectrogram, unit: sec
opt.apply_win_len = 0.5;    % apply DFT window for short events estimation, unit: sec
opt.wavelet_scale = 10;     % spectrum resolution for wavelet transform of ca2+ signal
% define frequency bands by paired labels and ranges
opt.band_freq = {[0.3, 1.5], [1.5, 4], [4, 8], [8, 13], [13, 25], [30, 100]};
opt.band_label = {'slow', 'delta','theta', 'alpha', 'beta', 'gamma'};
% define frequency bands for ca2+ signal by paired labels and ranges
opt.ca2_band_freq = {[0.1, 1], [1, 2], [2, 3], [3, 4], [1.5, 4], [4, 5]};
opt.ca2_band_label = {'slow', 'delta_1', 'delta_2', 'delta_3', 'delta', 'high'};
% bin width for calculate bandpower timecourse, unit: sec
opt.band_power_bin = 60;

% parameters to get delta burst event
% time window to smooth delta amplitude, unit: sec
opt.delta_smooth_window = 0;
% detect method: {'anes_2gm', 'bsl_rms'}
% "anes_2gm" fits delta amplitude peak under anesthesia with
% 2-gaussian model, final threshold defined as 2-fold mean value of the
% lower one. "bsl_rms" define final threshold as 3-fold root mean squared
% baseline delta amplitude.
opt.delta_burst_method = 'anes_2gm';   
% align method: {'first', 'highest', 'first-highest', 'first-skip'}
% - "first" align delta burst to their first local maximum. 
% - "highest" align to the maximum in a delta burst period.
% - "first-highest" align delta burst to their first threshold-passing
% local maximum, if not exist, use the maximum.
% - "first-skip" align delta burst to their first threshold-passing
% local maximum, if not exist, fill NaN.
opt.delta_burst_align = 'first-skip';
% bin width for calculate delta burst (duration/interval) timecourse, unit: sec
opt.delta_burst_bin = 60;    

% max lag for cross-correlation between delta and other bands
opt.delta_ccg_lag = 4;               
% time window for ca2+ singal/spike histogram/other related sources around
% delta peaks, unit: sec
opt.peri_event_window = (-3:0.01:3);    
% time window for optogenetic stimulus related electrophysiological signal, unit: sec
opt.peri_led_window = (-1:0.02:1);

% merge spikes belonging to the same {'cluster', 'channel', 'group'} 
opt.spike_merge_by = 'cluster';
% hyperpolarization period before burst beginning, unit: ms
opt.hyperpolarization_ms = 200;
% the maximum initial interval among burst spikes, unit: ms
opt.beginning_isi_ms = 20;
% the maximum final interval among burst spikes, unit: ms
opt.end_isi_ms = 60;
% minimum number of spikes within a burst
opt.minimum_sequence = 2;

% coincidence window for population jitter-based synchrony, unit: msec
opt.sync_window = 50;
% time lag range for cross-correlogram, unit: msec
opt.ccg_max_lag = 500;
% gap between each time lag, unit: msec
opt.ccg_tbin = 10;
% jitter spike time within a narrow window for permutation, unit: msec
opt.jitter_window = 30;

% export figures of burst raster, which takes a long time (~10 seconds per burst)
opt.plot_details = false;

%% main analysis
% add matlab path
addpath(genpath(script_path));

% add working folder
if ~isfolder(working_path); mkdir(working_path); end

% start parallel for acceleration
gcp;

% do single recording analysis of delta burst
num_ds = length(ds);
ops = cell(num_ds, 1);
for nn = 1:6
    %      
    tn = tic;
    fprintf('%.2d|%.2d %s ... ', nn, num_ds, ds(nn).LFPPath);
    [~, file_name] = fileparts(ds(nn).LFPPath);
    wp = fullfile(working_path, file_name);
    op = analyze_single_recording(ds(nn), opt, wp);
    ops{nn} = op;
    fprintf('totally %gs:\n', toc(tn));
    disp(op);
end
ops_all = cat(1, ops{:});
if iscell(ops_all{1})
    ops_all = cellfun(@(x)(x{1}), ops_all, 'UniformOutput', false); 
end

% finish, and kill the parallel pool
delete(gcp);

%% main statistics
% add results of single recordings manually
ops = ops_all(1:8);
% ops = { ...       
% 
%     };




% define task name for the list of results (output paths), default: 'stats'
task_name = 'stats';

% add matlab path
addpath(genpath(script_path));

% do statistics for the whole dataset
tn = tic;
fprintf('statistics ... ')
wp = fullfile(working_path, task_name);
[op, stats] = summarize_multi_outputs(ops, wp);
fprintf('%gs.\n', toc(tn));
disp(wp);

%%
% open.fig file
fig = open('');
set(fig, 'Visible', 'on');        
