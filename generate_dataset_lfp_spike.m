%% load data files to dataset (LFP-Spike)
% 'new': create a new data structure; 'append': add a recording to the
% existing data structure
read_mode = 'append';
% number of data file, positive integer
file_id = 1;
% absolute path of lfp (.ns2) file
lfp_path = '.ns2';
% absolute path of sorted (.kwik) file
spike_path = '.kwik';

% set events timestamps/markers
ev_mrk = struct;
% time of injection, unit: sec
ev_mrk.t_med = 0; 
% required periods refer to recording onset, unit: sec
ev_mrk.t_bsl = [100, 400];        % baseline 
ev_mrk.t_anes = [1200, 1800];     % anesthesia 
ev_mrk.t_emer = [35, 40]*60;     % emergence
% time of optogenetical pulse, default: nan, unit: sec
ev_mrk.t_led = 1650.323;

% electrode information             
grp_ch = {[1, 3, 5, 7], [2, 4, 6, 8], [9, 11, 13, 15], [10, 12], [14, 16]};
grp_label = { 'A', 'B', 'C', 'EMG', 'EEG'};
% electrodes (channels and clusters belonging to) to be excluded from recording statistics
ex_ch = [10, 12, 14, 16];
% clusters to be excluded from recording statistics
% ex_cluster = [1,11,13,18,21,22,23,25,26,30,31,33,34,36,37,41,45,49,50,51,55,58,62,65,67,69];
ex_cluster = [];

% load lfp data and get signal, sampling rate, timestamps
da = openNSx('read', lfp_path);
lfp_sig = double(da.Data') * da.ElectrodesInfo(1).Resolution;   % data matrix shape: timepoints x channels
lfp_fs = da.MetaTags.SamplingFreq;
lfp_ts = (0:1:size(lfp_sig, 1)-1)' / lfp_fs;

% load spike data and get sampling rate, timestamps, cluster ids, 1-index
% channel-cluster mapping
[spike_fs, spike_ts, spike_cluster, cluster_ch] = ...
    map_sorter_results(spike_path);

% get group label of each channel
ch_grp = arrayfun(@(x, y)(repmat(x, length(y{1}), 1)), ...
    grp_label, grp_ch, 'UniformOutput', false);
ch_grp = cat(2, num2cell(cat(2, grp_ch{:})'), cat(1, ch_grp{:}));

% add data files to structure
if strcmp(read_mode, 'new') || ~exist('ds', 'var'); ds = struct; end
ds(file_id).LFPPath = lfp_path;
ds(file_id).LFPSig = lfp_sig;
ds(file_id).LFPFs = lfp_fs;
ds(file_id).LFPTs = lfp_ts;
ds(file_id).ChGrp = ch_grp;
ds(file_id).SpkFs = spike_fs;
ds(file_id).SpkTs = spike_ts;
ds(file_id).SpkCluster = spike_cluster;
ds(file_id).ClusterCh = cluster_ch;
ds(file_id).ExCh = ex_ch;
ds(file_id).ExCluster = ex_cluster;
ds(file_id).EvMrk = ev_mrk;

%% save data structure
dataset_path = '\ds.mat';
save(dataset_path, 'ds', '-v7.3');

%% load data structure
dataset_path = '\test.mat';
load(dataset_path, 'ds');

%% convert sorting result to .mat
% absolute path of sorted (.kwik) file
kwikfile = '.kwik';
lines = convert_sorting_results(kwikfile);
% copy the last command line on window, paste to (python-installed)
% powershell and run.
