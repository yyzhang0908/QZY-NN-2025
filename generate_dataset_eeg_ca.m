%% load data files to dataset (EEG-Ca2)
% 'new': create a new data structure; 'append': add a recording to the
% existing data structure
read_mode = 'append';
% number of data file, positive integer
file_id = 1;
% absolute path of eeg file
eeg_path = '.mat';
% absolute path of calcium signal file
ca2_path = '.mat';

% set events timestamps/markers
ev_mrk = struct;
% time of injection, unit: sec
ev_mrk.t_med = 0; 
% required periods refer to recording onset, unit: sec
ev_mrk.t_bsl = [10, 150];        % baseline 
ev_mrk.t_anes = [1400, 1600];     % anesthesia 
ev_mrk.t_emer = [2400, 2700];     % emergence  
% time of optogenetical pulse, default: nan, unit: sec
ev_mrk.t_led = nan;
% ds(1).EvMrk = ev_mrk;

% eletrode information
grp_ch = {1};
grp_label = {'EEG'};

% load eeg data and get the required information
load(eeg_path);
% get eeg signal, unit: uV
eeg_sig = data(:) * 1e6;
% eeg sampling rate
eeg_fs = 1000;
% eeg timestamps refer to injection, unit: sec
eeg_ts = (0:1:length(eeg_sig)-1)' / eeg_fs;

% get group label of each channel
ch_grp = arrayfun(@(x, y)(repmat(x, length(y{1}), 1)), ...
    grp_label, grp_ch, 'UniformOutput', false);
ch_grp = cat(2, num2cell(cat(2, grp_ch{:})'), cat(1, ch_grp{:}));

if ~isempty(ca2_path)
    % load ca2+ data and get the required information
    load(ca2_path);
    % get preprocessed calcium signal, unit: dF/F
%     cal1 = cal1(:)/median(cal1)*100;  % cal470
%     cal0 = cal0(:)/median(cal0)*100;  % cal405
%     b = [cal0, cal0*0+1]\cal1;
%     psth1 = cal1-(cal0*b(1)+b(2));    % preprocess

    ca2_sig = data(:)/100;    
    % ca2_sig(1) = [];
%     ca2+ signal sampling rate
    ca2_fs = 1/median(diff(times));
%     ca2+ timestamps refer to injection, unit: sec
    ca2_ts = times(:);
    % ca2_ts(1) = [];
    
%     ca2_fs = 100;
%     ca2_ts = (0:1/ca2_fs:times(end))';
%     ca2_sig = interp1(times(:), data(:), ca2_ts, 'spline');
else
    [ca2_sig, ca2_fs, ca2_ts] = deal([]);
end

% add data files to structure
if strcmp(read_mode, 'new') || ~exist('ds', 'var'); ds = struct; end
ds(file_id).LFPPath = eeg_path;
ds(file_id).LFPSig = eeg_sig;
ds(file_id).LFPFs = eeg_fs;
ds(file_id).LFPTs = eeg_ts;
ds(file_id).ChGrp = ch_grp;
ds(file_id).Ca2Path = ca2_path;
ds(file_id).Ca2Sig = ca2_sig;
ds(file_id).Ca2Fs = ca2_fs;
ds(file_id).Ca2Ts = ca2_ts;
ds(file_id).EvMrk = ev_mrk;

%% save data structure
dataset_path = '\ds.mat';
save(dataset_path, 'ds', '-v7.3');

%% load data structure
dataset_path = '\ds.mat';
load(dataset_path, 'ds');
