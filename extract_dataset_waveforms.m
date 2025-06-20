%%
% individual dataset generated by generate_dataset_lfp_spike.m
dataset = ds(4);
% absolute path to save figure of waveforms
workpath = '';
% time for extraction refer to recording onset, unit: sec
time_range = ds(4).EvMrk.t_led + [-300, 300];
% waveforms of clusters to plot, DEFAULT: all clusters
cluster_list = [];
% number of timestamps before each spike peak to extract
n_before = 16;
% number of timestamps after each spike peak to extract
n_after = 16;
% gap between two channel waveforms to plot along y-axis
yscale_uv = 600;

subextract_klusta_waveforms(dataset, workpath, time_range, cluster_list, ...
    n_before, n_after, yscale_uv);


