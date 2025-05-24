function [spike_fs, spike_ts, spike_cluster, cluster_ch] = map_sorter_results(spike_path)
spike_dir = fileparts(spike_path);
load(fullfile(spike_dir, 'mapping.mat'), 'chamap_ns6');
fid = dir(fullfile(spike_dir, '*.prb'));
prb_path = fullfile(fid(1).folder, fid(1).name);

fid = fopen(prb_path, 'r');
prb_str = fread(fid, inf, '*char')';
fclose(fid);

ch_str = regexp(prb_str, "(?<='channels':)\S*(?=,\s*'graph')", 'match');
ngrp = length(ch_str);
[spike_ts, spike_cluster, cluster_ch, is_good] = deal(cell(ngrp, 1));

spike_fs = h5readatt(spike_path, '/recordings/0/', 'sample_rate');
offset = 0;
for gg = 1:ngrp
    if isfile(strrep(fullfile(spike_dir, '.phy', num2str(gg-1), 'memcache', ...
            'phy.apps.kwik.gui.get_best_channels.mat'), '\', '\\'))
        load(strrep(fullfile(spike_dir, '.phy', num2str(gg-1), 'memcache', ...
            'phy.apps.kwik.gui.get_best_channels.mat'), '\', '\\'), ...
            'cluster_channel');
    elseif isfile(strrep(fullfile(spike_dir, '.phy', 'memcache', ...
            'phy.apps.kwik.gui.get_best_channels.mat'), '\', '\\'))
        load(strrep(fullfile(spike_dir, '.phy', 'memcache', ...
            'phy.apps.kwik.gui.get_best_channels.mat'), '\', '\\'), ...
            'cluster_channel');
    else
        continue;
    end
    orig = [ones(size(cluster_channel, 1), 1) * (gg-1), cluster_channel(:, 1)];
    eval(sprintf('[~, iloc] = ismember(%s, chamap_ns6(:, 2));', ch_str{gg}));
    chid = chamap_ns6(iloc(cluster_channel(:, 2) + 1), 1);
    if cluster_channel(1) == 0
         cluster_channel(1) = 1;
    end  
    cluster_ch{gg} = double(cat(2, cluster_channel(:, 1)+offset, chid, orig)); 

    spike_ts{gg} = double(h5read(spike_path, ...
        ['/channel_groups/', num2str(gg-1), '/spikes/time_samples'])) / spike_fs;
    clus = double(h5read(spike_path, ...
        ['/channel_groups/', num2str(gg-1), '/spikes/clusters/main']));
    clus(clus == 0) = 1;
    spike_cluster{gg} = clus + offset;

    info = h5info(spike_path, ['/channel_groups/', num2str(gg-1), '/clusters/main']);
    label = arrayfun(@(x)(x.Attributes.Value), info.Groups);
    cid = arrayfun(@(x)(str2double(regexp(x.Name, ...
        '(?<=/main/)\d+', 'match', 'once'))), info.Groups);
    [is_good{gg}, ic] = ismember(orig(:, 2), cid);
    is_good{gg}(ic > 0) = label(ic(ic > 0)) > 0;

    offset = offset + size(orig, 1);
end
spike_ts = cat(1, spike_ts{:});
spike_cluster = cat(1, spike_cluster{:});
cluster_ch = cat(1, cluster_ch{:});

nspk = arrayfun(@(x)(sum(spike_cluster == x)), cluster_ch(:, 1));
is_good = cat(1, is_good{:}) ~= 0 & nspk > 0;

good_clu = cluster_ch(is_good, 1);
[~, iloc] = ismember(spike_cluster, good_clu);
good_spk = iloc ~= 0;
spike_cluster = iloc(good_spk);
spike_ts = spike_ts(good_spk);
cluster_ch = cluster_ch(is_good, :);
cluster_ch(:, 1) = (1:1:length(good_clu))';
end