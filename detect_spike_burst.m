function [is_burst, idx_beg, idx_end] = detect_spike_burst(spkt, ...
    max_init_isi, max_final_isi, min_num_spk, min_pre_isi)
if max_final_isi >= min_pre_isi
    error('The maximum interval among burst spikes must smaller than hyperpolarization period.');
end
nspk = length(spkt);
isit = diff(spkt);

idx_beg = find([nan; isit] >= min_pre_isi/1000 & [isit; nan] < max_init_isi/1000);
np = length(idx_beg);
idx_end = nan(np, 1);
for ii = 1:np
    if idx_beg(ii) > (nspk - min_num_spk)
        continue;
    end
    for jj = idx_beg(ii)+1:nspk-1
        if isit(jj) >= max_final_isi/1000
            break;
        end
    end
    idx_end(ii) = jj;
end
isin = (idx_end - idx_beg + 1) >= min_num_spk;
idx_beg = idx_beg(isin);
idx_end = idx_end(isin);

is_burst = false(size(spkt));
if any(isin)
    idx = arrayfun(@(x, y)(x:1:y), idx_beg, idx_end, 'UniformOutput', false);
    idx = cat(2, idx{:})';
    is_burst(idx) = true;
end
end