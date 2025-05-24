function [lag_corrs, lag_ts, pk_cc, pk_lag] = time_mask_crosscorr( ...
    ref_sig, lag_sig, lag_max, fs, t_mask)
lags = (round(-lag_max*fs):1:round(lag_max*fs))';
lag_ts = lags / fs;
lag_n = length(lags);
lag_corrs = zeros(lag_n, size(lag_sig, 2));
if length(find(diff([false; t_mask])==1))~=1 
    ref_sig(~t_mask) = nan;
    lag_sig(~t_mask, :) = nan;
    parfor ii = 1:lag_n
        if lags(ii) >= 0
            lag_corrs(ii, :) = corr(ref_sig(lags(ii)+1:end), ...
                lag_sig(1:end-lags(ii), :), 'rows', 'pairwise');
        else
            lag_corrs(ii, :) = corr(ref_sig(1:end+lags(ii)), ...
                lag_sig(-lags(ii)+1:end, :), 'rows', 'pairwise');
        end
    end
else
    is_out = ref_sig>prctile(ref_sig(t_mask), 99.9) | ...
        ref_sig<prctile(ref_sig(t_mask), 0.1);
    ref_sig(is_out) = median(ref_sig(t_mask));
    for bb = 1:size(lag_sig, 2) 
        tmp = lag_sig(t_mask, bb);
        is_out = tmp>prctile(tmp, 99.9) | tmp<prctile(tmp, 0.1);
        tmp(is_out) = median(tmp);
        lag_corrs(:, bb) = xcov(ref_sig(t_mask), tmp, ...
            round(lag_max * fs), 'coeff');
    end
end
lag_corrs = flip(lag_corrs, 1);

[pk_cc, pk_lag] = deal(zeros(1, size(lag_sig, 2)));
for bb = 1:size(lag_sig, 2) 
    [~, ipk] = max(lag_corrs(:, bb));
    pk_cc(bb) = lag_corrs(ipk, bb);
    pk_lag(bb) = lag_ts(ipk);
end
end

