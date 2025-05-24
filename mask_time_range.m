function out = mask_time_range(ts, trange, val)
if isempty(ts); out = []; end

if isempty(trange); trange = [0, inf]; end
nrg = size(trange, 1);

tmask = false(size(ts));
for nn = 1:nrg
    tmask(ts >= trange(nn, 1) & ts <= trange(nn, 2)) = true;
    if nargin == 3
        if isempty(val)
            out = [];
        else
            out = val(tmask, :);
        end
    else
        out = tmask;
    end
end