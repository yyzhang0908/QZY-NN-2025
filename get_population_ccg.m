function [C, Z, P, T] = get_population_ccg(spkt, mlag, dt, jw, trg, niter)
if nargin < 5 || isempty(trg); trg = [0, inf]; end
if nargin < 6; niter = 100; end

mskspkt = cellfun(@(x)(mask_time_range(x, trg, x)), spkt, 'UniformOutput', false);
nu = length(mskspkt);
uns = cellfun(@(x)(numel(x)), mskspkt);

T = 0:dt:mlag;
T = [-T(end:-1:2), T];
tw = (min(T)-dt/2):dt:(max(T)+dt/2);
nt = numel(T);

[C, Z, P] = deal(nan(nu, nu, nt));

for ii = 1:nu
    spkt1 = mskspkt{ii} * 1000;
    for jj = 1:nu
        if ii == jj; continue; end
        spkt2 = mskspkt{jj} * 1000;

        K = ccg_raw(spkt1, spkt2, tw);
        Ks = zeros(niter, nt);
        rng(1024);
        parfor kk = 1:niter
            jitt = (rand(uns(jj), 1) - 0.5) * 2 * 1000 * jw;
            Ks(kk, :) = ccg_raw(spkt1, spkt2 + jitt, tw);
        end

        C(ii, jj, :) = K / uns(ii) / uns(jj) * 1000;
        Z(ii, jj, :) = (K - mean(Ks, 1)) ./ std(Ks, 1, 1);
        P(ii, jj, :) = mean((K - Ks) >= 0, 1) * 100;
    end
end

C = reshape(C, [], nt);
Z = reshape(Z, [], nt);
Z(~isfinite(Z)) = nan;
P = reshape(P, [], nt);
end

function cch = ccg_raw(st1, st2, t)
nt = numel(t);
cch = zeros(1, nt - 1);
n2 = numel(st2);
for nn = 1:n2
    cch = cch + histcounts(st1 - st2(nn), t);
end
end