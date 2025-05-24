function sync = get_population_jbsi(spkt, sw, trg)
if nargin < 3; trg = [0, inf]; end
mskspkt = cellfun(@(x)(mask_time_range(x, trg, x)), spkt, 'UniformOutput', false);

nu = length(mskspkt);
uns = cellfun(@(x)(numel(x)), mskspkt);

jbsi = deal(nan(nu, nu));
parfor xx = 1:nu
    for yy = 1:nu
        if uns(xx) > uns(yy) || uns(xx) == 0 || xx == yy; continue; end
        jbsi(xx, yy) = jbsi_sync(mskspkt{xx}, mskspkt{yy}, sw);
    end
end

sync = mean(jbsi, 'all', 'omitnan');
end

function [jbsi, zs] = jbsi_sync(a, b, sw)
na = length(a);
nb = length(b);
sw = sw / 1000;
jw = 2 * sw;

t1 = b - sw;
t2 = b + sw;
while true
    isout = find(t1(2:end) <= t2(1:end-1));
    if ~any(isout)
        break;
    end
    t1(isout + 1) = [];
    t2(isout) = [];
end
seg = [t1, t2];
ns = size(seg, 1);

t1 = a - jw;
t2 = a + jw;
psi = zeros(na, 1);
for ii = 1:na
    for jj = 1:ns
        if seg(jj, 2) < t1(ii)
            continue;
        elseif  seg(jj, 1) > t2(ii)
            break;
        else
            psi(ii) = psi(ii) + (min(t2(ii), seg(jj, 2)) - max(t1(ii), seg(jj, 1)));
        end
    end
end
psi = psi / (jw * 2);

nc = 0;
for ii = 1:na
    for jj = 1:nb
        if abs(a(ii) - b(jj)) <= sw
            nc = nc + 1;
            break;
        end
    end
end

ave =  sum(psi);
sd = sqrt(sum(psi .* (1 - psi)));
zs = (nc - ave) / sd;
jbsi = 2 * (nc - ave) / na;
end
