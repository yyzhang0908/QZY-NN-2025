function  bp = integrate_psd_bandpower(sxx, fxx, freqs)
sz = size(sxx);
if sz(1) == length(fxx)
    nt = sz(2);
else
    nt = sz(1);
    sxx = sxx';
end

if length(unique(diff(fxx))) == 1
    method = 'trap1';
else
    method = 'trap2';
end

nb = length(freqs);
bp = zeros(nt, nb);
for bb = 1:nb
%     bp(:, bb) = bandpower(sxx, fxx, freqs{bb}, 'psd');
    isf = fxx >= min(freqs{bb}) & fxx < max(freqs{bb});
    switch method
        case 'trap1'
            bp(:, bb) = trapz(fxx(isf), sxx(isf, :), 1)/ range(freqs{bb}); 
                       
        case 'trap2'
            bp(:, bb) = trapz(fxx(isf), sxx(isf, :), 1) / range(freqs{bb});            
    end
end
end
