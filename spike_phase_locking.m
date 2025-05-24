function [spk_phase, phase_index] = spike_phase_locking(spkt, sig, ts, tev)
phis = angle(hilbert(sig));
is_win = false(size(spkt));
for kk = 1:size(tev, 1)
    is_win = (spkt>tev(kk, 1) & spkt<tev(kk, 2)) | is_win;
end
spk_phase = interp1(ts, phis, spkt(is_win));
[mrl, raylz, pvalue, mu, kappa] = circular_statistics(spk_phase);
phase_index = [mrl, raylz, pvalue, mu, kappa];
end