function [env, fsig] = get_band_envelope(sig, fs, freqs)
ord = 6;
method = 'butter';
nb = length(freqs);
nt = length(sig);
[env, fsig] = deal(zeros(nt, nb));
for bb = 1:nb
    h = design(fdesign.bandpass('N,F3dB1,F3dB2', ...
        ord, freqs{bb}(1), freqs{bb}(2), fs), method);
    fsig(:, bb) = filtfilt(h.SOSMatrix, h.ScaleValues, sig);
    env(:, bb) = abs(hilbert(fsig(:, bb)));
end
end