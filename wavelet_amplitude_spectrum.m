function [pxx, fxx] = wavelet_amplitude_spectrum(sig, ts, fs, freqr, nscales)
ts_new = (min(ts):1/round(fs):max(ts))';
sig_new = interp1(ts, sig, ts_new, 'spline', 0);
[pxx, fxx] = cwt(sig_new, 'amor', round(fs), 'FrequencyLimits', freqr, ...
    'VoicesPerOctave', nscales);
fxx = flipud(fxx);
pxx = flipud(abs(pxx))';
pxx = interp1(ts_new, pxx, ts, 'spline', 0);
pxx = (pxx).^2;
end