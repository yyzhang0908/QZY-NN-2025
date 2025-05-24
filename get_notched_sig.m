function fsig = get_notched_sig(sig, fs)
ord = 6;
F0 = 50;
Q = 20;
method = 'butter';
h = design(fdesign.notch(ord, F0, Q, fs), method);
fsig = filtfilt(h.SOSMatrix, h.ScaleValues, sig);
end