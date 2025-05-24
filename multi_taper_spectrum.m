function [pxxs, fxx, txx] = multi_taper_spectrum(sig, twin, tstep, fs, freqr, ts)
if nargin<5 || isempty(freqr)
    freqr = [-inf, inf];
end

if nargin<6 || isempty(ts)
    toff = 0;
else
    toff = ts(1);
end

nt = fix(length(sig)/fs);
ntx = fix(nt/twin);
nts = ntx*fs*twin;
nstep = round(twin/tstep);
tstep = twin/nstep;
sig = cat(1, sig, zeros(nts+ceil((nstep-1)*tstep*fs)-length(sig), 1));
txx = (1:ntx*nstep)*tstep+twin/2-tstep-toff;

nfft = 2^nextpow2(twin*fs);
fxx = linspace(0, fs/2, nfft/2+1);
infreq = fxx>=freqr(1) & fxx<=freqr(2);
fxx = fxx(infreq);
nfreq = sum(infreq);

tapers = permute(dpss(twin*fs, 3, 5)*sqrt(fs), [1, 3, 2]);
pxxs = zeros(nfreq, nstep, ntx);
for kk = 1:nstep
    idx = (1:nts)+round((kk-1)*tstep*fs);
    pxx = fft(reshape(sig(idx), [], ntx).*tapers, nfft)/fs;
    pxx = pxx(infreq,:,:);
    pxxs(:,kk,:) = mean(conj(pxx).*pxx, 3);
end
pxxs = pxxs(:,:);

end