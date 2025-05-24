function [pxxs, fxx, txx] = step_window_spectrum(sig, twin, tstep, fs, nwin, freqr, ts)
if nargin<6 || isempty(freqr)
    freqr = [-inf, inf];
end

if nargin<7 || isempty(ts)
    toff = 0;
else
    toff = ts(1);
end

nt = fix(length(sig)/fs);
ntx = fix(nt/twin);
nts = ntx*fs*twin;
nstep = round(twin/tstep);
tstep = twin/nstep;
sig = cat(1, sig, zeros(nts+(nstep-1)*tstep*fs-length(sig), 1));
txx = (1:ntx*nstep)*tstep+twin/2-tstep-toff;

nfft = 2^nextpow2(twin*fs);
fxx = linspace(0, fs/2, nfft/2+1);
infreq = fxx>=freqr(1) & fxx<=freqr(2);
fxx = fxx(infreq);
nfreq = sum(infreq);

npad = round((twin*fs-nwin)/2);
winv =  padarray(rectwin(nwin), npad, 0, 'both');

pxxs = zeros(nfreq, nstep, ntx);
parfor kk = 1:nstep
    idx = (1:nts)+round((kk-1)*tstep*fs);
    pxx = pwelch(reshape(sig(idx), [], ntx), winv, 0, nfft, fs);
    pxxs(:,kk,:) = pxx(infreq, :);
end
pxxs = pxxs(:,:);

end