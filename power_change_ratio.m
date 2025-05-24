num_ds = length(stats);
Ca2BandpowerBurst_ch = cell(num_ds, 1);
for nn = 1:num_ds
   Ca2SpectrumBurstCh = (stats(nn).Ca2SpectrumInterPer)./(stats(nn).Ca2SpectrumBSLPer)-1;
   Ca2BandpowerBurst_ch = integrate_psd_bandpower(Ca2SpectrumBurstCh, ca2_fxx, opt.ca2_band_freq); 
   A{nn, 1} = Ca2BandpowerBurst_ch;
end
%% 
num_ds = length(stats);
bandpower_anes_ch = cell(num_ds, 1);
for nn = 1:num_ds
   spectrum_anes_ch = (stats(nn).SpectrumAnes)./(stats(nn).SpectrumBSL)-1;
   bandpower_anes_ch = integrate_psd_bandpower(spectrum_anes_ch, spectrum_fxx, opt.band_freq); 
   A{nn, 1} = bandpower_anes_ch;
end