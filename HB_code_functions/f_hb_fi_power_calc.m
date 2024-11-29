function [FI_all_bandpower, FI_wd, FI_powtot, FI_pxx] = f_hb_fi_power_calc(satSignal, satFreq, respEventList)
% This function calculates different power parameters 
% according to Filchenko method

satSignal_nonnan = satSignal(~isnan(satSignal));

% Calculate the duration of the recording
duration_rec = (length(satSignal_nonnan)/satFreq)/(60*60); % in hours
duration_rec_original = (length(satSignal)/satFreq)/(60*60); % in hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple power analysis of the signal
% Select only finite elements
isBad=~isfinite(satSignal_nonnan);
satSignal_pa = satSignal_nonnan(~isBad);

% All bandpower
FI_all_bandpower = bandpower(satSignal_pa,satFreq,[0 satFreq/3]); % unit: %^2/Hz

% The width of the frequency band that contains 99% of the power of the signal
[wd,lo,hi,power] = obw(satSignal_pa,satFreq);
FI_wd = wd*1000; % Fs
% Total power
FI_powtot = power/0.99; % unit: %/Hz

% pwelch
segmentLength = 3*satFreq;
noverlap = satFreq;
pxx = pwelch(satSignal_pa,segmentLength,noverlap);
FI_pxx = mean(pxx); % unit: %/Hz

end


