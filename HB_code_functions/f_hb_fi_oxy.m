function [FI_peaks_yes, FI_peak_valid, FI_peak_index, FI_mean_peak_width, FI_mean_peak_height, ...
    FI_peak_area_sum, FI_peak_duration_mean, FI_peak_duration_sum, ...
    pks, locs, w, p] = f_hb_fi_oxy(satSignal, satFreq, ...
    FI_peak_prominence, FI_peak_distance, ...
    FI_time_wakesat, FI_TT, FI_treshold_deviation)
% This function calculates different hypoxic burden parameters 
% based on peak detection according to Filchenko method

%% Define baseline
FI_peaks_yes = [];
FI_peak_valid = [];
FI_peak_index = [];
FI_mean_peak_width = [];
FI_mean_peak_height = [];
FI_peak_area_sum = [];
FI_peak_duration_mean = [];
FI_peak_duration_sum = [];

satSignal_nonnan = satSignal(~isnan(satSignal));

% Calculate the duration of the recording
duration_rec = (length(satSignal_nonnan)/satFreq)/(60*60); % in hours
duration_rec_original = (length(satSignal)/satFreq)/(60*60); % in hours

%% Calculate AUC around peaks
% Find peaks
neg_satSignal = -satSignal;
[pks, locs, w, p] = findpeaks(neg_satSignal, satFreq, 'MinPeakDistance', FI_peak_distance, 'MinPeakProminence', FI_peak_prominence);
% locs are in seconds from the start since it is specified

% ! Continue with the rest of the code only if any peaks were detected
if length(pks) > 0

    % Set the marker for peak presence
    FI_peaks_yes = 1;

    % Calculate basic peak parameters
    npks = length(pks);
    FI_peak_index = npks/duration_rec; % peak index
    FI_mean_peak_width = mean(w); % seconds
    FI_mean_peak_height = mean(p);

    % Calculate AUC
    [FI_peak_valid, FI_peak_area_sum, FI_peak_duration_mean, FI_peak_duration_sum] = f_auc_fi_peaks_calc(satSignal, ...
        satFreq, locs, FI_TT, FI_time_wakesat, FI_treshold_deviation);

else
    disp("No peaks detected")
    % Set the marker for peak presence
    FI_peaks_yes = 0;
    FI_peak_valid = 0;
    FI_peak_index = 0;
    FI_mean_peak_width = 0;
    FI_mean_peak_height = 0;
    FI_peak_area_sum = 0;
    FI_peak_duration_mean = 0;

end

end
