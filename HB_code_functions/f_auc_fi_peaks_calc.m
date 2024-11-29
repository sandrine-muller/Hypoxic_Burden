function [FI_peak_valid, FI_peak_area_sum, FI_peak_duration_mean, FI_peak_duration_sum] = f_auc_fi_peaks_calc(satSignal, satFreq, locs, FI_TT, FI_time_wakesat, FI_treshold_deviation)
% This function is used to detect and describe signal peaks
% according to Filchenko

satSignal_nonnan = satSignal(~isnan(satSignal));

% Calculate the duration of the recording
duration_rec = length(satSignal_nonnan)/satFreq/60/60;
npks = height(locs);

% Find wake saturation 
treshold_mid = mode(satSignal_nonnan(1:(FI_time_wakesat*satFreq)));
treshold_min = treshold_mid - FI_treshold_deviation;
treshold_max = treshold_mid + FI_treshold_deviation;

% Convert locs to signal points
locs_point = locs*satFreq;

% In TT-sec window find the closest treshold_mid value to the peak
search_window = FI_TT*satFreq;

% In TT-sec window find the treshold value to the peak 
peak_treshold = [];
for ll = 1:npks
    sig_start = locs_point(ll) - search_window;
    if sig_start < 1
        sig_start = 1;
    else
        %
    end
    sig_end = locs_point(ll) + search_window;
    if sig_end > length(satSignal)
        sig_end = length(satSignal);
    else
        %
    end
    satTemp = satSignal(sig_start:sig_end);
    satTemp = satTemp(satTemp >= treshold_min & treshold_max >= satTemp);
    peak_treshold_this_peak = mode(satTemp(~isnan(satTemp)));
    peak_treshold = [peak_treshold; peak_treshold_this_peak];
end

% Calculate AUC per peak
peak_auc = [];
peak_duration = [];
for ll = 1:npks
    % Extract the signal in a window around the peak
    sig_start = locs_point(ll) - search_window;
    if sig_start < 1
        sig_start = 1;
    else
        %
    end
    sig_end = locs_point(ll) + search_window;
    if sig_end > length(satSignal)
        sig_end = length(satSignal);
    else
        %
    end

    satTemp = satSignal(sig_start:sig_end);
    satTemp_peak_loc = search_window+1;

    % Find the closest values over the peak_treshold to the peak location
    satTemp_before_peak = satTemp(1:(satTemp_peak_loc-1));

    start_point_to_peak = find(satTemp_before_peak >= peak_treshold(ll));
    start_point_to_peak = max(start_point_to_peak);

    satTemp_after_peak = satTemp((satTemp_peak_loc+1):length(satTemp));

    end_point_from_peak = find(satTemp_after_peak >= peak_treshold(ll));
    end_point_from_peak = min(end_point_from_peak)+search_window+1;

    % Select signal
    satTemp_calc = satTemp(start_point_to_peak:end_point_from_peak);
    if anynan(satTemp_calc) == 1
        % If there is an NaN, we do not assess this peak
        % Calculate AUC and duration
        this_peak_auc = NaN;
        this_peak_duration = NaN;
    else
        % Calculate AUC and duration
        this_peak_auc = trapz(1/(satFreq*60),-(satTemp_calc-peak_treshold(ll)));
        this_peak_duration = length(satTemp_calc)/satFreq;
    end

    % Add to common matrix
    peak_auc = [peak_auc;this_peak_auc];
    peak_duration = [peak_duration;this_peak_duration];
end

% Save the values with relation to the signal duration
FI_peak_area_sum = nansum(peak_auc)/duration_rec; %
FI_peak_duration_mean = nanmean(peak_duration(peak_duration > 0)); % seconds
FI_peak_duration_sum = (nansum(peak_duration)/60/60)/duration_rec*100; % recording time
FI_peak_valid = length(~isnan(peak_duration));

end