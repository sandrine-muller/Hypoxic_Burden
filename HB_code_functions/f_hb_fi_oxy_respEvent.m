function [FI_AHI, FI_peaks_respEvent_yes, FI_peak_respEvent_valid, FI_peak_respEvent_index, ...
    FI_mean_peak_respEvent_width, FI_mean_peak_respEvent_height, ...
    FI_peak_respEvent_area_sum, FI_peak_respEvent_duration_mean, ...
    FI_peak_respEvent_duration_sum] = f_hb_fi_oxy_respEvent(satSignal, satFreq, respEventList, ...
    FI_minimal_distance_between_events, FI_peak_prominence, FI_peak_distance, ...
    FI_peak_event_before, FI_peak_event_after, FI_time_wakesat, ...
    FI_TT, FI_treshold_deviation, pks, locs, w, p)
% This function calculates different hypoxic burden parameters 
% based on peaks around respiratory events according to Filchenko method

%% Define baseline
FI_AHI = [];
FI_peaks_respEvent_yes = [];
FI_peak_respEvent_valid = [];
FI_peak_respEvent_index = [];
FI_mean_peak_respEvent_width = [];
FI_mean_peak_respEvent_height = [];
FI_peak_respEvent_area_sum  = [];
FI_peak_respEvent_duration_mean = [];
FI_peak_respEvent_duration_sum = [];

satSignal_nonnan = satSignal(~isnan(satSignal));

% Calculate the duration of the recording
duration_rec = (length(satSignal_nonnan)/satFreq)/(60*60); % in hours
duration_rec_original = (length(satSignal)/satFreq)/(60*60); % in hours

%% Calculate AUC around event-related peaks
% Find event-related peaks
% Merge those respEvents together, which are closer
% than the minimal_distance_between_events

npks = length(pks);
% Calculate only if respEventList is not empty
if  npks > 0 & ~isempty(respEventList) & class(respEventList) ~= "double"
    
    % Merge respiratory events
    [respEventList_merged] = f_hb_merge_respEvent_list(respEventList, FI_minimal_distance_between_events);

    % Calculate AHI
    FI_AHI = height(respEventList_merged)/duration_rec_original;

    % Find is the peak is event-related or not
    pks_event_rel = zeros(npks,1);
    for ppp = 1:npks
        for e = 1:height(respEventList_merged)
            if locs(ppp) >= (respEventList_merged.starts_relative_sec(e) - FI_peak_event_before) & ...
                    (respEventList_merged.ends_relative_sec(e) + FI_peak_event_after) >= locs(ppp)
                pks_event_rel(ppp) = 1;
                continue
            else
                %
            end
        end
    end

    % Calculate parameters for respEvent-related peaks
    % It is done only if there are respEvent-related peaks
    % sum(pks_event_rel) > 0
    if sum(pks_event_rel) > 0
        FI_peaks_respEvent_yes = 1;
        pks_respEvent = pks(logical(pks_event_rel));
        locs_respEvent = locs(logical(pks_event_rel));
        w_respEvent = w(logical(pks_event_rel));
        p_respEvent = p(logical(pks_event_rel));

        % Calculate basic respEvent-related peak parameters
        npks_respEvent = length(pks_respEvent);
        FI_peak_respEvent_index = npks_respEvent/duration_rec; % peak index
        FI_mean_peak_respEvent_width = nanmean(w_respEvent);
        FI_mean_peak_respEvent_height = nanmean(p_respEvent);

        % Calculate AUC
        [FI_peak_respEvent_valid, FI_peak_respEvent_area_sum, ...
            FI_peak_respEvent_duration_mean, FI_peak_respEvent_duration_sum] = f_auc_fi_peaks_calc(satSignal, satFreq, ...
            locs_respEvent, FI_TT, FI_time_wakesat, FI_treshold_deviation);

    else
        FI_peak_respEvent_valid = 0;
        FI_peaks_respEvent_yes = 0;
        FI_peak_respEvent_area_sum = NaN;
        FI_peak_respEvent_duration_mean = NaN;
        FI_peak_respEvent_duration_sum = NaN;
    end % respEvent-related peaks detected

else
    FI_peak_respEvent_valid = 0;
    FI_peaks_respEvent_yes = 0;
    FI_peak_respEvent_area_sum = [];
    FI_peak_respEvent_duration_mean = [];
    FI_peak_respEvent_duration_sum = [];
end

end
