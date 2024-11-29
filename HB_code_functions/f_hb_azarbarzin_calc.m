function [AZB_AHI, AZB_nevents_max, AZB_max_area_sum, AZB_max_duration_mean, ...
    AZB_max_duration_sum, AZB_nevents_mode, AZB_mode_area_sum, AZB_mode_duration_mean, ...
    AZB_mode_duration_sum, Error_AZB_val] = f_hb_azarbarzin_calc(satSignal, satFreq, respEventList, ...
    AZB_window_find_max_beforeRE_ends, AZB_minimal_distance_between_events, AZB_window_beforeRE, ...
    AZB_window_afterRE, AZB_mode_low, AZB_mode_high, AZB_pp)
% This function is used to calculate the hypoxic burden
% the method according to Azabarzin

AZB_AHI = [];
AZB_nevents_max = [];
AZB_max_area_sum = [];
AZB_max_duration_mean = [];
AZB_max_duration_sum = [];
AZB_nevents_mode = [];
AZB_mode_area_sum = [];
AZB_mode_duration_mean = [];
AZB_mode_duration_sum = [];
Error_AZB_val = [];

satSignal_nonnan = satSignal(~isnan(satSignal));

% Calculate the duration of the recording
duration_rec = length(satSignal_nonnan)/satFreq/60/60; % in hours
duration_rec_original = length(satSignal)/satFreq/60/60; % in hours

% Merge those respEvents together, which are closer
% than the minimal_distance_between_events
respEventList_merged = respEventList;
tresh_val = 0;

n_events = height(respEventList_merged);

if n_events > 1
    % Loop over single events until no events are left with the short gap
    while tresh_val < AZB_minimal_distance_between_events
        n_events = height(respEventList_merged);
        for rr = 1:(n_events-1)
            respEventList_merged.dif_val(rr) = respEventList_merged.ends_relative_sec(rr+1) - respEventList_merged.starts_relative_sec(rr);
        end
        for rr = 1:(n_events-1)
            if respEventList_merged.dif_val(rr) < AZB_minimal_distance_between_events
                respEventList_merged.ends_relative_sec(rr) = respEventList_merged.ends_relative_sec(rr+1);
                respEventList_merged.delete(rr+1) = 1;
            else
                respEventList_merged.delete(rr+1) = 0;
            end
        end
        respEventList_merged.respEvduration_rec = respEventList_merged.ends_relative_sec - respEventList_merged.starts_relative_sec;
        respEventList_merged = respEventList_merged(respEventList_merged.delete == 0,:);
        n_events = height(respEventList_merged)-1;
        for rr = 1:(n_events-1)
            respEventList_merged.dif_val(rr) = respEventList_merged.ends_relative_sec(rr+1) - respEventList_merged.starts_relative_sec(rr);
        end
        tresh_val = min(respEventList_merged.dif_val(1:n_events));
    end % while
else
    respEventList_merged.delete = 0;
    respEventList_merged.dif_val = AZB_window_find_max_beforeRE_ends-1;
end

n_events = height(respEventList_merged);

% Create a vector with the points of satSignal, where the event ends
endVecNew = int64(respEventList.ends_relative_sec(1:height(respEventList_merged))*satFreq);
startVecNew = int64(respEventList.starts_relative_sec(1:height(respEventList_merged))*satFreq);

% Select only those events that are within the saturation signal
keep_events = endVecNew<length(satSignal);
endVecNew = endVecNew(keep_events);
startVecNew = startVecNew(keep_events);
respEventList_merged = respEventList_merged(keep_events,:);

% Calculate AHI
AZB_AHI = height(respEventList_merged)/duration_rec_original;

% Find treshold_vals (max und mode)
% in AZB_window_find_max_beforeRE_ends for each respEvent
AZB_max_treshold = [];
if n_events == 1;
    rr = 1;
    if respEventList_merged.starts_relative_sec > AZB_window_find_max_beforeRE_ends;
        satPosit_start = (endVecNew(rr)-AZB_window_find_max_beforeRE_ends*satFreq);
        satPosit_end = endVecNew(rr);
        if satPosit_start < 1 | satPosit_end > length(satSignal)
            Error_AZB_val = 1;
        else
            satSignTmp = satSignal(satPosit_start:satPosit_end);
            max_sat = nanmax(satSignTmp);
            AZB_max_treshold = [AZB_max_treshold; max_sat];
        end
    else
        satPosit_start = 1;
        satPosit_end = endVecNew(rr);
        if satPosit_end > length(satSignal)
            max_sat = NaN;

        else
            satSignTmp = satSignal(satPosit_start:satPosit_end);
            max_sat = nanmax(satSignTmp);
        end
        AZB_max_treshold = [AZB_max_treshold; max_sat];
    end
else
    for rr = 1:height(respEventList_merged)
        if respEventList_merged.starts_relative_sec(rr) > AZB_window_find_max_beforeRE_ends
            satPosit_start = (endVecNew(rr)-AZB_window_find_max_beforeRE_ends*satFreq);
            satPosit_end = endVecNew(rr);
            % Check if the position of off limits
            if satPosit_start < 1 | satPosit_end > length(satSignal)
                max_sat = NaN;
            else
                satSignTmp = satSignal(satPosit_start:satPosit_end);
                max_sat = nanmax(satSignTmp);
            end
            AZB_max_treshold = [AZB_max_treshold; max_sat];
        else
            % Check if the position of off limits
            satPosit_start = 1;
            satPosit_end = endVecNew(rr);
            if satPosit_end > length(satSignal)
                max_sat = NaN;
            else
                satSignTmp = satSignal(satPosit_start:satPosit_end);
                max_sat = nanmax(satSignTmp);
            end
            AZB_max_treshold = [AZB_max_treshold; max_sat];
        end
    end
end
% if there are NaNs in the max_treshold, these events could not be assessed
AZB_nevents_max = length(~isnan(AZB_max_treshold));

AZB_mode_treshold = [];
if n_events == 1;
    rr = 1;
    if respEventList_merged.starts_relative_sec > AZB_window_find_max_beforeRE_ends;
        satPosit_start = (endVecNew(rr)-AZB_window_find_max_beforeRE_ends*satFreq);
        satPosit_end = endVecNew(rr);
        if satPosit_start < 1 | satPosit_end > length(satSignal)
            Error_AZB_val = 1;
        else
            satSignTmp = satSignal(satPosit_start:satPosit_end);
            mode_sat = f_hb_azabarzin_find_mode(satSignTmp, AZB_mode_low, AZB_mode_high);
            AZB_mode_treshold = [AZB_mode_treshold; mode_sat];
        end
    else
        satPosit_start = 1;
        satPosit_end = endVecNew(rr);
        if satPosit_end > length(satSignal)
            mode_sat = NaN;

        else
            satSignTmp = satSignal(satPosit_start:satPosit_end);
            mode_sat = mode(satSignTmp);
        end
        AZB_mode_treshold = [AZB_mode_treshold; mode_sat];
    end
else
    for rr = 1:height(respEventList_merged)
        if respEventList_merged.starts_relative_sec(rr) > AZB_window_find_max_beforeRE_ends
            satPosit_start = (endVecNew(rr)-AZB_window_find_max_beforeRE_ends*satFreq);
            satPosit_end = endVecNew(rr);
            % Check if the position of off limits
            if satPosit_start < 1 | satPosit_end > length(satSignal)
                mode_sat = NaN;
            else
                satSignTmp = satSignal(satPosit_start:satPosit_end);
                mode_sat = f_hb_azabarzin_find_mode(satSignTmp, AZB_mode_low, AZB_mode_high);
            end
            AZB_mode_treshold = [AZB_mode_treshold; mode_sat];
        else
            % Check if the position of off limits
            satPosit_start = 1;
            satPosit_end = endVecNew(rr);
            if satPosit_end > length(satSignal)
                mode_sat = NaN;
            else
                satSignTmp = satSignal(satPosit_start:satPosit_end);
                mode_sat = f_hb_azabarzin_find_mode(satSignTmp, AZB_mode_low, AZB_mode_high);
            end
            AZB_mode_treshold = [AZB_mode_treshold; mode_sat];
        end
    end
end

AZB_mode_treshold = [];
for rr = 1:height(respEventList_merged)
    if respEventList_merged.starts_relative_sec > AZB_window_find_max_beforeRE_ends
        satPositTemp = (endVecNew(rr)-AZB_window_find_max_beforeRE_ends*satFreq):endVecNew(rr);
        satSignTmp = satSignal(satPositTemp);
        if range(satSignTmp) > 0 & length(satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp)) > 1;
            % Find mode only in the physiological range
            satSignTmp = satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp);
            max_sat = mode(satSignTmp(~isnan(satSignTmp)));
        else
            max_sat = nanmean(satSignTmp);
        end
        AZB_mode_treshold = [AZB_mode_treshold; max_sat];
    else
        satPositTemp = 1:endVecNew(rr);
        satSignTmp = satSignal(satPositTemp);
        if range(satSignTmp) > 0 & length(satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp)) > 1;
            % Find mode only in the physiological range
            satSignTmp = satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp);
            max_sat = mode(satSignTmp(~isnan(satSignTmp)));
        else
            max_sat = nanmean(satSignTmp);
        end
        AZB_mode_treshold = [AZB_mode_treshold; max_sat];
    end
end
% if there are NaNs in the mode_treshold, these events could not be assessed
AZB_nevents_mode = length(~isnan(AZB_mode_treshold));

% Calculate AUC based on treshold for baseline
if AZB_nevents_max > 1
    Error_AZB_val = 0;
    % 1) Use maximal values as baseline
    treshold_vals = AZB_max_treshold;
    [AZB_max_area_sum, AZB_max_duration_mean, AZB_max_duration_sum] = f_auc_azarbarzin_calc(treshold_vals, satSignal, satFreq, respEventList_merged, startVecNew, endVecNew, AZB_window_beforeRE, AZB_window_afterRE, AZB_pp);
    % 2) Use mode values as baseline
    treshold_vals = AZB_mode_treshold;
    [AZB_mode_area_sum, AZB_mode_duration_mean, AZB_mode_duration_sum] = f_auc_azarbarzin_calc(treshold_vals, satSignal, satFreq, respEventList_merged, startVecNew, endVecNew, AZB_window_beforeRE, AZB_window_afterRE, AZB_pp);
else
    Error_AZB_val = 1;
end

end


