function [respEv_hrSegDuration, respEv_HR_loc, respEv_rrHRV_loc, ...
    respEv_SDNN_loc, respEv_SDSD_loc, respEv_RMSSD_loc, ...
    respEv_pNN50_loc, respEv_TRI_val_loc, respEv_TINN_val_loc, ...
    respEv_DFA_1_loc, respEv_DFA_2_loc, respEv_ApEn_loc, ...
    respEv_fvf_pLF, respEv_fvf_pHF, respEv_fvf_LFHFratio, ...
    respEv_fvf_VLF, respEv_fvf_LF, respEv_fvf_HF, ...
    respEv_SD1, respEv_SD2, ...
    respEv_SD1SD2ratio] = f_hb_fi_hrv_respEvent(respEventList, ...
    FI_minimal_distance_between_events, FI_time_wakesat, ...
    FI_TT, FI_treshold_deviation, ...
    hrFreq, hrSignal, hr_type, hr_qrs_loc, ...
    FI_hrv_respEvent_before, FI_hrv_respEvent_after)
% This function calculates HRV around respEvents according to Filchenko method

%% Define baseline
respEv_hrSegDuration = [];
respEv_HR_loc = [];
respEv_rrHRV_loc = [];
respEv_SDNN_loc = [];
respEv_SDSD_loc = [];
respEv_RMSSD_loc = [];
respEv_pNN50_loc = [];
respEv_TRI_val_loc = [];
respEv_TINN_val_loc = [];
respEv_DFA_1_loc = [];
respEv_DFA_2_loc = [];
respEv_ApEn_loc = [];
respEv_fvf_pLF = [];
respEv_fvf_pHF = [];
respEv_fvf_LFHFratio = [];
respEv_fvf_VLF = [];
respEv_fvf_LF = [];
respEv_fvf_HF = [];
respEv_SD1 = [];
respEv_SD2 = [];
respEv_SD1SD2ratio  = [];

% Calculate the duration of the recording
duration_rec = (length(hrSignal)/hrFreq)/(60*60); % in hours
duration_rec_original = (length(hrSignal)/hrFreq)/(60*60); % in hours

%% Calculate HRV around respEventList_merged
n_events = height(respEventList);

if ~isempty(respEventList) & n_events > 1 & class(respEventList) ~= "double"

    % Merge the respEventList
    [respEventList_merged] = f_hb_merge_respEvent_list(respEventList, FI_minimal_distance_between_events);

    % Calculate HRV per segment around each respEventList_merged
    tresh_val = 0;
    hrv_output_array_all = [];
    for ll = 1:height(respEventList_merged)
        if respEventList_merged.starts_relative_sec(ll) > FI_hrv_respEvent_before
            start_hr = respEventList_merged.starts_relative_sec(ll) - FI_hrv_respEvent_before;
            end_hr = respEventList_merged.ends_relative_sec(ll) + FI_hrv_respEvent_before;
            start_hr = int64(start_hr*hrFreq);
            end_hr = int64(end_hr*hrFreq);
            if end_hr > length(hrSignal)
                end_hr = length(hrSignal);
            else
                %
            end
            hrSeg = hrSignal(start_hr:end_hr);
            [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
            hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
        else
            start_hr = respEventList_merged.starts_relative_sec(ll);
            end_hr = respEventList_merged.ends_relative_sec(ll) + FI_hrv_respEvent_before;
            start_hr = int64(start_hr*hrFreq);
            end_hr = int64(end_hr*hrFreq);
            if end_hr > length(hrSignal)
                end_hr = length(hrSignal);
            else
                %
            end
            hrSeg = hrSignal(start_hr:end_hr);
            [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
            hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
        end
    end
    hrv_output_array_mean = nanmean(hrv_output_array_all);
    % If there is only one event
elseif ~isempty(respEventList) & n_events == 1 & ~isnan(respEventList{1,1})
    %elseif ~isempty(respEventList) & n_events == 1 & ~isnan(respEventList)
    hrv_output_array_all = [];

    respEventList_merged = respEventList;

    if respEventList_merged.starts_relative_sec > FI_hrv_respEvent_before
        start_hr = respEventList_merged.starts_relative_sec - FI_hrv_respEvent_before;
        end_hr = respEventList_merged.ends_relative_sec + FI_hrv_respEvent_before;
        start_hr = int64(start_hr*hrFreq);
        end_hr = int64(end_hr*hrFreq);
        if end_hr > length(hrSignal)
            end_hr = length(hrSignal);
        else
            %
        end
        hrSeg = hrSignal(start_hr:end_hr);
        [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
        hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
    else
        start_hr = respEventList_merged.starts_relative_sec;
        end_hr = respEventList_merged.ends_relative_sec + FI_hrv_respEvent_before;
        start_hr = int64(start_hr*hrFreq);
        end_hr = int64(end_hr*hrFreq);
        if end_hr > length(hrSignal)
            end_hr = length(hrSignal);
        else
            %
        end
        hrSeg = hrSignal(start_hr:end_hr);
        [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
        hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
    end

    hrv_output_array_mean = hrv_output_array_all;
else
    hrv_output_array_mean = cell(1,21);
end

%% Set the values (they are in a specific order)
respEv_hrSegDuration = hrv_output_array_mean(1);
respEv_HR_loc = hrv_output_array_mean(2);
respEv_rrHRV_loc = hrv_output_array_mean(3);
respEv_SDNN_loc = hrv_output_array_mean(4);
respEv_SDSD_loc = hrv_output_array_mean(5);
respEv_RMSSD_loc = hrv_output_array_mean(6);
respEv_pNN50_loc = hrv_output_array_mean(7);
respEv_TRI_val_loc = hrv_output_array_mean(8);
respEv_TINN_val_loc = hrv_output_array_mean(9);
respEv_DFA_1_loc = hrv_output_array_mean(10);
respEv_DFA_2_loc = hrv_output_array_mean(11);
respEv_ApEn_loc = hrv_output_array_mean(12);
respEv_fvf_pLF = hrv_output_array_mean(13);
respEv_fvf_pHF = hrv_output_array_mean(14);
respEv_fvf_LFHFratio = hrv_output_array_mean(15);
respEv_fvf_VLF = hrv_output_array_mean(16);
respEv_fvf_LF = hrv_output_array_mean(17);
respEv_fvf_HF = hrv_output_array_mean(18);
respEv_SD1 = hrv_output_array_mean(19);
respEv_SD2 = hrv_output_array_mean(20);
respEv_SD1SD2ratio = hrv_output_array_mean(21);

end
