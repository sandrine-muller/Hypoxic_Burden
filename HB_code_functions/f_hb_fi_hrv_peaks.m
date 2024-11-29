function [peak_hrSegDuration, peak_HR_loc, peak_rrHRV_loc, ...
    peak_SDNN_loc, peak_SDSD_loc, peak_RMSSD_loc, ...
    peak_pNN50_loc, peak_TRI_val_loc, peak_TINN_val_loc, ...
    peak_DFA_1_loc, peak_DFA_2_loc, peak_ApEn_loc, ...
    peak_fvf_pLF, peak_fvf_pHF, peak_fvf_LFHFratio, ...
    peak_fvf_VLF, peak_fvf_LF, peak_fvf_HF, ...
    peak_SD1, peak_SD2, peak_SD1SD2ratio] = f_hb_fi_hrv_peaks(pks, locs, ...
    FI_peak_prominence, FI_peak_distance, ...
    FI_peak_event_before, FI_peak_event_after, FI_time_wakesat, ...
    FI_TT, FI_treshold_deviation, ...
    hrFreq, hrSignal, hr_type, hr_qrs_loc, ...
    FI_hrv_peaks_before, FI_hrv_peaks_after)
% This function calculates HRV around peaks according to Filchenko method

%% Define baseline
peak_hrSegDuration = [];
peak_HR_loc = [];
peak_rrHRV_loc = [];
peak_SDNN_loc = [];
peak_SDSD_loc = [];
peak_RMSSD_loc = [];
peak_pNN50_loc = [];
peak_TRI_val_loc = [];
peak_TINN_val_loc = [];
peak_DFA_1_loc = [];
peak_DFA_2_loc = [];
peak_ApEn_loc = [];
peak_fvf_pLF = [];
peak_fvf_pHF = [];
peak_fvf_LFHFratio = [];
peak_fvf_VLF = [];
peak_fvf_LF = [];
peak_fvf_HF = [];
peak_SD1 = [];
peak_SD2 = [];
peak_SD1SD2ratio = [];

% Calculate the duration of the recording
duration_rec = (length(hrSignal)/hrFreq)/(60*60); % in hours
duration_rec_original = (length(hrSignal)/hrFreq)/(60*60); % in hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Calculate HRV around peaks
% Only if there is an HR channel
npks = length(pks);
if npks > 1
    % Calculate HRV per segment around each peak
    hrv_output_array_all = [];
    for ll = 1:npks
        start_hr = locs(ll) - FI_hrv_peaks_before;
        if start_hr < 1
            start_hr = 1;
        else
            %
        end
        end_hr = locs(ll) + FI_hrv_peaks_after;
        if end_hr > (length(hrSignal)/hrFreq)
            end_hr = length(hrSignal)/hrFreq;
        else
            %
        end
        start_hr = int64(start_hr*hrFreq);
        end_hr = int64(end_hr*hrFreq);
        hrSeg = hrSignal(start_hr:end_hr);
        [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
        hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
    end
    hrv_output_array_mean = nanmean(hrv_output_array_all);
else
    hrv_output_array_all = [];

    start_hr = locs - FI_hrv_peaks_before;
    if start_hr < 1
        start_hr = 1;
    else
        %
    end
    end_hr = locs + FI_hrv_peaks_after;
    if end_hr > (length(hrSignal)/hrFreq)
        end_hr = length(hrSignal)/hrFreq;
    else
        %
    end
    start_hr = int64(start_hr*hrFreq);
    end_hr = int64(end_hr*hrFreq);
    hrSeg = hrSignal(start_hr:end_hr);
    [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc);
    hrv_output_array_all = [hrv_output_array_all;hrv_output_array];
    hrv_output_array_mean = hrv_output_array_all;
end

% Set the values
peak_hrSegDuration = hrv_output_array_mean(1);
peak_HR_loc = hrv_output_array_mean(2);
peak_rrHRV_loc = hrv_output_array_mean(3);
peak_SDNN_loc = hrv_output_array_mean(4);
peak_SDSD_loc = hrv_output_array_mean(5);
peak_RMSSD_loc = hrv_output_array_mean(6);
peak_pNN50_loc = hrv_output_array_mean(7);
peak_TRI_val_loc = hrv_output_array_mean(8);
peak_TINN_val_loc = hrv_output_array_mean(9);
peak_DFA_1_loc = hrv_output_array_mean(10);
peak_DFA_2_loc = hrv_output_array_mean(11);
peak_ApEn_loc = hrv_output_array_mean(12);
peak_fvf_pLF = hrv_output_array_mean(13);
peak_fvf_pHF = hrv_output_array_mean(14);
peak_fvf_LFHFratio = hrv_output_array_mean(15);
peak_fvf_VLF = hrv_output_array_mean(16);
peak_fvf_LF = hrv_output_array_mean(17);
peak_fvf_HF = hrv_output_array_mean(18);
peak_SD1 = hrv_output_array_mean(19);
peak_SD2 = hrv_output_array_mean(20);
peak_SD1SD2ratio = hrv_output_array_mean(21);

end
