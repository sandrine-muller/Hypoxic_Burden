function [hrv_output_array] = f_hb_fi_hrv(hrFreq, hrSeg, hr_type, hr_qrs_loc)
% This function calculates the HRV based on the selected signal segment

% Settings for the automated beat detection
load(hr_qrs_loc)
s = hr_type;
Fs = hrFreq; % set your sampling frequency
Beat_min = qrs_settings.Beat_min(s);
Beat_max = qrs_settings.Beat_max(s);
wl_tma = ceil(qrs_settings.wl_tma(s)*Fs);
wl_we  = ceil(qrs_settings.wl_we(s,:).*Fs);
d_fs = Fs;

% Define sig_waveform
sig_waveform = hrSeg;

% Calculate the duration
hrSegDuration = length(sig_waveform)/Fs;

% Heart beat detection
Ann = [];
seg = ceil(length(sig_waveform)/(300*Fs));
if seg>2
    for i=0:seg
        sig_waveform_tmp = sig_waveform(max(300*Fs*i-10*Fs,1):min(300*Fs*(i+1),length(sig_waveform)));
        if sum(isnan(sig_waveform_tmp)) ~= length(sig_waveform_tmp)
            Ann_tmp = singleqrs(sig_waveform_tmp,Fs,'downsampling',d_fs,'Beat_min',Beat_min,'Beat_max',Beat_max,'wl_tma',wl_tma,'wl_we',wl_we);
            Ann = [Ann; Ann_tmp+max(300*Fs*i-10*Fs,1)];
        end
    end
    Ann = Ann(Ann>0 & Ann<=length(sig_waveform));
    Ann = unique(sort(Ann));
    Ann(diff(Ann)<.05*Fs)=[];
else
    Ann = singleqrs(sig_waveform,Fs,'downsampling',d_fs,'Beat_min',Beat_min,'Beat_max',Beat_max,'wl_tma',wl_tma,'wl_we',wl_we);
end
Ann = Ann/Fs;

% RR intervals and filtering of artifacts
RR = diff(Ann);
RR_filt = HRV.RRfilter(RR,20);
RR_loc = RR_filt;

% HRV cannot be computed
hrv_cannot_be_computed = [];
if length(RR_loc) - sum(isnan(RR_loc)) > 7 
    hrv_cannot_be_computed = 0;
else
    hrv_cannot_be_computed = 1;
end

if hrv_cannot_be_computed == 0; 
    % Compute local HRV measures with NAs
    HR_loc  = HRV.HR(RR_loc,0);
    rrHRV_loc = HRV.rrHRV(RR_loc);
    SDNN_loc  = HRV.SDNN(RR_loc,0)*1000;
    SDSD_loc  = HRV.SDSD(RR_loc,0)*1000;
    RMSSD_loc = HRV.RMSSD(RR_loc,0)*1000;
    pNN50_loc = HRV.pNN50(RR_loc,0)*100;
    TRI_val_loc = HRV.TRI(RR_loc);
    TINN_val_loc = HRV.TINN(RR_loc);

    % Compute local HRV measures without NAs
    % Remove NAs
    nNA_RR_loc = rmmissing(RR_loc);

    % Compute local HRV measures
    DFA_loc = HRV.DFA(nNA_RR_loc);
    DFA_1_loc = DFA_loc(1);
    DFA_2_loc = DFA_loc(2);
    ApEn_loc = HRV.ApEn(nNA_RR_loc);
    [fv_pLF,fv_pHF,fv_LFHFratio,fv_VLF,fv_LF,fv_HF] = HRV.fft_val(nNA_RR_loc,60,Fs,'linear');
    [fvf_pLF,fvf_pHF,fvf_LFHFratio, fvf_VLF, fvf_LF, fvf_HF, f,Y,NFFT] = HRV.fft_val_fun(nNA_RR_loc,Fs,'linear');
    [SD1,SD2,SD1SD2ratio] = HRV.returnmap_val(RR,0,1,0);

    % Create output array
    hrv_output_array = [hrSegDuration, HR_loc, rrHRV_loc, SDNN_loc, SDSD_loc, RMSSD_loc, pNN50_loc, TRI_val_loc, TINN_val_loc, DFA_1_loc, DFA_2_loc, ApEn_loc, fvf_pLF, fvf_pHF, fvf_LFHFratio, fvf_VLF, fvf_LF, fvf_HF, SD1, SD2, SD1SD2ratio];
else
    hrv_output_array = nan(1,21);
end

end