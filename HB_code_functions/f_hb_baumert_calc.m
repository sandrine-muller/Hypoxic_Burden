function [areaBSum, durationBMean, durationBSum, ...
    areaBStartSum, durationBStartMean, durationBStartSum, ...
    areaB3pSum, durationB3pMean, durationB3pSum] = f_hb_baumert_calc(satSignal, satFreq, ...
    desatStartOK, resatEnd, desatResatOK, Baumert_cutoff_level, Baumert_o2Threshold)
% This function is used to calculate the hypoxic burden
% according to Baumert method

% To get %*min/hr as a unit, we set the step-size in "trapz" to 1/(satFreq*60)
% For each of the calculated parameters, we substract the cutoff level 
% (90%, start, 3% drop) from the signal of the event, and calculate everything that is
% below 0 (taking the absolute of that value)
% Finally, we divide the value by the total hours of the recording

satSignal_nonnan = satSignal(~isnan(satSignal));

duration_rec = (length(satSignal_nonnan)/satFreq)/(60*60);
resatNumber = length(resatEnd);

% 1) Hypoxic Burden Baumert
areaB = zeros(1, resatNumber-1);
durationB = zeros(1, resatNumber-1);
for m = 1:resatNumber
    satSignTmp = satSignal(desatStartOK(m):resatEnd(m));
    if anynan(satSignTmp) == 0
        satSignTmp(satSignTmp > Baumert_cutoff_level) = Baumert_cutoff_level; % substract the cut-off, anything below the curve is our region of interest
        satSignTmp = satSignTmp-Baumert_cutoff_level;
        % Check the flat signal
        if range(satSignTmp) == 0
            durationB(m) = 0;
            areaB(m) = 0;
        else
            durationB(m) = length(satSignTmp)/satFreq;
            areaB(m) = abs(trapz(1/(satFreq*60),-satSignTmp));
        end
    else
        durationB(m) = NaN;
        areaB(m) = NaN;
    end
end
areaBSum = nansum(areaB)/duration_rec; 
durationBMean = nanmean(durationB); % in seconds
durationBSum = (nansum(durationB)/60/60)/duration_rec*100; % in % from the time

% 2) Hypoxic Burden B-start
% This an alternative version is the Baumert algorithm
% Event detection is the same, but we measure the area 
% under the curve compared to the starting level of the event
areaBStart = zeros(1, resatNumber-1);
durationBStart = zeros(1, resatNumber-1);
for m = 1:resatNumber
    satSignTmp = satSignal(desatStartOK(m):resatEnd(m));
    if anynan(satSignTmp) == 0
        sig_tresh = satSignal(desatStartOK(m));
        satSignTmp(satSignTmp > sig_tresh) = sig_tresh; % substract the starting saturation of this event
        satSignTmp = satSignTmp - sig_tresh;
        if range(satSignTmp) == 0
            durationBStart(m) = 0;
            areaBStart(m) = 0;
        else
            durationBStart(m) = length(satSignTmp)/satFreq;
            areaBStart(m) = abs(trapz(1/(satFreq*60),satSignTmp(satSignTmp<=0)));
        end
    else
        durationBStart(m) = NaN;
        areaBStart(m) = NaN;
    end
end
areaBStartSum = nansum(areaBStart)/duration_rec;
durationBStartMean = nanmean(durationBStart); % in seconds
durationBStartSum = (nansum(durationBStart)/60/60)/duration_rec*100; % in % from the time

% 3) Hypoxic Burden Baumert_o2Threshold
% This follows the logic from Watanabe et al.:
% only consider the desaturations which are greater than 3%
areaB3p = zeros(1, resatNumber-1);
durationB3p = zeros(1, resatNumber-1);
for m = 1:resatNumber
    satSignTmp = satSignal(desatStartOK(m):resatEnd(m));
    if anynan(satSignTmp) == 0
        if range(satSignTmp) > Baumert_o2Threshold
            satSignTmp = satSignTmp - max(satSignTmp); %substract the starting saturation - 3.0% of this event
            durationB3p(m) = length(satSignTmp)/satFreq;
            areaB3p(m) = abs(trapz(1/(satFreq*60),-satSignTmp));
        else
            durationB3p(m) = NaN;
            areaB3p(m) = NaN;
        end
    else
        durationB3p(m) = NaN;
        areaB3p(m) = NaN;
    end
end
areaB3pSum = nansum(areaB3p)/duration_rec;
durationB3pMean = nanmean(durationB3p); % in seconds
durationB3pSum = (nansum(durationB3p)/60/60)/duration_rec*100; % in % from the time

end