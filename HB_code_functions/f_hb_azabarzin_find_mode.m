function [mode_sat] = f_hb_azabarzin_find_mode(satSignTmp, AZB_mode_low, AZB_mode_high)
% This functions finds mode only in a physiological range defined 
% by AZB_mode_low and AZB_mode_high

mode_sat = [];
if range(satSignTmp) > 0 & length(satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp)) > 1;
    satSignTmp = satSignTmp(satSignTmp > AZB_mode_low & AZB_mode_high > satSignTmp);
    mode_sat = mode(satSignTmp(~isnan(satSignTmp)));
else
    mode_sat = nanmean(satSignTmp);
end

end