function [new_datetime] = convert_time_RP(old_time)
% This function converts time stamps to a standartized format
      new_datetime = datetime(double(old_time)/1e7,'ConvertFrom','epochtime','Epoch','1-Jan-0001','Format','dd-MMM-yyyy HH:mm:ss.SSSSSSSSS');
end
