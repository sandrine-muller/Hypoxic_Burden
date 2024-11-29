function [hdr, signal, satFreq, satChannel, satSignal, hr_yes, hrFreq, ...
    hrChannel, hr_type, hrSignal, Error_satSignal_missing] = f_process_edf(filename)
% This function extract main properties and signals from the edf.

%filename = fullFileNames{k}
hdr = [];
signal = [];
satFreq = [];
satChannel = [];
satSignal = [];
hr_yes = [];
hrFreq = [];
hrChannel = [];
hr_type = [];
hrSignal = [];
Error_satSignal_missing = [];

% Reads the edf file and returns the (prefiltered) saturation signal
[signal, hdr] = read_edf(filename);
% [signal, hdr] = edf_read(filename);

% Check for missing satSignal
% Get info from the header
%maxFreq = max(hdr.Fs);
% Saturation
satChannel = find(contains(hdr.name, {'Saettigung' 'Entsaettigung' 'SAT' 'Saturation' 'SaO2' 'SpO2' 'SpO2BB'}),1);

% Check for missing satSignal
if ~isempty(satChannel)
    % If the signal was found
    Error_satSignal_missing = 0;
    satChannelname = hdr.name(satChannel);
    satFreq = hdr.Fs(satChannel);

    % Read the saturation signal
    satSignal = signal.(satChannelname{1,1});

else
    Error_satSignal_missing = 1;
end

% Heart rate
hrChannel_1 = find(contains(hdr.name, {'Pulse Waveform' 'PulseWaveform'}));
hrChannel_2 = find(contains(hdr.name, {'Pleth'}));
hrChannel_3 = find(contains(hdr.name, {'EKG' 'ECG'}));

% Select correct HR channel
if length(hrChannel_1) > 0
    hrChannel = hrChannel_1;
elseif length(hrChannel_2) > 0
    hrChannel = hrChannel_2;
elseif length(hrChannel_3) > 0
    hrChannel = hrChannel_3;
else
    hrChannel = 'NoHR';
end

% If no channel, replace with NA
% Otherwise, extract values
if hrChannel == 'NoHR';
    hr_yes = 0;
    hrChannelname = NaN;
    hrFreq = NaN;
    hrSignal = NaN;
    hr_type = NaN;
else
    hr_yes = 1;
    hrChannelname = hdr.name(hrChannel);
    hrFreq = hdr.Fs(hrChannel);
    hrChannelname = strrep(hrChannelname,' ','');
    hrSignal = signal.(hrChannelname{1,1});
    % Find the type of the channel
    if contains(hrChannelname, {'EKG' 'ECG'})
        hr_type = 1; % for ECG
    else
        hr_type = 2; % for PWV 
    end
end

end


