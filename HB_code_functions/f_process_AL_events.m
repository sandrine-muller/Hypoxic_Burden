function [respEventList, Error_AL_no_event_channels] = f_process_AL_events(hdr, signal)
% This function extracts respiratory events from an Apnea Link file.
% The events are stored in the .EDF files directly, as a binary-coded channel.
% This function finds and imports respiratory event signals
% that are coded in different columns.

Error_AL_no_event_channels = [];
% Here we have to add a checkpoint to ensure we dispose from events
eventTypes = find(contains(hdr.name, {'Hypopnoe' 'ObstruktiveApno' 'ZentraleApnoe' 'GemischteApnoe' 'Unklassifizierte'}),1);
if  numel(eventTypes) == 0
    disp("no Events channels found in AL edf")
    respEventList = NaN;
    Error_AL_no_event_channels = 1;
else
    Error_AL_no_event_channels = 0;
    
    hypChannel = find(contains(hdr.name, {'Hypopnoe'}),1);
    oApChannel = find(contains(hdr.name, {'Obstruktive Apno'}),1);
    cApChannel = find(contains(hdr.name, {'Zentrale Apnoe'}),1);
    mApChannel = find(contains(hdr.name, {'Gemischte Apnoe'}),1);
    uApChannel = find(contains(hdr.name, {'Unklassifizierte'}),1);

    % Get the signal of the channels
    hypChannelname = hdr.name(hypChannel);
    hypChannelFreq = hdr.Fs(hypChannel);
    hypSignal = signal.(hypChannelname{1,1});

    oApChannelname = hdr.name(oApChannel);
    oApChannelname = strrep(oApChannelname, ' ', '');
    oApChannelFreq = hdr.Fs(oApChannel);
    oApSignal = signal.(oApChannelname{1,1});

    cApChannelname = hdr.name(cApChannel);
    cApChannelname = strrep(cApChannelname, ' ', '');
    cApChannelFreq = hdr.Fs(cApChannel);
    cApSignal = signal.(cApChannelname{1,1});

    mApChannelname = hdr.name(mApChannel);
    mApChannelname = strrep(mApChannelname, ' ', '');
    mApChannelFreq = hdr.Fs(mApChannel);
    mApSignal = signal.(mApChannelname{1,1});

    uApChannelname = hdr.name(uApChannel);
    uApChannelname = strrep(uApChannelname, ' ', '');
    uApChannelFreq = hdr.Fs(uApChannel);
    uApSignal = signal.(uApChannelname{1,1});

    % Create a table with the events
    % Merge all events together
    respSignal = (hypSignal+oApSignal+cApSignal+mApSignal+uApSignal);

    % Check if respSignal is empty
    if sum(respSignal)>0
        % If not empty
        % Create matrixes with relative indexes
        original_event_index = find(respSignal == 1);
        find_borders_of_original_segments = diff(original_event_index);
        index_of_the_borders = find(find_borders_of_original_segments>1);

        % Run the loop to retrieve the original index
        all_pairs = [];
        for i = 1:length(index_of_the_borders)
            if i == 1
                seg_start = original_event_index(1);
                seg_end = index_of_the_borders(i);
                seg_end = original_event_index(seg_end);
                %plot(respSignal(seg_start:seg_end))
                %plot(respSignal(seg_start-50:seg_end+50))
            elseif i == length(index_of_the_borders)
                seg_start = index_of_the_borders(i-1)+1;
                seg_start = original_event_index(seg_start);
                seg_end = original_event_index(length(original_event_index));
            else
                seg_start = index_of_the_borders(i-1)+1;
                seg_start = original_event_index(seg_start);
                seg_end = index_of_the_borders(i);
                seg_end = original_event_index(seg_end);
            end
            add_pair = [seg_start seg_end];
            all_pairs = [all_pairs; add_pair];
        end

        % Convert the values from points to seconds
        all_pairs = all_pairs/hypChannelFreq;

        % Create a table with the respEventList
        respEventList = array2table(all_pairs);
        respEventList.Properties.VariableNames = ["starts_relative_sec" "ends_relative_sec"];

    else % if the signal is empty
        respEventList = [];
        disp("no AL Events detected in file")
    end
end
end