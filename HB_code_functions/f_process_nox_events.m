function [respEventList, AHI_ann, ODI_ann, ...
    Error_RP_no_event_file, Error_RP_no_events] = f_process_nox_events(edfPath, filename, systemUsed, satSignal, satFreq, StartRec, EndRec)
% This function processed respiratory events 
% from the respiratory polygraphy (Noxturnal). 
% Event data is stored in a file called "DATA.ndb" 
% and can be parsed with "mksqlite".
% The code creates an array ("respEvent") with the respiratory data

%% Find event file
eventFile = dir(fullfile(edfPath, 'Data.ndb'));
if size(eventFile)==0 % go get MARS annotation info
    
    
end
Error_RP_no_event_file = [];
Error_RP_no_events = [];
respEventList = NaN;
AHI_ann = NaN;
ODI_ann = NaN;

% Check if it exists
if numel(eventFile) == 0
    respEventList = NaN;
    AHI_ann = NaN;
    ODI_ann = NaN;
    disp("no Events file found ")
    Error_RP_no_event_file = 1;
    Error_RP_no_events = NaN;
else
    Error_RP_no_event_file = 0;
    Error_RP_no_events = 0;

    %% Open database
    if systemUsed == "Mac"
        db = mksqlite('open', fullfile(edfPath, 'Data.ndb'));
    else
        db = mksqlite1('open', fullfile(edfPath, 'Data.ndb'));
    end

    % Extract AHI and ODI of the recording
    sqlQueryEvents = 'SELECT * FROM internal_property';
    if systemUsed == "Mac"
        rec_properties = mksqlite(sqlQueryEvents);
    else
        rec_properties = mksqlite1(sqlQueryEvents);
    end

    rec_properties = struct2table(rec_properties);

    AHI_ann = find(contains(rec_properties.key, {'AHI'}));
    if ~isempty(AHI_ann)
        AHI_ann = rec_properties.value(AHI_ann);
    else
        AHI_ann = NaN;
    end
    
    ODI_ann = find(contains(rec_properties.key, {'ODI'}));
    if ~isempty(ODI_ann)
        ODI_ann = rec_properties.value(ODI_ann);
    else
        ODI_ann = NaN;
    end

    %% Extract events
    % the table that we are interested in is called:
    % "temporary_scoring_marker", so let's get the data from there
    sqlQueryEvents = 'SELECT * FROM temporary_scoring_marker';
    % if you want to display all the available tables, use:
    %sqlQueryTables = 'SELECT * FROM sqlite_master WHERE type="table"';
    % read the relevant data into a cell array called "eventArray"
    if systemUsed == "Mac"
        try
            eventArray = mksqlite(sqlQueryEvents);
        catch
            disp("table not found in database");
            sqlQueryEvents = 'SELECT * FROM scoring_marker';
            eventArray = mksqlite(sqlQueryEvents);
        end
        mksqlite('close')
    else
        try
            eventArray = mksqlite1(sqlQueryEvents);
        catch
            disp("table not found in database");
            sqlQueryEvents = 'SELECT * FROM scoring_marker';
            eventArray = mksqlite1(sqlQueryEvents);
        end
        mksqlite1('close')
    end
    
    % we select "hypopnea, apnea-central, apnea-obstructive --> common
    % element is "pnea", so grep for that
    respEventList = eventArray(contains({eventArray.type}, "pnea" ),:);

    % If there are no events detected
    if isempty(respEventList);
        Error_RP_no_events = 1;
    elseif class(respEventList) == 'double';
        Error_RP_no_events = 1;
    else
        % I assume that the "is_deleted" variable is "0" for not deleted and "1" is deleted
        % so let's remove all 1's
        respEventList = respEventList([respEventList.is_deleted] == 0, :);

        % Convert the absolute time of the trace to the time of its beginning
        if ~isempty(respEventList)

            % The time of the analysis of the events is limited
            % by the time of the saturation trace
            start_point = StartRec;
            end_point = EndRec;

            % Convert the format of the time stamps of the events
            if sum(size(respEventList)) > 2
                respEventList = struct2table(respEventList);
                for ll = 1:height(respEventList)
                    time_stamp = respEventList.starts_at(ll);
                    respEventList.starts_conv(ll) = f_convert_time_RP(time_stamp);
                    time_stamp = respEventList.ends_at(ll);
                    respEventList.ends_conv(ll) = f_convert_time_RP(time_stamp);
                end
            else
                ll = 1;
                resp_tab = table(respEventList.id);
                resp_tab.starts_at = respEventList.starts_at;
                resp_tab.ends_at = respEventList.ends_at;
                resp_tab.type = respEventList.type;
                resp_tab.location = respEventList.location;
                resp_tab.delete = 0;
                resp_tab.is_deleted = 0;
                resp_tab.key_id = respEventList.key_id;
                time_stamp = resp_tab.starts_at(ll);
                resp_tab.starts_conv(ll) = f_convert_time_RP(time_stamp);
                time_stamp = resp_tab.ends_at(ll);
                resp_tab.ends_conv(ll) = f_convert_time_RP(time_stamp);
                respEventList = resp_tab;
            end

            % Keep respiratory events between analysis start and end points
            toDelete = respEventList.starts_conv < start_point | respEventList.ends_conv > end_point;
            respEventList(toDelete,:) = [];

            % Create columns indicating event start and end
            % with regards to the time stamp of the analysis start
            for ll = 1:height(respEventList)
                respEventList.starts_relative_sec(ll) = seconds(diff(datetime([start_point;respEventList.starts_conv(ll)])));
                respEventList.ends_relative_sec(ll) = seconds(diff(datetime([start_point;respEventList.ends_conv(ll)])));
            end

            % Select columns with events start and event end
            respEventList = respEventList(:,["starts_relative_sec","ends_relative_sec"]);

        else
            % this should catch the situation 
            % where no events have been recorded
            respEventList = [];
            Error_RP_no_events = 1;
            disp("no RP Events detected")
        end
    end
end
end
