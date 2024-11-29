
% Version 20230906 IF, MD, SB
% System: MacOS, Matlab R2022b

%%%% This script is intended to calculate hypoxic burden from resp. polygraphy edf files (Apnealink and Noxturnal).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a path and preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

% Which groups of parameters should be calculated?
simple_AUC_yes = 1;
Baumert_yes = 1;
AZB_yes = 1;
FI_power_yes = 1;
FI_peaks_yes = 1;
FI_hrv_peaks_yes = 1;
FI_hrv_respEv_yes = 1;

% Add preferences and dependancies
filePath = 'C:\Users\admin-mullersa\Documents\data\test_HB\PA16953'; % Folder where you store the files
code_folder = 'C:\Users\admin-mullersa\Documents\code\Hypoxic_Burden_Bern\Current_231201\Current_231201\'; % Folder where you store the code
systemUsed = "Windows"; % "Windows" or "Mac"
path_save_output_matrix = 'C:\Users\admin-mullersa\Documents\data\test_HB\out\'; % Folder where you want to save the matrix
save_raw_output = 0; % 1 - yes, 0 - no (recommended)
path_save_raw_output = 'C:\Users\admin-mullersa\Documents\data\test_HB\out\save_raw'; % Folder where you want to save raw output

% Automatically add toolboxes
pathMksqlite = [];
if systemUsed == "Mac"
    pathMksqlite = strcat(code_folder, 'mksqlite-master');
else
    pathMksqlite = strcat(code_folder, 'mksqlite-2.11-win64');
end

pathFunctions = strcat(code_folder, 'HB_code_functions');
path_hrv_toolbox = strcat(code_folder, 'HRV-master');
if systemUsed == "Mac"
    hr_qrs_loc = strcat(path_hrv_toolbox, filesep, 'qrs_settings.mat'); % Location of the qrs file for the HRV
else
    hr_qrs_loc = strcat(path_hrv_toolbox, filesep, 'qrs_settings.mat'); % Location of the qrs file for the HRV
end

% Define parameters for the removal of signal outliers
% If the signal includes less than n_sat_val unique values,
% it is likely to be defect
n_sat_val = 4; % values
% If the signal stays the same for more than ssec
% than it is likely to be an artifact and it should be removed
ssec = 180; % seconds
% If the segment changes more than 4% (change_range) within 1 second (change_time)
% replace with NaN
change_time = 1; % seconds
change_range = 4; % range in absolute % saturation
% The signal should not change faster than filt_cutoff_diff
filt_cutoff_diff = 5; % higheest
filt_low_range = 40; % lowest possible value of the signal in absolute % saturation
filt_high_range = 100; % highest possible value of the signal in absolute % saturation
nan_min_dur = 5; % minimal duration of NaN in seconds to replace
% If 85% of the lonely segment is under XX% saturation - remove it
% XX% is signal mean - segment_tresh
segment_time = 85; % time in %
segment_tresh = 7;  % saturation in %
% segment_tresh = 80; % saturation, now it is defined in the code as the mean - 10
segments_keep = 120; % minimal duration of the segment to keep in seconds
segments_all_low = 90; % if all the values of the isolated segments are there - remove
% minimal time of the clean signal to check if it is enought to continue
% the analysis
time_clean = 20; % in min

% Define parameters for respEvents
min_resp_event_dur = 10; % seconds

% Define parameters for the detection of saturations
o2Threshold = 3.0; % Threshold set to > 3% drop
desat_duration = 50; % Maximal search window in seconds for desaturation in seconds
% Associated resaturation happens within 100 s and reached at least 2/3 of the drop level
resat_duration = 100; % Maximal search window in seconds for resaturation in seconds (Azarbarzin uses 100s)
resat_drop_tresh = 2/3; % Drio level

% Define parameters for the calculation of HB
simple_AUC_level = 90;

% Baumert
Baumert_cutoff_level = 90; % Saturation line in %
Baumert_o2Threshold = 3.0; % Threshold set to > 3% drop

% AZB
AZB_window_find_max_beforeRE_ends = 100; % time window to search for treshold values before the end of the event
AZB_minimal_distance_between_events = 20; % minimal time between respiratory events in seconds
AZB_window_beforeRE = 2; % time window in seconds to calculate AUC before the respiratory event start
AZB_window_afterRE = 60; % time in seconds to calculate AUC after the respiratory event start
AZB_mode_low = 91; % % of absolute saturation - minimal possible values for mode treshold
AZB_mode_high = 95; % % of absolute saturation - maximal possible values for mode treshold
AZB_pp = 30; % seach time window in seconds with regards to the event duration

% FI
FI_minimal_distance_between_events = 20; % minimal time between respiratory events in seconds
FI_peak_prominence = 3; % drop at least 3% of absolute saturation on either side before the signal attains a higher value
FI_peak_distance = 20; % minimal time between the peaks in seconds
FI_peak_event_before = 2; % time in seconds for the location of the associated peak before the start of the the event
FI_peak_event_after = 20; % time in seconds for the location of the associated peak after the start of the the event
FI_time_wakesat = 10; % time in seconds to detect wake saturation from the start of the recording
FI_treshold_deviation = 4; % potential "normal"/"physiological" deviation of the saturation for "wake" values
FI_TT = 100; % time window in seconds find the closest treshold_mid value to the peak
FI_hrv_peaks_before = FI_peak_distance; % HRV window before peak
FI_hrv_peaks_after = FI_peak_distance; % HRV window after peak
FI_hrv_respEvent_before = FI_minimal_distance_between_events; % HRV window before respEvent
FI_hrv_respEvent_after = FI_minimal_distance_between_events; % HRV window after respEvent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path
restoredefaultpath
addpath(genpath(pathMksqlite))
addpath(genpath(filePath))
addpath(genpath(pathFunctions))
addpath(genpath(path_hrv_toolbox))

% Define string separator
if systemUsed == "Mac"
    ssep = '/';
else
    ssep = '\';
end

%% Step A. Create an array of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***
satContents_fileinfo = {'Number', 'Filepath', 'RecordID', 'Visit', 'Device'};
% ***
satContents_error = {'Error_with_record', 'Error_satSignal_missing', ...
    'Error_satSignal_defect', 'Error_no_clean_satSignal', 'Error_extract_respEvents', ...
    'Error_RP_no_event_file', 'Error_RP_no_events', 'Error_AL_no_event_channels', ...
    'Error_no_respEvents', 'Error_process_respEvents', 'respEvents_number', 'Error_extract_sat', ...
    'Error_basic_sat', 'Error_basic_sat_AUC_level', 'Error_HB_BMT', 'Error_AZB_respEvents', 'Error_AZB_val', 'Error_undefined'};
% ***
satContents_timeinfo = {'satFreq', 'hr_yes', 'hr_type', 'hrFreq', ...
    'satSignal_duration', 'satSignal_nonnan_duration', 'arti_percent', ...
    'StartRec', 'EndRec', 'Recording_duration_ann', ...
    'AHI_ann','ODI_ann', 'MeanRespEventDur', 'AHI_calc_satSignal'};
% ***
satContents_basic_parameters = {'Mean_Sat', 'Lowest_Sat', ...
    'Mean_Drop', 'DesatEventLength', ...
    'T90_percent', 'T90_min', 'ODI_calc_satSignal', 'AUC_nonspecific_drift'};
% ***
satContents_AUC_level = {'AUC_under_level'};
% ***
satContents_hb_Baumert = {'areaBSum', 'durationBMean', 'durationBSum', ...
    'areaBStartSum', 'durationBStartMean', 'durationBStartSum', ...
    'areaB3pSum', 'durationB3pMean', 'durationB3pSum'};
% ***
satContents_hb_AZB = {'AZB_AHI', ...
    'AZB_max_area_sum', 'AZB_max_duration_mean', 'AZB_max_duration_sum', ...
    'AZB_mode_area_sum', 'AZB_mode_duration_mean', 'AZB_mode_duration_sum'};
% ***
satContents_hb_FI_power = {'FI_all_bandpower', 'FI_wd', 'FI_powtot', 'FI_pxx'};
% ***
satContents_hb_FI_peaks = {'FI_peaks_yes', 'FI_peak_valid', 'FI_peak_index', 'FI_mean_peak_width', 'FI_mean_peak_height', ...
    'FI_peak_area_sum', 'FI_peak_duration_mean', 'FI_peak_duration_sum', ...
    'FI_AHI', 'FI_peak_respEvent_valid', 'FI_peaks_respEvent_yes', 'FI_peak_respEvent_index', ...
    'FI_mean_peak_respEvent_width', ' FI_mean_peak_respEvent_height', ...
    'FI_peak_respEvent_area_sum', 'FI_peak_respEvent_duration_mean', ...
    'FI_peak_respEvent_duration_sum'};
% ***
satContents_hb_FI_hrv_peaks = {'peak_hrSegDuration', 'peak_HR_loc', 'peak_rrHRV_loc', ...
    'peak_SDNN_loc', 'peak_SDSD_loc', 'peak_RMSSD_loc', ...
    'peak_pNN50_loc', 'peak_TRI_val_loc', 'peak_TINN_val_loc', ...
    'peak_DFA_1_loc', 'peak_DFA_2_loc', 'peak_ApEn_loc', ...
    'peak_fvf_pLF', 'peak_fvf_pHF', 'peak_fvf_LFHFratio', ...
    'peak_fvf_VLF', 'peak_fvf_LF', 'peak_fvf_HF', ...
    'peak_SD1', 'peak_SD2', 'peak_SD1SD2ratio'};
% ***
satContents_hb_FI_hrv_respEv = {'respEv_hrSegDuration', 'respEv_HR_loc', 'respEv_rrHRV_loc', ...
    'respEv_SDNN_loc', 'respEv_SDSD_loc', 'respEv_RMSSD_loc', ...
    'respEv_pNN50_loc', 'respEv_TRI_val_loc', 'respEv_TINN_val_loc', ...
    'respEv_DFA_1_loc', 'respEv_DFA_2_loc', 'respEv_ApEn_loc', ...
    'respEv_fvf_pLF', 'respEv_fvf_pHF', 'respEv_fvf_LFHFratio', ...
    'respEv_fvf_VLF', 'respEv_fvf_LF', 'respEv_fvf_HF', ...
    'respEv_SD1', 'respEv_SD2', 'respEv_SD1SD2ratio'};
% *** Sum up
satContents = [satContents_fileinfo satContents_error ...
    satContents_timeinfo ...
    satContents_basic_parameters];
% Add selected parameters
if simple_AUC_yes == 1;
    satContents = [satContents satContents_AUC_level];
else
    %
end
if Baumert_yes == 1;
    satContents = [satContents satContents_hb_Baumert];
else
    %
end
if AZB_yes == 1;
    satContents = [satContents satContents_hb_AZB];
else
    %
end
if FI_power_yes == 1;
    satContents = [satContents satContents_hb_FI_power];
else
    %
end
if FI_peaks_yes == 1;
    satContents = [satContents satContents_hb_FI_peaks];
else
    %
end
if FI_hrv_peaks_yes == 1;
    satContents = [satContents satContents_hb_FI_hrv_peaks];
else
    %
end
if FI_hrv_respEv_yes == 1;
    satContents = [satContents satContents_hb_FI_hrv_respEv];
else
    %
end

%% Step B. Create a list of files to process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a list of all edf files in the current folder, or subfolders of it.
fds = fileDatastore(filePath, 'ReadFcn', @edfread, 'FileExtensions',{'.edf'}, 'IncludeSubfolders', true);
fullFileNames = fds.Files;

% Discard .edfplus filenames. They are consistently named '_edfplus.EDF'
str = find(contains(fullFileNames, {'edfplus' 'edflplus' 'edfpflus'}));
fullFileNames(str,:)=[];

% Number the files
numFiles = length(fullFileNames);

%% Step C. Create a matrix to records values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
satData = cell(numFiles, length(satContents));

% Create raw_output structure if needed
if save_raw_output == 1
    rawContents = {'Number', 'Filepath', ...
        'RecordID', 'Visit', ...
        'satChannel', 'satSignal', ...
        'hrChannel', 'hrSignal', ...
        'respEventList'};
    rawData = cell(length(rawContents),1);
    rawData = cell2struct(rawData,rawContents);
else
    %
end

%% Step D. Save the workspace
% The other variables will be deleted at the end of the file loop
% Open a progress bar
waitHandle = waitbar(0,'Please wait...', 'Name', 'Calculating hypoxic burden');
waitbar_marker = [];
k = [];
keepvars = who;
keepvars = [keepvars; 'keepvars'];
keepvars = [keepvars; 'k'];

% Code check
% load('numbers_to_follow.mat')
% numbers_to_follow = numbers_to_follow(~isnan(numbers_to_follow));

%% Step E. Process the files (loop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : numFiles
    % for g = 1 : length(numbers_to_follow)
    %     k = numbers_to_follow(g);
    % end for k=1
    try
        % Update progress bar
        waitbar_marker = k;
        waitbar(waitbar_marker/numFiles, waitHandle, sprintf('Recording %d of %d', waitbar_marker, numFiles))
        fprintf('Now reading edf file %s\n', fullFileNames{k});

        % Define baseline errors
        Error_with_record = [];
        Error_satSignal_missing = [];
        Error_satSignal_defect = [];
        Error_no_clean_satSignal = [];
        Error_extract_respEvents = [];
        Error_RP_no_event_file = [];
        Error_RP_no_events = [];
        Error_AL_no_event_channels = [];
        Error_no_respEvents = [];
        Error_process_respEvents = [];
        respEvents_number = [];
        Error_extract_sat = [];
        Error_basic_sat = [];
        Error_HB_BMT = [];
        Error_AZB_respEvents = [];
        Error_AZB_val = [];
        Error_undefined = [];

        %% Step E1. Retrieve the file information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract filename and save values in the array
        inds = find(contains(satContents, 'Number'));
        satData(k,inds) = {k}; % number of the loop
        inds = find(contains(satContents, 'Filepath'));
        satData(k,inds) = {fullFileNames(k)}; % path of the file (including the filename!!)
        [edfPath, edfName, ext] = fileparts(fullFileNames{k});

        % Extract RecordID based on the name of the directoty
        [record_id, Error_with_record] = f_extract_record(edfName, edfPath);
        inds = find(contains(satContents, 'RecordID'));
        satData(k,inds) = record_id;
        inds = find(contains(satContents, 'Error_with_record'));
        satData(k,inds) = {Error_with_record};

        % Extract Visit
        [visit] = f_extract_visit(edfName, edfPath);
        inds = find(contains(satContents, 'Visit'));
        satData(k,inds) = {visit};

        %% Step E2. Read the EDFs and extract signals and their properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename = fullFileNames{k};
        [hdr, signal, satFreq, satChannel, satSignal, ...
            hr_yes, hrFreq, hrChannel, hr_type, hrSignal, ...
            Error_satSignal_missing] = f_process_edf(filename);

        % NB: In SDSO dataset, there is one recording where the ECG
        % channel is wrongly indicated. Therefore, it should be overwritten
        obs_id = strcat(record_id{1,1}, visit);
        all_obs = {'C01S197V1', 'C01S393V1'};
        if any(contains(all_obs, obs_id));
            hr_yes = 0;
            hrChannel = 'NoHR';
            hrChannelname = NaN;
            hrFreq = NaN;
            hrSignal = NaN;
            hr_type = NaN;
        else
            %
        end

        % Write to array
        inds = find(contains(satContents, 'Error_satSignal_missing'));
        satData(k,inds) = {Error_satSignal_missing};
        inds = find(contains(satContents, 'satFreq'));
        satData(k,inds) = {satFreq};
        inds = find(contains(satContents, 'hr_yes'));
        satData(k,inds) = {hr_yes};
        inds = find(contains(satContents, 'hrFreq'));
        satData(k,inds) = {hrFreq};
        inds = find(contains(satContents, 'hr_type'));
        satData(k,inds) = {hr_type};

        %% Step E3. Find the start and the end of the recording
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        info = edfinfo(filename);
        date_rec = strsplit(info.StartDate, '.');
        date_rec = strcat('20', date_rec{3}, '.', date_rec{2}, '.', date_rec{1});
        StartRec = convertCharsToStrings(date_rec) + ' ' + strrep(info.StartTime,'.',':') + '.0';
        StartRec = datetime(StartRec, 'InputFormat', 'yyyy.MM.dd HH:mm:ss.S', 'Format', 'dd MMM yyyy HH:mm:ss');
        if Error_satSignal_missing == 0;
            EndRec = seconds(length(satSignal)/satFreq);
            EndRec = StartRec + EndRec;
            inds = find(contains(satContents, 'Recording_duration_ann'));
            Recording_duration_ann = minutes(EndRec - StartRec);
            satData(k,inds) = {Recording_duration_ann}; % min.
        elseif hr_yes == 1;
            EndRec = seconds(length(hrSignal)/hrFreq);
            EndRec = StartRec + EndRec;
            inds = find(contains(satContents, 'Recording_duration_ann'));
            Recording_duration_ann = minutes(EndRec - StartRec);
            satData(k,inds) = {Recording_duration_ann}; % min.
        else
            EndRec = [];
        end
        inds = find(contains(satContents, 'StartRec'));
        satData(k,inds) = {StartRec};
        inds = find(contains(satContents, 'EndRec'));
        satData(k,inds) = {EndRec};
        inds = find(contains(satContents, 'satFreq'));
        satData(k,inds) = {satFreq};

        %% Step E4. Read and process respiratory events
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This part does not run if there is a saturation signal missing
        disp("Extracting respiratory events...")
        device_inds = find(contains(satContents, 'Device'));
        % 1) Determine the file type and extract the respiratory events from the files

        % NB: There are records with a defect event file
        % In this case do not process the respiratory events and
        % assume that there are no events
        obs_id = strcat(record_id{1,1}, visit);
        all_obs = {'C01S348V1', 'C01S325V1', 'C01S354V1', ...
            'C01S377V1', 'C01S379V1', 'C01S403V1', ...
            'C01S418V1', 'C01S420V1', 'C01S434V1', ...
            'C01S441V1', 'C01S453V1', 'C01S398V1'};
        if any(contains(all_obs, obs_id));
            respEventList = [];
            Error_extract_respEvents = 1;
        else

            if Error_satSignal_missing == 0;
                % For RP
                if sum(contains(hdr.name, {'SAT' 'Saturation' 'SpO2BB' 'SpO2 B-B' 'SaO2' 'SpO2'})) > 0 % check if RP % SM added saturation signal name
                    satData(k,device_inds) = {'RP'};
                    % process the DATA.ndb --> mksqlite function
                    [respEventList, AHI_ann, ODI_ann, ...
                        Error_RP_no_event_file, Error_RP_no_events] = f_process_nox_events(edfPath, filename, systemUsed, satSignal, satFreq, StartRec, EndRec);
                    % Write the information to the array
                    inds = find(contains(satContents, 'AHI_ann'));
                    satData(k,inds) = {AHI_ann};
                    inds = find(contains(satContents, 'ODI_ann'));
                    satData(k,inds) = {ODI_ann};
                    inds = find(contains(satContents, 'Error_RP_no_event_file'));
                    satData(k,inds) = {Error_RP_no_event_file};
                    inds = find(contains(satContents, 'Error_RP_no_events'));
                    satData(k,inds) = {Error_RP_no_events};
                    % For AL
                elseif sum(contains(hdr.name, {'Saettigung' 'Entsaettigung'})) > 0 % check if AL.
                    satData(k,device_inds) = {'AL'};
                    % Process by reading the scored events in the edf file
                    [respEventList, Error_AL_no_event_channels] = f_process_AL_events(hdr, signal);
                    inds = find(contains(satContents, 'Error_AL_no_event_channels'));
                    satData(k,inds) = {Error_AL_no_event_channels};
                else
                    satData(k,inds) = {'unknown type'}; % shouldn't occur, but if it does --> check!
                end

                % Error if the extraction failed
                if exist('respEventList', 'var'),
                    Error_extract_respEvents = 0;
                else
                    Error_extract_respEvents = 1;
                end
                inds = find(contains(satContents, 'Error_extract_respEvents'));
                satData(k,inds) = {Error_extract_respEvents};

                % Error both for RP and AL:
                % If there are any events
                if Error_RP_no_event_file == 1
                    Error_no_respEvents = 1;
                elseif Error_AL_no_event_channels == 1
                    Error_no_respEvents = 1;
                elseif Error_RP_no_events == 1
                    Error_no_respEvents = 1;
                else
                    Error_no_respEvents = 0;
                end
                inds = find(contains(satContents, 'Error_no_respEvents'));
                satData(k,inds) = {Error_no_respEvents};

                % 2) Process respiratory events
                % Sort rows respEventList
                if Error_no_respEvents == 0
                    if class(respEventList) ~= "double"
                        if ~isempty(respEventList)
                            % Sort rows respEventList
                            respEventList = sortrows(respEventList);
                            % Check the minimal duration of events
                            respEventList.respEvDuration = respEventList.ends_relative_sec - respEventList.starts_relative_sec;
                            respEventList = respEventList(respEventList.respEvDuration > min_resp_event_dur, :);
                            % Delete events that are later than the saturation signal
                            to_delete = int64(respEventList.ends_relative_sec) > (length(satSignal)/satFreq);
                            respEventList = respEventList(~to_delete,:);
                            % Calculate properties resp. events and save to array
                            inds = find(contains(satContents, 'MeanRespEventDur'));
                            satData(k,inds) = {mean(respEventList.respEvDuration)};
                            % Save AHI calculated based on respEvents
                            % in relation to the satSignal duration
                            inds = find(contains(satContents, 'AHI_calc_satSignal'));
                            AHI_calc_satSignal = height(respEventList)/((length(satSignal)/satFreq)/60/60); % per hour
                            satData(k,inds) = {AHI_calc_satSignal};
                            % Save respEvent number
                            respEvents_number = height(respEventList);
                            inds = find(contains(satContents, 'respEvents_number'));
                            satData(k,inds) = {respEvents_number};
                            % Error
                            Error_process_respEvents = 0;
                        else
                            inds = find(contains(satContents, 'MeanRespEventDur'));
                            satData(k,inds) = {0};
                            inds = find(contains(satContents, 'AHI_calc_satSignal'));
                            satData(k,inds) = {0};
                            Error_process_respEvents = 0;
                            respEvents_number = 0;
                            inds = find(contains(satContents, 'respEvents_number'));
                            satData(k,inds) = {respEvents_number};
                        end
                    else
                        Error_process_respEvents = 1;
                    end
                else
                    %
                end

                inds = find(contains(satContents, 'Error_process_respEvents'));
                satData(k,inds) = {Error_process_respEvents}; % min.
            else
                %
            end

        end % Defect record

        %% Step E5. Signal analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Part E5-1: Saturation signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This runs only if the saturation signal is not missing
        if Error_satSignal_missing == 0;
            %% E5-1.1) Remove artifacts from the saturation signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp("Removing artfacts in a saturation signal ...")
            % satSignal_original = satSignal;
            [new_satSignal, Error_satSignal_defect] = f_process_satTrace(satSignal, satFreq, ...
                n_sat_val, ssec, change_time, change_range, ...
                filt_cutoff_diff, filt_low_range, filt_high_range, ...
                nan_min_dur, segment_tresh, segment_time, segments_keep, segments_all_low);

            index_error = find(contains(satContents, 'Error_satSignal_defect'));
            satData(k,index_error) = {Error_satSignal_defect}; % is not very meaningful, and should be changed manually

            satSignal = new_satSignal;
            satSignal_nonnan = satSignal(~isnan(satSignal));

            %% E5-1.2) Validity check
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check if there is anything left from the clean signal
            if isempty(satSignal_nonnan)
                Error_no_clean_satSignal = 1;
            elseif length(satSignal_nonnan)/satFreq < time_clean*60
                Error_no_clean_satSignal = 1;
                disp("No clean saturation signal ...")
            else
                Error_no_clean_satSignal = 0;
                disp("Saturation signal is cleaned ...")
            end
            index_error = find(contains(satContents, 'Error_no_clean_satSignal'));
            satData(k,index_error) = {Error_no_clean_satSignal}; % is not very meaningful, and should be changed manually

            if Error_no_clean_satSignal == 0

                % Write the information to the array
                inds = find(contains(satContents, 'satSignal_duration'));
                satSignal_duration = length(satSignal)/satFreq/60; % min
                satData(k,inds) = {satSignal_duration}; % min.
                satSignal_nonnan_duration = length(satSignal_nonnan)/satFreq/60; % min
                satSignal_nonnan_duration_h = satSignal_nonnan_duration/60; % hours
                arti_percent = 100-length(satSignal_nonnan)/length(satSignal)*100;
                inds = find(contains(satContents, 'arti_percent'));
                satData(k,inds) = {arti_percent};
                inds = find(contains(satContents, 'satSignal_nonnan_duration'));
                satData(k,inds) = {satSignal_nonnan_duration};

                %% E5-1.3) Calculate the area under the curve according to the defined treshold
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % To get %*min/hr as a unit, we set the step-size in "trapz" to 1/(satFreq*60)
                % for each of the variants, we substract the cutoff level (90%, start, 3%
                % drop ) from the signal of the event, and calculate everything that is
                % below 0 (taking the absolute of that value)
                % finally, we divide the value by the total hours of the recording.
                if simple_AUC_yes == 1;
                    disp("Calculating AUC (treshold) ...")
                    signal_sAUC = satSignal;
                    signal_sAUC(signal_sAUC > simple_AUC_level) = simple_AUC_level;
                    signal_sAUC = signal_sAUC - simple_AUC_level;
                    signal_sAUC = signal_sAUC(~isnan(signal_sAUC));
                    AUC_under_level = (trapz(1/(satFreq*60),-signal_sAUC))/satSignal_nonnan_duration_h; % NB this is the area, NOT the duration!!
                    inds = find(contains(satContents, 'AUC_under_level'));
                    satData(k,inds) = {AUC_under_level}; %*min/hr
                    Error_basic_sat_AUC_level = 0;
                    inds = find(contains(satContents, 'Error_basic_sat_AUC_level'));
                    satData(k,inds) = {Error_basic_sat_AUC_level}; % min.
                else
                    %
                end

                %% E5-1.4) Detect desaturations and resaturations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp("Extracting desaturation/resaturation events ...")
                % This function returns a vector with desaturations (all) and
                % resaturations. the desatResatOK logical vector indicates which of the detected desaturations
                % has an associated resaturation (of at least 2/3 of the drop level
                % within 100s)
                [desatEvent, resatEvent, desatResatOK] = f_detect_sat_events(satSignal_nonnan, satFreq, o2Threshold, desat_duration, resat_duration, resat_drop_tresh);

                % Create a single vector
                desatStart = find(diff(desatEvent) == 1);
                desatEnd = find(diff(desatEvent) == -1);
                resatEnd = find(diff(resatEvent) == -1);
                desatStartOK = desatStart(desatResatOK == 1);
                desatEndOK = desatEnd(desatResatOK == 1);
                desatDur = ((resatEnd - desatStartOK)/satFreq);

                % Convert to a single satEventvector
                % From this desat-resat vector we'll create one single resp-event
                % vector, coded -1 for desat (only associated events!!), 1 for resat
                % and 0 for the rest
                satEvent = resatEvent;
                for  j = 1:length(desatResatOK)
                    if desatResatOK(j) == 1
                        satEvent(desatStart(j):desatEnd(j)) = -1 ;
                    end
                end
                Error_extract_sat = 0;
                inds = find(contains(satContents, 'Error_extract_sat'));
                satData(k,inds) = {Error_extract_sat}; % min.

                %% E5-1.5) Calculate basic saturation parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                disp("Calculating basic saturation parameters ...")
                % If breathing information exists (either from an "events" file that
                % we parsed at point 2., or which was already included in the edf
                % (ie. with Apnealink files) --> also calculate acc. azarbarzin method)
                % Mean Saturation
                inds = find(contains(satContents, 'Mean_Sat'));
                satData(k,inds) = {mean(satSignal_nonnan)};
                % Lowest Saturation
                inds = find(contains(satContents, 'Lowest_Sat'));
                satData(k,inds) = {min(satSignal_nonnan)};
                % Mean desturation drop
                inds = find(contains(satContents, 'Mean_Drop'));
                satData(k,inds) = {nanmean(satSignal(desatStartOK) - satSignal(desatEndOK))};
                % Desaturation event duration in seconds
                inds = find(contains(satContents, 'DesatEventLength'));
                satData(k,inds) = {mean(desatDur)};
                % Percent time under 90%
                inds = find(contains(satContents, 'T90_percent'));
                satData(k,inds) = {length(satSignal_nonnan(satSignal_nonnan<=90))/length(satSignal_nonnan)*100};
                % Absolute time under 90% in minutes
                inds = find(contains(satContents, 'T90_min'));
                satData(k,inds) = {length(satSignal_nonnan(satSignal_nonnan<=90))/satFreq/60};
                % ODI
                resatNumber = length(desatResatOK(desatResatOK == 1));
                ODI_calc_satSignal = resatNumber/satSignal_nonnan_duration_h;
                inds = find(contains(satContents, 'ODI_calc_satSignal'));
                satData(k,inds) = {ODI_calc_satSignal};

                %% E5-1.6) Hypoxic burden acc. Baumert
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% E5-1.6.1) Calculate the area under the curve which is not associated with desurations/resaturations
                % (i.e., the unspecific drift in the Baumert paper)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % To get %*min/hr as a unit, we set the step-size in "trapz" to 1/(satFreq*60)
                % for each of the variants, we substract the cutoff level (90%, start, 3%
                % drop ) from the signal of the event, and calculate everything that is
                % below 0 (taking the absolute of that value)
                % finally, we divide the value by the total hours of the recording.
                if Baumert_yes == 1;
                    disp("Hypoxic burden acc. Baumert...")
                    disp("Computing ...")
                    satSign_nonevent = satSignal(satEvent == 0); % no desuration/resaturation
                    satSign_nonevent(satSign_nonevent > 90) = 90;
                    satSign_nonevent = satSign_nonevent - 90;
                    satSign_nonevent = satSign_nonevent(~isnan(satSign_nonevent));
                    AUC_nonspecific_drift = (trapz(1/(satFreq*60),-satSign_nonevent))/satSignal_nonnan_duration_h; % NB this is the area, NOT the duration!!
                    inds = find(contains(satContents, 'AUC_nonspecific_drift'));
                    satData(k,inds) = {AUC_nonspecific_drift}; %*min/hr
                    Error_basic_sat = 0;
                    inds = find(contains(satContents, 'Error_basic_sat'));
                    satData(k,inds) = {Error_basic_sat}; % min.

                    %% E5-1.6.2) Calculate the area under the curve which is associated with desurations/resaturations
                    % Here we can use different cutoffs to
                    % determine the hypoxic burden: 90% (like in the paper), starting
                    % saturation or end saturation. Here we collect all of them
                    [areaBSum, durationBMean, durationBSum, ...
                        areaBStartSum, durationBStartMean, durationBStartSum, ...
                        areaB3pSum, durationB3pMean, durationB3pSum] = f_hb_baumert_calc(satSignal, satFreq, ...
                        desatStartOK, resatEnd, desatResatOK, Baumert_cutoff_level, Baumert_o2Threshold);
                    % Write the information to the array
                    % %min./h
                    inds = find(contains(satContents, 'areaBSum'));
                    satData(k,inds) = {areaBSum};
                    % sec
                    inds = find(contains(satContents, 'durationBMean'));
                    satData(k,inds) = {durationBMean};
                    % %time
                    inds = find(contains(satContents, 'durationBSum'));
                    satData(k,inds) = {durationBSum};
                    % %min./h
                    inds = find(contains(satContents, 'areaBStartSum'));
                    satData(k,inds) = {areaBStartSum};
                    % sec
                    inds = find(contains(satContents, 'durationBStartMean'));
                    satData(k,inds) = {durationBStartMean};
                    % %time
                    inds = find(contains(satContents, 'durationBStartSum'));
                    satData(k,inds) = {durationBStartSum};
                    % %min./h
                    inds = find(contains(satContents, 'areaB3pSum'));
                    satData(k,inds) = {areaB3pSum};
                    % sec
                    inds = find(contains(satContents, 'durationB3pMean'));
                    satData(k,inds) = {durationB3pMean};
                    % %time
                    inds = find(contains(satContents, 'durationB3pSum'));
                    satData(k,inds) = {durationB3pSum};

                    Error_HB_BMT = 0;
                    inds = find(contains(satContents, 'Error_HB_BMT'));
                    satData(k,inds) = {Error_HB_BMT}; % min.
                else
                    %
                end % Baumert_yes == 1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% E5-1.7) Hypoxic burden acc. Azarbarzin
                if AZB_yes == 1;
                    disp("Hypoxic burden acc. Azarbarzin ...")
                    Error_AZB_respEvents = [];
                    % It can be calculated only if there are respiratory events
                    if ~isempty(respEventList) & class(respEventList) ~= "double"
                        disp("Computing ...")
                        [AZB_AHI, AZB_nevents_max, AZB_max_area_sum, ...
                            AZB_max_duration_mean, AZB_max_duration_sum, ...
                            AZB_nevents_mode, AZB_mode_area_sum, AZB_mode_duration_mean, ...
                            AZB_mode_duration_sum, Error_AZB_val] = f_hb_azarbarzin_calc(satSignal, satFreq, respEventList, ...
                            AZB_window_find_max_beforeRE_ends, ...
                            AZB_minimal_distance_between_events, AZB_window_beforeRE, ...
                            AZB_window_afterRE, AZB_mode_low, ...
                            AZB_mode_high, AZB_pp);
                        % Write to the array
                        % /h
                        inds = find(contains(satContents, 'AZB_AHI'));
                        satData(k,inds) = {AZB_AHI};
                        % %min./h
                        inds = find(contains(satContents, 'AZB_max_area_sum'));
                        satData(k,inds) = {AZB_max_area_sum};
                        % sec
                        inds = find(contains(satContents, 'AZB_max_duration_mean'));
                        satData(k,inds) = {AZB_max_duration_mean};
                        % %time
                        inds = find(contains(satContents, 'AZB_max_duration_sum'));
                        satData(k,inds) = {AZB_max_duration_sum};
                        % %min./h
                        inds = find(contains(satContents, 'AZB_mode_area_sum'));
                        satData(k,inds) = {AZB_mode_area_sum};
                        % sec
                        inds = find(contains(satContents, 'AZB_mode_duration_mean'));
                        satData(k,inds) = {AZB_mode_duration_mean};
                        % %time
                        inds = find(contains(satContents, 'AZB_mode_duration_sum'));
                        satData(k,inds) = {AZB_mode_duration_sum};
                        Error_AZB_respEvents = 0;
                        inds = find(contains(satContents, 'Error_AZB_val'));
                        satData(k,inds) = {Error_AZB_val};
                    else
                        disp("Cannot be computed ...")
                        Error_AZB_respEvents = 1;
                    end
                    inds = find(contains(satContents, 'Error_AZB_respEvents'));
                    satData(k,inds) = {Error_AZB_respEvents}; % min.
                else
                    %
                end % AZB_yes == 1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% E5-1.8) Hypoxic burden acc. Filchenko
                if FI_power_yes == 1;
                    disp("Hypoxic burden acc. Filchenko (power) ...")
                    disp("Computing ...")
                    %% E5-1.8.1) Calculate signal power
                    [FI_all_bandpower, FI_wd, FI_powtot, FI_pxx] = f_hb_fi_power_calc(satSignal, satFreq, respEventList);
                    % Write to the array
                    % %^2/Hz
                    inds = find(contains(satContents, 'FI_all_bandpower'));
                    satData(k,inds) = {FI_all_bandpower};
                    % Fs
                    inds = find(contains(satContents, 'FI_wd'));
                    satData(k,inds) = {FI_wd};
                    % %/Hz
                    inds = find(contains(satContents, 'FI_powtot'));
                    satData(k,inds) = {FI_powtot};
                    % %^2/Hz
                    inds = find(contains(satContents, 'FI_pxx'));
                    satData(k,inds) = {FI_pxx};
                else
                    %
                end % FI_power_yes == 1;

                %% E5-1.8.2) Find and describe peaks
                if FI_peaks_yes == 1;
                    disp("Hypoxic burden acc. Filchenko (peaks) ...")
                    [FI_peaks_yes, FI_peak_valid, FI_peak_index, FI_mean_peak_width, FI_mean_peak_height, ...
                        FI_peak_area_sum, FI_peak_duration_mean, FI_peak_duration_sum, ...
                        pks, locs, w, p] = f_hb_fi_oxy(satSignal, satFreq, ...
                        FI_peak_prominence, FI_peak_distance, ...
                        FI_time_wakesat, FI_TT, FI_treshold_deviation);
                    % Write to the array
                    % bin
                    inds = find(contains(satContents, 'FI_peaks_yes'));
                    satData(k,inds) = {FI_peaks_yes};
                    % number of valid peals
                    inds = find(contains(satContents, 'FI_peak_valid'));
                    satData(k,inds) = {FI_peak_valid};
                    % /h
                    inds = find(contains(satContents, 'FI_peak_index'));
                    satData(k,inds) = {FI_peak_index};
                    % seconds
                    inds = find(contains(satContents, 'FI_mean_peak_width'));
                    satData(k,inds) = {FI_mean_peak_width};
                    % %
                    inds = find(contains(satContents, 'FI_mean_peak_height'));
                    satData(k,inds) = {FI_mean_peak_height};
                    % %min./h
                    inds = find(contains(satContents, 'FI_peak_area_sum'));
                    satData(k,inds) = {FI_peak_area_sum};
                    % seconds
                    inds = find(contains(satContents, 'FI_peak_duration_mean'));
                    satData(k,inds) = {FI_peak_duration_mean};
                    % %time
                    inds = find(contains(satContents, 'FI_peak_duration_sum'));
                    satData(k,inds) = {FI_peak_duration_sum};
                    %% E5-1.8.3) Find and describe peaks, associated with respiratory events
                    % It can be calculated only if there are respiratory events
                    if ~isempty(respEventList) & class(respEventList) ~= "double"
                        [FI_AHI, FI_peaks_respEvent_yes, FI_peak_respEvent_valid, FI_peak_respEvent_index, ...
                            FI_mean_peak_respEvent_width, FI_mean_peak_respEvent_height, ...
                            FI_peak_respEvent_area_sum, FI_peak_respEvent_duration_mean, ...
                            FI_peak_respEvent_duration_sum] = f_hb_fi_oxy_respEvent(satSignal, satFreq, respEventList, ...
                            FI_minimal_distance_between_events, FI_peak_prominence, FI_peak_distance, ...
                            FI_peak_event_before, FI_peak_event_after, FI_time_wakesat, ...
                            FI_TT, FI_treshold_deviation, pks, locs, w, p);
                        % Write to the array
                        inds = find(contains(satContents, 'FI_AHI'));
                        satData(k,inds) = {FI_AHI};
                        % binary
                        inds = find(contains(satContents, 'FI_peaks_respEvent_yes'));
                        satData(k,inds) = {FI_peaks_respEvent_yes};
                        % number of valid peals
                        inds = find(contains(satContents, 'FI_peak_respEvent_valid'));
                        satData(k,inds) = {FI_peak_respEvent_valid};
                        % /h
                        inds = find(contains(satContents, 'FI_peak_respEvent_index'));
                        satData(k,inds) = {FI_peak_respEvent_index};
                        % seconds
                        inds = find(contains(satContents, 'FI_mean_peak_respEvent_width'));
                        satData(k,inds) = {FI_mean_peak_respEvent_width};
                        % %
                        inds = find(contains(satContents, 'FI_mean_peak_respEvent_height'));
                        satData(k,inds) = {FI_mean_peak_respEvent_height};
                        % %min./h
                        inds = find(contains(satContents, 'FI_peak_respEvent_area_sum'));
                        satData(k,inds) = {FI_peak_respEvent_area_sum};
                        % seconds
                        inds = find(contains(satContents, 'FI_peak_respEvent_duration_mean'));
                        satData(k,inds) = {FI_peak_respEvent_duration_mean};
                        % %time
                        inds = find(contains(satContents, 'FI_peak_respEvent_duration_sum'));
                        satData(k,inds) = {FI_peak_respEvent_duration_sum};
                    else
                        %
                    end % 8.3
                else
                    %
                end % FI_peaks_yes == 1;
            else
                % Do nothing since no clean satSignal
            end % if satSignal is clean
        else
            % Do nothing since no satSignal present
        end % if satSignal is present

        %% Part E5-2: Heart rate variability
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% E5-2.1: Calculate HRV around peaks
        if FI_hrv_peaks_yes == 1
            if exist('hr_yes') & exist('FI_peak_valid');
                disp("Hypoxic burden acc. Filchenko (HRV around peaks) ...")
                important_variables = [Error_satSignal_missing Error_no_clean_satSignal Error_satSignal_defect];
                if length(important_variables) == 3 & sum(important_variables) == 3;
                elseif hr_yes == 1 & FI_peak_valid > 1
                    disp("Computing ...")
                    [peak_hrSegDuration, peak_HR_loc, peak_rrHRV_loc, ...
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
                        FI_hrv_peaks_before, FI_hrv_peaks_after);
                    % Write to the array
                    % seconds
                    inds = find(contains(satContents, 'peak_hrSegDuration'));
                    satData(k,inds) = {peak_hrSegDuration};
                    % %time
                    inds = find(contains(satContents, 'peak_HR_loc'));
                    satData(k,inds) = {peak_HR_loc};
                    % ms
                    inds = find(contains(satContents, 'peak_rrHRV_loc'));
                    satData(k,inds) = {peak_rrHRV_loc};
                    % ms
                    inds = find(contains(satContents, 'peak_SDNN_loc'));
                    satData(k,inds) = {peak_SDNN_loc};
                    % ms
                    inds = find(contains(satContents, 'peak_SDSD_loc'));
                    satData(k,inds) = {peak_SDSD_loc};
                    % ms
                    inds = find(contains(satContents, 'peak_RMSSD_loc'));
                    satData(k,inds) = {peak_RMSSD_loc};
                    % %
                    inds = find(contains(satContents, 'peak_pNN50_loc'));
                    satData(k,inds) = {peak_pNN50_loc};
                    % AU
                    inds = find(contains(satContents, 'peak_TRI_val_loc'));
                    satData(k,inds) = {peak_TRI_val_loc};
                    % AU
                    inds = find(contains(satContents, 'peak_TINN_val_loc'));
                    satData(k,inds) = {peak_TINN_val_loc};
                    % AU
                    inds = find(contains(satContents, 'peak_DFA_1_loc'));
                    satData(k,inds) = {peak_DFA_1_loc};
                    % AU
                    inds = find(contains(satContents, 'peak_DFA_2_loc'));
                    satData(k,inds) = {peak_DFA_2_loc};
                    % ms2
                    inds = find(contains(satContents, 'peak_ApEn_loc'));
                    satData(k,inds) = {peak_ApEn_loc};
                    % %
                    inds = find(contains(satContents, 'peak_fvf_pLF'));
                    satData(k,inds) = {peak_fvf_pLF};
                    % %
                    inds = find(contains(satContents, 'peak_fvf_pHF'));
                    satData(k,inds) = {peak_fvf_pHF};
                    % AU
                    inds = find(contains(satContents, 'peak_fvf_LFHFratio'));
                    satData(k,inds) = {peak_fvf_LFHFratio};
                    % ms2
                    inds = find(contains(satContents, 'peak_fvf_VLF'));
                    satData(k,inds) = {peak_fvf_VLF};
                    % ms2
                    inds = find(contains(satContents, 'peak_fvf_LF'));
                    satData(k,inds) = {peak_fvf_LF};
                    % ms2
                    inds = find(contains(satContents, 'peak_fvf_HF'));
                    satData(k,inds) = {peak_fvf_HF};
                    % ms
                    inds = find(contains(satContents, 'peak_SD1'));
                    satData(k,inds) = {peak_SD1};
                    % ms
                    inds = find(contains(satContents, 'peak_SD2'));
                    satData(k,inds) = {peak_SD2};
                    % AU
                    inds = find(contains(satContents, 'peak_SD1SD2ratio'));
                    satData(k,inds) = {peak_SD1SD2ratio};
                else
                    disp("Cannot be computed ...")
                end
            else
                disp("Cannot be computed ...")
            end
        else
            %
        end% FI_hrv_peaks_yes == 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% E5-2.2: Calculate HRV around respiratory events
        if FI_hrv_respEv_yes == 1 & exist('respEventList');
            disp("Hypoxic burden acc. Filchenko (HRV around respEvents) ...")
            if hr_yes == 0 | isempty(respEventList) | class(respEventList) == "double"
                disp("Cannot be computed ...")
            elseif ~isnan(respEventList{1,1})
                disp("Computing ...")
                [respEv_hrSegDuration, respEv_HR_loc, respEv_rrHRV_loc, ...
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
                    FI_hrv_respEvent_before, FI_hrv_respEvent_after);
                % Write to the array
                % seconds
                inds = find(contains(satContents, 'respEv_hrSegDuration'));
                satData(k,inds) = {respEv_hrSegDuration};
                % %time
                inds = find(contains(satContents, 'respEv_HR_loc'));
                satData(k,inds) = {respEv_HR_loc};
                % ms
                inds = find(contains(satContents, 'respEv_rrHRV_loc'));
                satData(k,inds) = {respEv_rrHRV_loc};
                % ms
                inds = find(contains(satContents, 'respEv_SDNN_loc'));
                satData(k,inds) = {respEv_SDNN_loc};
                % ms
                inds = find(contains(satContents, 'respEv_SDSD_loc'));
                satData(k,inds) = {respEv_SDSD_loc};
                % ms
                inds = find(contains(satContents, 'respEv_RMSSD_loc'));
                satData(k,inds) = {respEv_RMSSD_loc};
                % %
                inds = find(contains(satContents, 'respEv_pNN50_loc'));
                satData(k,inds) = {respEv_pNN50_loc};
                % AU
                inds = find(contains(satContents, 'respEv_TRI_val_loc'));
                satData(k,inds) = {respEv_TRI_val_loc};
                % AU
                inds = find(contains(satContents, 'respEv_TINN_val_loc'));
                satData(k,inds) = {respEv_TINN_val_loc};
                % AU
                inds = find(contains(satContents, 'respEv_DFA_1_loc'));
                satData(k,inds) = {respEv_DFA_1_loc};
                % AU
                inds = find(contains(satContents, 'respEv_DFA_2_loc'));
                satData(k,inds) = {respEv_DFA_2_loc};
                % ms2
                inds = find(contains(satContents, 'respEv_ApEn_loc'));
                satData(k,inds) = {respEv_ApEn_loc};
                % %
                inds = find(contains(satContents, 'respEv_fvf_pLF'));
                satData(k,inds) = {respEv_fvf_pLF};
                % %
                inds = find(contains(satContents, 'respEv_fvf_pHF'));
                satData(k,inds) = {respEv_fvf_pHF};
                % AU
                inds = find(contains(satContents, 'respEv_fvf_LFHFratio'));
                satData(k,inds) = {respEv_fvf_LFHFratio};
                % ms2
                inds = find(contains(satContents, 'respEv_fvf_VLF'));
                satData(k,inds) = {respEv_fvf_VLF};
                % ms2
                inds = find(contains(satContents, 'respEv_fvf_LF'));
                satData(k,inds) = {respEv_fvf_LF};
                % ms2
                inds = find(contains(satContents, 'respEv_fvf_HF'));
                satData(k,inds) = {respEv_fvf_HF};
                % ms
                inds = find(contains(satContents, 'respEv_SD1'));
                satData(k,inds) = {respEv_SD1};
                % ms
                inds = find(contains(satContents, 'respEv_SD2'));
                satData(k,inds) = {respEv_SD2};
                % AU
                inds = find(contains(satContents, 'respEv_SD1SD2ratio'));
                satData(k,inds) = {respEv_SD1SD2ratio};
            else
                %
            end
        else
        end % FI_hrv_respEv_yes == 1;
        % End Step 5

        %% Step E6. Record raw data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write .mat with raw data if needed
        if save_raw_output == 1
            % Create time stamp
            format longG
            t = now;
            d = datetime(t,'ConvertFrom','datenum');
            DateString = datestr(d);
            DateString = strrep(DateString, '-', '');
            DateString = strrep(DateString, ' ', '_');
            DateString = strrep(DateString, ':', '');

            % Create the structure for the raw data
            rawData.Number = k; % number of the loop
            rawData.Filepath = fullFileNames(k); % path of the file (including the filename!!)
            inds_main = find(contains(satContents, 'RecordID'));
            RecordID = satData(k,inds_main);
            rawData.RecordID = RecordID{1,1};
            inds_main = find(contains(satContents, 'Visit'));
            Visit = satData(k,inds_main);
            rawData.Visit = Visit{1,1};
            rawData.satChannel = satChannel;
            rawData.satSignal = satSignal;
            if hr_yes == 1
                rawData.hrChannel = hrChannel;
                rawData.hrSignal = hrSignal;
            else
                %
            end
            rawData.respEventList = respEventList;

            % Write table as .mat
            tab_name = strcat(path_save_raw_output, ssep, ...
                rawData.RecordID, '_', rawData.Visit, '_', DateString, '.mat');
            save(char(tab_name), 'rawData')
        else
            %
        end % raw output

        %% Step E7. Delete temporary variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clearvars('-except',keepvars{:})

    catch
        % This covers all errors
        % Files have to be manually checked
        inds = find(contains(satContents, 'Error_undefined'));
        satData(k,inds) = {'Undefined error in file'};
        clearvars('-except',keepvars{:})
    end % catch

end

%% Save output
% Create time stamp
format longG
t = now;
d = datetime(t,'ConvertFrom','datenum');
DateString = datestr(d);
DateString = strrep(DateString, '-', '');
DateString = strrep(DateString, ' ', '_');
DateString = strrep(DateString, ':', '');

% Save the output file
% Create table
satTable = cell2table(satData, 'VariableNames', satContents);

% Name table
tab_name = strcat(path_save_output_matrix, 'HB_output_', DateString, '.csv');
writetable(satTable, tab_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
