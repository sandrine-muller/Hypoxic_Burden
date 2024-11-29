function [record_id, Error_with_record] = f_extract_record(edfName, edfPath)
% This function extracts record_id based on the path
% or the edf name.
% This algorithm works for eSATIS and SDSO.
Error_with_record = [];
record_id = [];

% Select meaningful directories
pathParts = strsplit(edfPath, filesep);
first_topDir = pathParts(length(pathParts));
second_topDir = pathParts(length(pathParts)-1);

% Search for specific patterns
pattern_search = {'C01S', 'C02S', 'C03S', 'C04S', 'C05S', 'C06S'};
if find(contains(edfName, pattern_search))
  correct_answer = '';
    % Iterate through each pattern
    for i = 1:length(pattern_search)
        pattern = pattern_search{i};
        % Check if the pattern is present in edfName
        if ~isempty(strfind(edfName, pattern))
            correct_answer = pattern;
            break; % Exit the loop if a match is found
        end
    end
    pattern = correct_answer + digitsPattern;
    record_id = extract(edfName, pattern);
    Error_with_record = 0;  
elseif find(contains(first_topDir, pattern_search)) % next we try the first top directory
    res = first_topDir{1,1};
    res = strrep(res, 'AL_', '');
    res = strrep(res, 'RP_', '');
    if contains(res, '_v');
        res = extractBefore(res, '_v');
    else
    end
    if contains(res, '_V');
        res = extractBefore(res, '_V');
    else
    end
    record_id = res;
    record_id = {record_id};
elseif find(contains(second_topDir, pattern_search)) % next we try the second top directory
    correct_answer = '';
    % Iterate through each pattern
    for i = 1:length(pattern_search)
        pattern = pattern_search{i};
        % Check if the pattern is present in edfName
        if ~isempty(strfind(edfName, pattern))
            correct_answer = pattern;
            break; % Exit the loop if a match is found
        end
    end
    pattern = correct_answer + digitsPattern;
    record_id = extract(topDir, pattern);
    Error_with_record = 0; 
else
    record_id = {edfName};
    % Find error index
    Error_with_record = 1; % is not very meaningful, and should be changed manually
end

end