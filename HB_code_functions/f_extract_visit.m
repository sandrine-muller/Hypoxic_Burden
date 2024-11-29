function [visit] = f_extract_visit(edfName, edfPath)
% This function extracts visit from the edfName or edfPath
% It works for eSATIS and SDSO
visit = [];

% Check filepath
fpath = edfPath;
if contains(fpath, 'visit')
    fpath = extractAfter(fpath,'visit');
    fpath = regexp(fpath,'[0-9]','match');
    fpath = fpath{1,1};
elseif contains(fpath, '_V')
    fpath = extractAfter(fpath,'_V');
    fpath = regexp(fpath,'[0-9]','match');
    fpath = fpath{1,1};
else
    fpath = 'Unknown'; % assume that this is visit 1
    % fpath = '0'; % assume that the visit is unknown
end

% Check filename
if find(contains(edfName, {'V1' 'visit_1'}))
    visit = '1';
elseif find(contains(edfName, {'V3' 'visit_3'}))
    visit = '3';
elseif find(contains(edfName, {'V4' 'visit_4' 'V5'}))
    visit = '4';
elseif find(contains(edfName, {'V6'}))
    visit = '6';
else
    visit = 'Unknown'; % Indicates files that have to be renamed --> most likely Apnealinks --> always V1
end

% Set priority
if fpath ~= visit & visit == 'Unknown'
    visit = fpath;
else
    %
end

% Format
if fpath == 'Unknown' & visit == 'Unknown'
    visit = strcat('V1_', visit);
else
    visit = strcat('V', visit);
end

end