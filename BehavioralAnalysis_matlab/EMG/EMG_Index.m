function [M] = EMG_Index(xds, muscle_groups)

%% File Description:

% This function finds the indices of relevant EMG in an XDS file based on
% the muscle groups you are interested in looking at
%
% -- Inputs --
% xds: the XDS file in question
% muscle_groups: 'Flex', 'Ext', Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All'

%% If there is no EMG

if ~xds.has_EMG
    disp('No EMG in this file');
    M = [];
    return
end

%% Define the muscles 

muscle_names = strings;

if strcmp(muscle_groups, 'Flex')
    muscle_names(1) = 'FCR';
    muscle_names(2) = 'FCU';
end

if strcmp(muscle_groups, 'Ext')
    muscle_names(1) = 'ECR';
    muscle_names(2) = 'ECU';
end

if strcmp(muscle_groups, 'Uln_Dev')
    muscle_names(1) = 'FCU';
    muscle_names(2) = 'ECU';
end

if strcmp(muscle_groups, 'Rad_Dev')
    muscle_names(1) = 'FCR';
    muscle_names(2) = 'ECR';
end

if strcmp(muscle_groups, 'Both')
    muscle_names(1) = 'FCR';
    muscle_names(2) = 'FCU';
    muscle_names(3) = 'ECR';
    muscle_names(4) = 'ECU';
end

if strcmp(muscle_groups, 'Grasp')
    muscle_names(1) = 'FDP';
    muscle_names(2) = 'FDS';
end

if strcmp(muscle_groups, 'Pinch')
    muscle_names(1) = '1DI';
    muscle_names(2) = 'APB';
    muscle_names(3) = 'FPB';
end

if strcmp(muscle_groups, 'Custom')
    muscle_names(1) = 'FCR';
    muscle_names(2) = 'FDP';
    muscle_names(3) = 'FDS';
    muscle_names(4) = 'EDC';
end

% If muscle_names is still blank
if (muscle_names == "")
    muscle_names(1) = muscle_groups;
end

% Find the indices of the muscles of interest
muscle_idx = struct([]);
cc = 0;
for ii = 1:length(muscle_names)
    muscle_idx{ii,1} = find(contains(xds.EMG_names, muscle_names(ii)));
    % Find how many muscles there are
    cc = cc + length(muscle_idx{ii,1});
end

% Concatenate the indices
M = zeros(cc,1);
cc = 1;
for ii = 1:length(muscle_idx)
    for jj = 1:length(muscle_idx{ii})
        M(cc) = muscle_idx{ii,1}(jj);
        cc = cc + 1;
    end
end

if strcmp(muscle_groups, 'All')
    M = 1:length(xds.EMG_names);
    M = M';
end

