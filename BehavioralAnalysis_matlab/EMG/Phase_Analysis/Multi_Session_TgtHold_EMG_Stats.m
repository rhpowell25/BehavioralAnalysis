
%% Loading the morning and afternoon files
close all
clear
clc

% Monkey Name
Monkey = 'Tot';

% Drug
Drug_Choice = 'Lex';

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

% Only use trial tasks selected
disp('Powergrasp Only')
powergrasp_idx = find(contains(Tasks, 'PG'));
Dates = Dates(powergrasp_idx);
Tasks = Tasks(powergrasp_idx);

% Sorted or unsorted (1 vs 0)
Sorted = 1;

EMG_Names = struct([]);
TgtHold_EMG_p_values = struct([]);
TgtHold_EMG_perc_changes = struct([]);
Morn_TgtHold_EMG = struct([]);
Noon_TgtHold_EMG = struct([]);

%% Loop through each xds file
for xx = 1:length(Dates)

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx}, Tasks{xx}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx}, Tasks{xx}, Sorted, 'Noon');
    
    % Process the xds files
    Match_The_Targets = 0;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);
    
    %% Plotting Parameters 

    % Zero? (1 = Yes, 0 = No)
    zero_EMG = 1;
    zero_method = 'Prctile';
    
    % Normalize? (1 = Yes, 0 = No)
    norm_EMG = 1;
    norm_perctile = 95;
    
    %% EMG
    
    % All, Flex, Exten, Both, or Custom?
    muscle_groups = 'Grasp';
    
    EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
    EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

    %% EMG Statistics
    % Target hold phase
    [EMG_Names{xx}, TgtHold_EMG_p_values{xx}, TgtHold_EMG_perc_changes{xx}, Morn_TgtHold_EMG{xx}, ~, Noon_TgtHold_EMG{xx}, ~] = ...
        TgtHold_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
        EMG_Zero_Factor, EMG_Norm_Factor, 1, 0);
    
end
