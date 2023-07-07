
%% File Description:

% This function plots and saves behavioral summaries for the morning &
% afternoon xds files. These plots include reaction time box plots, trial
% averaged EMG plots, EMG statistics, consecutive EMG plots, & force plots

%% Loading the morning and afternoon files
clear
close all
clc

% Monkey Name
Monkey = 'Tot';

% Drug
Drug_Choice = 'Caff';

% Load the file information
[Dates, Tasks, ~] = File_Details(Monkey, Drug_Choice);

% Only use trial tasks selected
disp('Powergrasp Only')
powergrasp_idx = find(contains(Tasks, 'PG'));
Dates = Dates(powergrasp_idx);
Tasks = Tasks(powergrasp_idx);

% Sorted or unsorted (1 vs 0)
Sorted = 1;

%% Loop through each xds file
for xx = 1:length(Dates)

    % Load the xds files
    xds_morn = Load_XDS(Monkey, Dates{xx}, Tasks{xx}, Sorted, 'Morn');
    xds_noon = Load_XDS(Monkey, Dates{xx}, Tasks{xx}, Sorted, 'Noon');
    
    % Process the xds files
    Match_The_Targets = 0;
    [xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);
    
    %% Save directory
    
    % Define the save folder
    save_folder = strcat(Dates{xx}, '_', Tasks{xx});
    
    % Define the save directory
    drug_save_dir = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Figures\', Monkey, '_', Drug_Choice, '\');
    trial_save_dir = strcat(drug_save_dir, save_folder, '\Behavioral Figs\');
    if ~exist(trial_save_dir, 'dir')
        mkdir(fullfile(trial_save_dir));
    end
    
    %% Plotting & Saving Parameters 
    
    % Select the unit of interest ('elec1_1, #, 'All', or nan)
    unit_name = nan;
    
    % Decide whether or not to plot (1 = Yes; 0 = No)
    Plot_Figs = 1;
    % Save the figures to desktop? ('pdf', 'png', 'fig', 0 = No)
    Save_Figs = 'png';
    
    % Zero? (1 = Yes, 0 = No)
    zero_EMG = 1;
    zero_method = 'Prctile';
    
    % Normalize? (1 = Yes, 0 = No)
    norm_EMG = 1;
    norm_perctile = 95;
    
    % Number of Trials to Plot (# or 'All')
    trial_num = 10;
    
    %% Reaction Time Folders
    
    save_dir = strcat(trial_save_dir, 'Reaction Times\');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    Task_Metric_ViolinPlot(xds_morn, xds_noon, 'Rxn_Time', 0)
    
    % Save Figures
    for ii = 1:length(findobj('type','figure'))
        save_title = strcat('Reaction Time', Dates{xx}, '_', Monkey, '_', Tasks{xx});
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
        end
        close gcf
    end
    
    %% EMG Folders
    save_dir = strcat(trial_save_dir, 'Extrinsic Hand EMG\');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % All, Flex, Exten, Both, or Custom?
    muscle_groups = 'Grasp';
    
    EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
    EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);
    
    %% Consecutive Raw EMG
    
    % Morning
    ConsecTrialsRawEMG(xds_morn, trial_num, muscle_groups, 0)
    
    % Save Figures
    save_title = strcat('Morning', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Raw EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf
    
    % Afternoon
    ConsecTrialsRawEMG(xds_noon, trial_num, muscle_groups, 0)
    
    % Save Figures
    save_title = strcat('Afternoon', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Raw EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf
    
    %% Consecutive EMG
    
    % Morning
    ConsecTrialsEMG(xds_morn, EMG_Zero_Factor, EMG_Norm_Factor, trial_num, muscle_groups, 0)
    
    % Save Figures
    save_title = strcat('Morning', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf
    
    % Afternoon
    ConsecTrialsEMG(xds_noon, EMG_Zero_Factor, EMG_Norm_Factor, trial_num, muscle_groups, 0)
    
    % Save Figures
    save_title = strcat('Afternoon', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf

    %% EMG Statistics
    % Baseline phase
    [EMG_Names, Baseline_EMG_p_values, Baseline_EMG_perc_changes, ~, ~, ~, ~] = ...
        Baseline_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
        EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, 0);
    % Ramp phase
    [~, Ramp_EMG_p_values, Ramp_EMG_perc_changes, ~, ~, ~, ~] = ...
        Ramp_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
        EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, 0);
    % Target hold phase
    [~, TgtHold_EMG_p_values, TgtHold_EMG_perc_changes, ~, ~, ~, ~] = ...
        TgtHold_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
        EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, 0);
  
    stat_save_dir = strcat(save_dir, 'EMG_Stats\');
    if ~exist(stat_save_dir, 'dir')
        mkdir(stat_save_dir);
    end

    % Save Figures
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(subplot(211),'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    %% Define the Event
    
    for jj = 1:2
    
        if isequal(jj, 1)
            mkdir(fullfile(strcat(trial_save_dir, 'Extrinsic Hand EMG\Trial Go Cue\')));
            save_dir = strcat(trial_save_dir, 'Extrinsic Hand EMG\Trial Go Cue\');
    
            % Select the event to align to:
            event = 'trial_gocue';
        end
    
        if isequal(jj, 2)
            mkdir(fullfile(strcat(trial_save_dir, 'Extrinsic Hand EMG\Trial End\')));
            save_dir = strcat(trial_save_dir, 'Extrinsic Hand EMG\Trial End\');
    
            % Select the event to align to:
            event = 'trial_end';
        end

        %% Per Trial EMG

        % Morning
        Per_Trial_EMG(xds_morn, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(subplot(211),'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Morning', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        % Afternoon
        Per_Trial_EMG(xds_noon, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(subplot(211),'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Afternoon', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end

        %% Average EMG
    
        % Morning
        PlotEMG(xds_morn, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, Plot_Figs, 0);
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Morning', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        % Afternoon
        PlotEMG(xds_noon, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, Plot_Figs, 0);
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Afternoon', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        %% Overlap Average EMG
    
        OverlapPlotEMG(xds_morn, xds_noon, event, 1, 1, muscle_groups, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
    end % End of the event loop
    
    
    %% Force
    
    % Normalize Force? (1 = Yes, 0 = No, 'Convert')
    norm_force = 'Convert';
    Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force);
    
    force_YLims = ForceYLimit(xds_morn, xds_noon, norm_force);
    
    mkdir(fullfile(strcat(trial_save_dir, 'Force\')));
    save_dir = strcat(trial_save_dir, 'Force\');
    
    %% Consecutive Force
    
    % Morning
    ConsecPlotForce(xds_morn, Force_Norm_Factor, trial_num, 0)
    
    % Save Figures
    save_title = strcat('Morning', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Force', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf
    
    % Afternoon
    ConsecPlotForce(xds_noon, Force_Norm_Factor, trial_num, 0)
    
    % Save Figures
    save_title = strcat('Afternoon', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Force', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf
    
    %% Force Sensor plotting
    
    Target_Force_Traces(xds_morn, xds_noon, norm_force, 0)
    
    % Save Figures
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
        end
        close gcf
    end
    
    Decomposed_Force_Scatter(xds_morn, xds_noon, norm_force, 0)

    % Save Figures
    save_title = strcat(Dates{xx}, '_', Tasks{xx}, '_', Drug_Choice, '_', 'Decomposed Force');
    save_title = strrep(save_title, ':', '');
    save_title = strrep(save_title, '.', '_');
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
    end
    close gcf

    Force_Sensor_Ratios(xds_morn, xds_noon, norm_force, 0)

    % Save Figures
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    %% Force Statistics

    % Ramp phase
    [Ramp_Force_p_values, Ramp_Force_perc_changes, ...
        Morn_Ramp_Force, Err_Morn_Ramp_Force, Noon_Ramp_Force, Err_Noon_Ramp_Force] = ...
        Ramp_Force_Stats(xds_morn, xds_noon, norm_force, Plot_Figs, 0);

    % Target hold phase
    [TgtHold_Force_p_values, TgtHold_Force_perc_changes, ...
        Morn_TgtHold_Force, Err_Morn_TgtHold_Force, Noon_TgtHold_Force, Err_Noon_TgtHold_Force] = ...
        TgtHold_Force_Stats(xds_morn, xds_noon, norm_force, Plot_Figs, 0);
  
    stat_save_dir = strcat(save_dir, 'Force_Stats\');
    if ~exist(stat_save_dir, 'dir')
        mkdir(stat_save_dir);
    end

    % Save Figures
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(subplot(211),'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(stat_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end
    
    %% Define the Event
    
    for jj = 1:2
    
        if isequal(jj, 1)
            mkdir(fullfile(strcat(trial_save_dir, 'Force\Trial Go Cue\')));
            save_dir = strcat(trial_save_dir, 'Force\Trial Go Cue\');
    
            % Select the event to align to:
            event = 'trial_gocue';
        end
    
        if isequal(jj, 2)
            mkdir(fullfile(strcat(trial_save_dir, 'Force\Trial End\')));
            save_dir = strcat(trial_save_dir, 'Force\Trial End\');
    
            % Select the event to align to:
            event = 'trial_end';
        end
    
        %% Per Trial Force
    
        % Morning
        Per_Trial_Force(xds_morn, event, unit_name, Force_Norm_Factor, force_YLims, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Morning', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        % Afternoon
        Per_Trial_Force(xds_noon, event, unit_name, Force_Norm_Factor, force_YLims, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Afternoon', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        %% Average Force
    
        % Morning
        PlotForce(xds_morn, event, unit_name, Force_Norm_Factor, force_YLims, Plot_Figs, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Morning', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        % Afternoon
        PlotForce(xds_noon, event, unit_name, Force_Norm_Factor, force_YLims, Plot_Figs, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            save_title = strcat('Afternoon', {' '}, save_title);
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
        %% Overlap Average Force
    
        OverlapPlotForce(xds_morn, xds_noon, event, 0)
    
        % Save Figures
        for ii = 1:numel(findobj('type','figure'))
            fig_info = get(gca,'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            save_title = strrep(save_title, '/', '_');
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
            end
            close gcf
        end
    
    end % End of the event loop
end














