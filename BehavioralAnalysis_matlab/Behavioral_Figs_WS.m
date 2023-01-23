
%% File Description:

% This function plots and saves behavioral summaries for the morning &
% afternoon xds files. These plots include reaction time box plots, trial
% averaged EMG plots, EMG statistics, consecutive EMG plots, & cursor position plots

%% Loading the morning and afternoon files
clear
close all
clc

% Select The Date & Task To Analyze
Date = '20210610';
Task = 'WS';

[xds_morn, xds_noon] = Load_XDS(Date, Task, Process_XDS);

% Process the xds files
Match_The_Targets = 0;
[xds_morn, xds_noon] = Process_XDS(xds_morn, xds_noon, Match_The_Targets);


File_Name_morn = xds_morn.meta.rawFileName;

%% Save directory

% Date
File_Name = xds_noon.meta.rawFileName;
nondate_info = extractAfter(File_Name, '_');
Date = erase(File_Name, strcat('_', nondate_info));
% Monkey
nonmonkey_info = extractAfter(nondate_info, '_');
% Task
nontask_info = extractAfter(nonmonkey_info, '_');
Task = erase(nonmonkey_info, strcat('_', nontask_info));
% Drug
if contains(nontask_info, 'Caff')
    Drug = 'Caffeine';
end
if contains(nontask_info, 'Lex')
    Drug = 'Escitalopram';
end
if contains(nontask_info, 'Cyp')
    Drug = 'Cypro';
end
if contains(nontask_info, 'Con')
    Drug = 'Control';
end

% Define the save folder
save_folder = strcat(Date, '_', Task);

% Define the save directory
drug_save_dir = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Figures\', Monkey, '_', Drug, '\');
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

%% Find the number of targets in the morning and afternoon

[dir_idx_morn, unique_tgts_morn, ~, ~, ~, ~] = ... 
    WindowTrialGoCueFiringRate(xds_morn, 1);
[dir_idx_noon, unique_tgts_noon, ~, ~, ~, ~] = ... 
    WindowTrialGoCueFiringRate(xds_noon, 1);

% Check to see if both sessions use a consistent number of directions
if ~isequal(unique(dir_idx_morn), unique(dir_idx_noon))
    disp('Uneven Target Directions Between Morning & Afternoon');
    % Only use the info of target directions conserved between morn & noon
    shared_dir_idx = ismember(dir_idx_morn, dir_idx_noon);
    dir_idx_morn = dir_idx_morn(shared_dir_idx);
    unique_tgts_morn = unique_tgts_morn(shared_dir_idx);
    dir_idx_noon = dir_idx_noon(shared_dir_idx);
    unique_tgts_noon = unique_tgts_noon(shared_dir_idx);
end

% Check to see if both sessions use a consistent number of targets
if ~isequal(unique(unique_tgts_morn), unique(unique_tgts_noon))
    disp('Uneven Target Centers Between Morning & Afternoon');
    % Only use the info of target centers conserved between morn & noon
    shared_target_centers_idx = ismember(unique_tgts_morn, unique_tgts_noon);
    dir_idx_morn = dir_idx_morn(shared_target_centers_idx);
    unique_tgts_morn = unique_tgts_morn(shared_target_centers_idx);
    dir_idx_noon = dir_idx_noon(shared_target_centers_idx);
    unique_tgts_noon = unique_tgts_noon(shared_target_centers_idx);
end

%% Reaction Time Folders

save_dir = strcat(trial_save_dir, 'Reaction Times\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

Rxn_Time_BoxPlot(xds_morn, xds_noon, 0)

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

%% Normalize Cursor? (1 = Yes, 0 = No)
norm_cursor = 1;
Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_perctile, norm_cursor);

for kk = 1:length(dir_idx_morn)
    %% EMG Folders

    if isequal(dir_idx_morn(kk), 0) && strcmp(xds_morn.meta.hand, 'Left')
        emg_save_dir = strcat(trial_save_dir, 'Flexor EMG\');
        if ~exist(emg_save_dir, 'dir')
            mkdir(emg_save_dir);
        end

        % All, Flex, Exten, Uln_Dev, Rad_Dev, Grasp, or Custom?
        muscle_groups = 'Flex';

        EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
        EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

    end

    if isequal(dir_idx_morn(kk), 90)
        emg_save_dir = strcat(trial_save_dir, 'Rad Dev EMG\');
        if ~exist(emg_save_dir, 'dir')
            mkdir(emg_save_dir);
        end
        
        % All, Flex, Exten, Uln_Dev, Rad_Dev, Grasp, or Custom?
        muscle_groups = 'Rad_Dev';

        EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
        EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

    end

    if isequal(dir_idx_morn(kk), 180) && strcmp(xds_morn.meta.hand, 'Left')
        emg_save_dir = strcat(trial_save_dir, 'Extensor EMG\');
        if ~exist(emg_save_dir, 'dir')
            mkdir(emg_save_dir);
        end

        % All, Flex, Exten, Uln_Dev, Rad_Dev, Grasp, or Custom?
        muscle_groups = 'Exten';

        EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
        EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

    end

    if isequal(dir_idx_morn(kk), -90)
        emg_save_dir = strcat(trial_save_dir, 'Uln Dev EMG\');
        if ~exist(emg_save_dir, 'dir')
            mkdir(emg_save_dir);
        end

        % All, Flex, Exten, Uln_Dev, Rad_Dev, Grasp, or Custom?
        muscle_groups = 'Uln_Dev';

        EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);
        EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

    end  

    %% Consecutive Raw EMG

    % Morning
    ConsecTrialsRawEMG(xds_morn, trial_num, muscle_groups, 0)

    % Save Figures
    save_title = strcat('Morning', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Raw EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'fig')
    end
    close gcf

    % Afternoon
    ConsecTrialsRawEMG(xds_noon, trial_num, muscle_groups, 0)

    % Save Figures
    save_title = strcat('Afternoon', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, Raw EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'fig')
    end
    close gcf

    %% Consecutive EMG

    % Morning
    ConsecTrialsEMG(xds_morn, EMG_Zero_Factor, EMG_Norm_Factor, trial_num, muscle_groups, 0)

    % Save Figures
    save_title = strcat('Morning', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'fig')
    end
    close gcf

    % Afternoon
    ConsecTrialsEMG(xds_noon, EMG_Zero_Factor, EMG_Norm_Factor, trial_num, muscle_groups, 0)

    % Save Figures
    save_title = strcat('Afternoon', {' '}, num2str(trial_num), {' '}, 'Consecutive Succesful Trials, EMG', {' '}, muscle_groups);
    if ~strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), Save_Figs)
    end
    if strcmp(Save_Figs, 'All')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'png')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'pdf')
        saveas(gcf, fullfile(emg_save_dir, char(save_title)), 'fig')
    end
    close gcf

    %% Define the Event

    for jj = 1:2

        if isequal(jj, 1)
            save_dir = strcat(emg_save_dir, 'Aligned Trial Go Cue\');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            % Select the event to align to:
            event = 'trial_gocue';
        end

        if isequal(jj, 2)
            save_dir = strcat(emg_save_dir, 'Aligned Trial End\');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            % Select the event to align to:
            event = 'trial_end';
        end

        %% EMG Statistics
        if strcmp(event, 'trial_gocue')
            % Baseline EMG Statistics
            [~, ~] = Baseline_EMG_Stats(xds_morn, xds_noon, muscle_groups, EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, 0);
        end
        if strcmp(event, 'trial_end')
            % Trial End EMG Statistics
            [~, ~] = TgtHold_EMG_Stats(xds_morn, xds_noon, muscle_groups, EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, 0);
        end

        % Save Figures
        for ii = 1:length(findobj('type','figure'))
            fig_info = get(subplot(211),'title');
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

        %% Per Trial EMG

        % Morning
        Per_Trial_EMG(xds_morn, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0)

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
        Per_Trial_EMG(xds_noon, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0)

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

        OverlapPlotEMG(xds_morn, xds_noon, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0)

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

end % End of EMG loop


%% Cursor

cursor_pos_YLims = CursorPosYLimit(xds_morn, xds_noon);

curs_save_dir = strcat(trial_save_dir, 'Cursor\');
if ~exist(curs_save_dir, 'dir')
    mkdir(fullfile(strcat(trial_save_dir, 'Cursor\')));
end


%% Define the Event

for jj = 1:2

    if isequal(jj, 1)
        cursor_save_dir = strcat(curs_save_dir, 'Aligned Trial Go Cue\');
        if ~exist(cursor_save_dir, 'dir')
            mkdir(cursor_save_dir);
        end

        % Select the event to align to:
        event = 'trial_gocue';
    end

    if isequal(jj, 2)
        cursor_save_dir = strcat(curs_save_dir, 'Aligned Trial End\');
        if ~exist(cursor_save_dir, 'dir')
            mkdir(cursor_save_dir);
        end

        % Select the event to align to:
        event = 'trial_end';
    end

    %% Cursor Statistics
    if strcmp(event, 'trial_gocue')
        % Cursor Statistics
        [~, ~, ~, ~, ~, ~] = Baseline_CursorPos_Stats(xds_morn, xds_noon, Cursor_Norm_Factor, Plot_Figs, 0);

        % Save Figures
        for ii = 1:length(findobj('type','figure'))
            fig_info = get(subplot(211),'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
            end
            close gcf
        end

    end

    if strcmp(event, 'trial_end')
        % Cursor Statistics
        [~, ~, ~, ~, ~, ~] = TgtHold_CursorPos_Stats(xds_morn, xds_noon,  Cursor_Norm_Factor, Plot_Figs, 0);

        % Save Figures
        for ii = 1:length(findobj('type','figure'))
            fig_info = get(subplot(211),'title');
            save_title = get(fig_info, 'string');
            save_title = strrep(save_title, ':', '');
            save_title = strrep(save_title, '.', '_');
            if ~strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
            end
            if strcmp(Save_Figs, 'All')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
                saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
            end
            close gcf
        end

    end

    %% Per Trial Cursor

    % Morning
    Per_Trial_CursorPos(xds_morn, event, unit_name, cursor_pos_YLims, 0)

    % Save Figures
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        save_title = strcat('Morning', {' '}, save_title);
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    % Afternoon
    Per_Trial_CursorPos(xds_noon, event, unit_name, cursor_pos_YLims, 0)

    % Save Figures
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        save_title = strcat('Afternoon', {' '}, save_title);
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    %% Average Cursor

    % Morning
    PlotCursorPos(xds_morn, unit_name, event, cursor_pos_YLims, Plot_Figs, 0)

    % Save Figures
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        save_title = strcat('Morning', {' '}, save_title);
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    % Afternoon
    PlotCursorPos(xds_noon, unit_name, event, cursor_pos_YLims, Plot_Figs, 0)

    % Save Figures
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        save_title = strcat('Afternoon', {' '}, save_title);
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

    %% Overlap Average Cursor

    OverlapPlotCursor(xds_morn, xds_noon, event, 0)

    % Save Figures
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(cursor_save_dir, char(save_title)), 'fig')
        end
        close gcf
    end

end % End of the event loop














