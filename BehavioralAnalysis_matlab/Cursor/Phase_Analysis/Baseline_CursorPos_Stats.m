function [Baseline_CursorPos_p_values, Baseline_CursorPos_perc_changes, ...
    Morn_Baseline_CursorPos, Err_Morn_Baseline_CursorPos, Noon_Baseline_CursorPos, Err_Noon_Baseline_CursorPos] = ...
    Baseline_CursorPos_Stats(xds_morn, xds_noon, norm_cursor, Plot_Figs, Save_Figs)

%% File Description:

% This function finds changes in cursor position during the center hold 
% between two XDS files & plots a boxplot & individual trial traces.
% If you set Plot_Figs to 0, the figure will not be plotted.
% If you set Save_Figs to 0, the figure will not be saved to your desktop.
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% Cursor_Norm_Factor: the factor to normalize the cursor position
% Plot_Figs: 1 or 0
% Save_Figs: 'pdf', 'png', 'fig', or 0

%% Display the function being used
disp('Baseline Cursor Position Statistics:');

%% Basic Settings, some variable extractions, & definitions

% Font specifications
legend_font_size = 12;
title_font_size = 15;
font_name = 'Arial';

% Define the window for the baseline phase
time_before_gocue = 0.4;

% Do you want a baseline cursor for each target / direction combo? (1 = Yes, 0 = No)
per_dir_curs = 0;

%% Extract the target directions & centers
[target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
[target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

%% Check to see if both sessions use a consistent number of targets

% Find matching targets between the two sessions
[Matching_Idxs_Morn, Matching_Idxs_Noon] = ...
    Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon);

% Only use the info of target centers conserved between morn & noon
if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
    disp('Uneven Targets Between Morning & Afternoon');
    target_centers_morn = target_centers_morn(Matching_Idxs_Morn);
    target_centers_noon = target_centers_noon(Matching_Idxs_Noon);
    target_dirs_morn = target_dirs_morn(Matching_Idxs_Morn);
    target_dirs_noon = target_dirs_noon(Matching_Idxs_Noon);
end

%% Settings to loop through every target direction

% Counts the number of directions used
num_dirs = length(target_dirs_morn);

% Save Counter
if ~isequal(Save_Figs, 0)
    close all
    save_title = strings;
end

%% Begin the loop through all directions
for jj = 1:num_dirs

    %% Times for rewarded trials
    if ~isequal(per_dir_curs, 0)
        [rewarded_gocue_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_gocue');
        [rewarded_gocue_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_gocue');
    else
        [rewarded_gocue_time_morn] = EventAlignmentTimes(xds_morn, NaN, NaN, 'trial_gocue');
        [rewarded_gocue_time_noon] = EventAlignmentTimes(xds_noon, NaN, NaN, 'trial_gocue');
    end

    %% Define the output variables
    if jj == 1 && ~isequal(per_dir_curs, 0)
        Morn_Baseline_CursorPos = zeros(1, num_dirs);
        Err_Morn_Baseline_CursorPos = zeros(1, num_dirs);
        Noon_Baseline_CursorPos = zeros(1, num_dirs);
        Err_Noon_Baseline_CursorPos = zeros(1, num_dirs);
        Baseline_CursorPos_p_values = zeros(1, num_dirs);
        Baseline_CursorPos_perc_changes = zeros(1, num_dirs);
    elseif jj == 1 && ~isequal(per_dir_curs, 1)
        Morn_Baseline_CursorPos = NaN;
        Err_Morn_Baseline_CursorPos = NaN;
        Noon_Baseline_CursorPos = NaN;
        Err_Noon_Baseline_CursorPos = NaN;
        Baseline_CursorPos_p_values = NaN;
        Baseline_CursorPos_perc_changes = NaN;
    end

    %% Cursor position aligned to go cue

    idx_length = time_before_gocue / xds_morn.bin_width;

    % Cursor position during each succesful trial
    cursor_p_morn = struct([]); % Cursor position at the start
    for ii = 1:length(rewarded_gocue_time_morn)
        idx_morn = find(xds_morn.time_frame == rewarded_gocue_time_morn(ii));
        cursor_p_morn{ii, 1} = xds_morn.curs_p(idx_morn - idx_length : idx_morn, :);
    end
    
    cursor_p_noon = struct([]); % Cursor position at the start
    for ii = 1:length(rewarded_gocue_time_noon)
        idx_noon = find(xds_noon.time_frame == rewarded_gocue_time_noon(ii));
        cursor_p_noon{ii, 1} = xds_noon.curs_p(idx_noon - idx_length : idx_noon, :);
    end

    %% Find the vector sum of the cursor position
    z_cursor_p_morn = struct([]);
    for ii = 1:length(rewarded_gocue_time_morn)
        z_cursor_p_morn{ii,1} = sqrt(cursor_p_morn{ii,1}(:, 2).^2 + cursor_p_morn{ii, 1}(:, 1).^2);
    end

    z_cursor_p_noon = struct([]);
    for ii = 1:length(rewarded_gocue_time_noon)
        z_cursor_p_noon{ii,1} = sqrt(cursor_p_noon{ii,1}(:, 2).^2 + cursor_p_noon{ii, 1}(:, 1).^2);
    end
    
    %% Putting all succesful trials in one array
    all_trials_z_cursor_p_morn = zeros(length(z_cursor_p_morn{1,1}), length(rewarded_gocue_time_morn));
    for ii = 1:length(rewarded_gocue_time_morn)
        all_trials_z_cursor_p_morn(:,ii) = z_cursor_p_morn{ii, 1};
    end

    all_trials_z_cursor_p_noon = zeros(length(z_cursor_p_noon{1,1}), length(rewarded_gocue_time_noon));
    for ii = 1:length(rewarded_gocue_time_noon)
        all_trials_z_cursor_p_noon(:,ii) = z_cursor_p_noon{ii, 1};
    end

    %% Normalizing the average cursor position
    Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_cursor);
    
    all_trials_z_cursor_p_morn = all_trials_z_cursor_p_morn / Cursor_Norm_Factor*100;
    all_trials_z_cursor_p_noon = all_trials_z_cursor_p_noon / Cursor_Norm_Factor*100;

    %% Calculating average cursor position (Average per trial)
    per_trial_avg_z_cursor_p_morn = zeros(length(cursor_p_morn), 1);
    per_trial_avg_z_cursor_p_noon = zeros(length(cursor_p_noon), 1);
    for ii = 1:length(cursor_p_morn)
        per_trial_avg_z_cursor_p_morn(ii, 1) = mean(all_trials_z_cursor_p_morn(:, ii));
    end
    for ii = 1:length(cursor_p_noon)
        per_trial_avg_z_cursor_p_noon(ii, 1) = mean(all_trials_z_cursor_p_noon(:, ii));
    end

    %% Do the statistics on the cursor position
    [~, Baseline_CursorPos_p_values(1, jj)] = ...
        ttest2(per_trial_avg_z_cursor_p_morn, per_trial_avg_z_cursor_p_noon);
    % Average cursor position
    Morn_Baseline_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_morn);
    Noon_Baseline_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_noon);
    % Standard deviation
    Baseline_cursor_p_morn_std = std(per_trial_avg_z_cursor_p_morn);
    Baseline_cursor_p_noon_std = std(per_trial_avg_z_cursor_p_noon);
    % Standard error
    Err_Morn_Baseline_CursorPos(1, jj) = Baseline_cursor_p_morn_std / sqrt(length(per_trial_avg_z_cursor_p_morn));
    Err_Noon_Baseline_CursorPos(1, jj) = Baseline_cursor_p_noon_std / sqrt(length(per_trial_avg_z_cursor_p_noon));
    % Cursor position percent change
    if all(Cursor_Norm_Factor + 1)
        Baseline_CursorPos_perc_changes(1, jj) = ...
            (Noon_Baseline_CursorPos(1, jj) - Morn_Baseline_CursorPos(1, jj)) / 100;
    else
        Baseline_CursorPos_perc_changes(1, jj) = ...
            (Noon_Baseline_CursorPos(1, jj) - Morn_Baseline_CursorPos(1, jj)) / abs(Morn_Baseline_CursorPos(1, jj));
    end
    
    %% Plot the cursor position box & whisker plot
    if isequal(Plot_Figs, 1)

        % Top Plot
        % Combining the morning & afternoon into one matrix
        avg_cursor_p = [per_trial_avg_z_cursor_p_morn; per_trial_avg_z_cursor_p_noon];
        baseline_cursor_p_labels = zeros(length(avg_cursor_p),1);
        baseline_cursor_p_labels(length(per_trial_avg_z_cursor_p_morn) + 1:end) = 1;

        Baseline_cursor_p_mean(1, 1) = Morn_Baseline_CursorPos;
        Baseline_cursor_p_mean(1, 2) = Noon_Baseline_CursorPos;
        Baseline_cursor_p_std(1, 1) = Baseline_cursor_p_morn_std;
        Baseline_cursor_p_std(1, 2) = Baseline_cursor_p_noon_std;

        % Finding the min and max for the y-axis limit
        Baseline_y_min = Baseline_cursor_p_mean - Baseline_cursor_p_std;
        Baseline_y_max = Baseline_cursor_p_mean + Baseline_cursor_p_std;

        figure
        subplot(211) 
        hold on
        boxplot(avg_cursor_p, baseline_cursor_p_labels)
    
        % Setting the y-axis limits
        y_max = max(Baseline_y_max);
        y_min = min(Baseline_y_min);
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        if ~isequal(per_dir_curs, 0)
            title_string = sprintf('Baseline Wrist Position, %i°, TgtCenter at %0.1f', ...
                target_dirs_noon(jj), target_centers_morn(jj));
            title(title_string, 'FontSize', title_font_size)
        else
            title_string = 'Baseline Wrist Position';
            title(title_string, 'FontSize', title_font_size)
        end

        % Get the top subplot title for saving
        if ~isequal(Save_Figs, 0)
            save_title(jj) = title_string;
        end

        % Annotation of the p_value
        if round(Baseline_CursorPos_p_values(1, jj), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(Baseline_CursorPos_p_values(1, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(Baseline_CursorPos_p_values(1, jj), 3), 0)
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

        % Annotation of the percent change
        if ~isequal(round(Baseline_CursorPos_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(Baseline_CursorPos_perc_changes(1, jj), 2)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(Baseline_CursorPos_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'xticklabel',{[]})
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes,'FontName', font_name);
        
        % Bottom Left Plot
        per_trial_ymax_morn = max(all_trials_z_cursor_p_morn, [], 'All');
        per_trial_ymax_noon = max(all_trials_z_cursor_p_noon, [], 'All');
        per_trial_ymin_morn = min(all_trials_z_cursor_p_morn, [], 'All');
        per_trial_ymin_noon = min(all_trials_z_cursor_p_noon, [], 'All');
        y_max = max(per_trial_ymax_morn, per_trial_ymax_noon);
        y_min = min(per_trial_ymin_morn, per_trial_ymin_noon);

        subplot(223) 
        hold on
        % Set the title
        title('Morning', 'FontSize', title_font_size);
        for pp = 1:width(all_trials_z_cursor_p_morn)
            plot(all_trials_z_cursor_p_morn(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.175 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_z_cursor_p_morn)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;

        % Setting the y-axis limits
        ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])

        subplot(224) % Bottom Right Plot
        hold on
        % Set the title
        title('Afternoon', 'FontSize', title_font_size);
        for pp = 1:width(all_trials_z_cursor_p_noon)
            plot(all_trials_z_cursor_p_noon(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.625 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_z_cursor_p_noon)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;

        % Setting the y-axis limits
        ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes, 'FontName', font_name);

    end % End of the Plot if-statement

    % End the function if you only want one baseline cursor
    if ~isequal(per_dir_curs, 1) && isequal(jj, 2)
        return
    end
    
end % End of target loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        save_title(ii) = strrep(save_title(ii), ':', '');
        save_title(ii) = strrep(save_title(ii), 'vs.', 'vs');
        save_title(ii) = strrep(save_title(ii), 'mg.', 'mg');
        save_title(ii) = strrep(save_title(ii), 'kg.', 'kg');
        save_title(ii) = strrep(save_title(ii), '.', '_');
        save_title(ii) = strrep(save_title(ii), '/', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title(ii))), 'fig')
        end
        close gcf
    end
end
