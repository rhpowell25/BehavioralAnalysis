function [TgtHold_CursorPos_p_values, TgtHold_CursorPos_perc_changes, ...
    Morn_TgtHold_CursorPos, Err_Morn_TgtHold_CursorPos, Noon_TgtHold_CursorPos, Err_Noon_TgtHold_CursorPos] = ...
    TgtHold_CursorPos_Stats(xds_morn, xds_noon, norm_cursor, Plot_Figs, Save_File)

%% File Description:

% This function finds changes in cursor position during the target hold 
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
disp('TgtHold Cursor Position Statistics:');

%% Basic Settings, some variable extractions, & definitions

% Font specifications
legend_font_size = 12;
title_font_size = 15;
font_name = 'Arial';

% Define the window for the baseline phase
TgtHold_time = xds_morn.meta.TgtHold;

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
if ~isequal(Save_File, 0)
    close all
end

%% Begin the loop through all directions
for jj = 1:num_dirs

    %% Times for rewarded trials
    [rewarded_end_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_end');
    [rewarded_end_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_end');
    
    %% Define the output variables
    if jj == 1
        Morn_TgtHold_CursorPos = zeros(1, num_dirs);
        Err_Morn_TgtHold_CursorPos = zeros(1, num_dirs);
        Noon_TgtHold_CursorPos = zeros(1, num_dirs);
        Err_Noon_TgtHold_CursorPos = zeros(1, num_dirs);
        TgtHold_CursorPos_p_values = zeros(1, num_dirs);
        TgtHold_CursorPos_perc_changes = zeros(1, num_dirs);
    end

    %% Cursor position aligned to trial end
    
    idx_length = TgtHold_time / xds_morn.bin_width;

    % Cursor position during each succesful trial
    cursor_p_morn = struct([]); % Cursor position at the start
    for ii = 1:length(rewarded_end_time_morn)
        idx_morn = find(xds_morn.time_frame == rewarded_end_time_morn(ii));
        cursor_p_morn{ii, 1} = xds_morn.curs_p(idx_morn - idx_length : idx_morn, :);
    end
    
    cursor_p_noon = struct([]); % Cursor position at the start
    for ii = 1:length(rewarded_end_time_noon)
        idx_noon = find(xds_noon.time_frame == rewarded_end_time_noon(ii));
        cursor_p_noon{ii, 1} = xds_noon.curs_p(idx_noon - idx_length : idx_noon, :);
    end

    %% Find the vector sum of the cursor position
    z_cursor_p_morn = struct([]);
    for ii = 1:length(rewarded_end_time_morn)
        z_cursor_p_morn{ii,1} = sqrt(cursor_p_morn{ii,1}(:, 2).^2 + cursor_p_morn{ii, 1}(:, 1).^2);
    end

    z_cursor_p_noon = struct([]);
    for ii = 1:length(rewarded_end_time_noon)
        z_cursor_p_noon{ii,1} = sqrt(cursor_p_noon{ii,1}(:, 2).^2 + cursor_p_noon{ii, 1}(:, 1).^2);
    end
        
    %% Putting all succesful trials in one array
    all_trials_z_cursor_p_morn = zeros(length(z_cursor_p_morn{1,1}), length(rewarded_end_time_morn));
    for ii = 1:length(rewarded_end_time_morn)
        all_trials_z_cursor_p_morn(:,ii) = z_cursor_p_morn{ii, 1};
    end

    all_trials_z_cursor_p_noon = zeros(length(z_cursor_p_noon{1,1}), length(rewarded_end_time_noon));
    for ii = 1:length(rewarded_end_time_noon)
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

    %% Calculate the output variables
    [~, TgtHold_CursorPos_p_values(1, jj)] = ...
        ttest2(per_trial_avg_z_cursor_p_morn, per_trial_avg_z_cursor_p_noon);
    % Average Cursor Position
    Morn_TgtHold_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_morn);
    Noon_TgtHold_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_noon);
    % Standard deviation
    End_cursor_p_morn_std = std(per_trial_avg_z_cursor_p_morn);
    End_cursor_p_noon_std = std(per_trial_avg_z_cursor_p_noon);
    % Standard error
    Err_Morn_TgtHold_CursorPos(1, jj) = End_cursor_p_morn_std / sqrt(length(per_trial_avg_z_cursor_p_morn));
    Err_Noon_TgtHold_CursorPos(1, jj) = End_cursor_p_noon_std / sqrt(length(per_trial_avg_z_cursor_p_noon));
    % Cursor position percent change
    TgtHold_CursorPos_perc_changes(1, jj) = (Noon_TgtHold_CursorPos - Morn_TgtHold_CursorPos) / abs(Morn_TgtHold_CursorPos);

    %% Plot the target hold cursor position

    if isequal(Plot_Figs, 1)

        %% Plot the EMG box and whisker plot

        % Combining the morning & afternoon into one matrix
        avg_cursor_p = [per_trial_avg_z_cursor_p_morn; per_trial_avg_z_cursor_p_noon];
        TgtHold_cursor_p_labels = zeros(length(avg_cursor_p),1);
        TgtHold_cursor_p_labels(length(per_trial_avg_z_cursor_p_morn) + 1:end) = 1;

        End_cursor_p_mean(1, 1) = Morn_TgtHold_CursorPos(1, jj);
        End_cursor_p_mean(1, 2) = Noon_TgtHold_CursorPos(1, jj);
        End_cursor_p_std(1, 1) = End_cursor_p_morn_std;
        End_cursor_p_std(1, 2) = End_cursor_p_noon_std;

        % Finding the min and max for the y-axis limit
        End_y_min = End_cursor_p_mean - End_cursor_p_std;
        End_y_max = End_cursor_p_mean + End_cursor_p_std;

        figure
        subplot(211) % Top Plot
        hold on
        boxplot(avg_cursor_p, TgtHold_cursor_p_labels)
    
        % Setting the y-axis limits
        y_max = max(End_y_max);
        y_min = min(End_y_min);
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        Fig_Title = sprintf('TgtHold Wrist Position, %i°, TgtCenter at %0.1f', ...
            target_dirs_noon(jj), target_centers_noon(jj));
        title(Fig_Title, 'FontSize', title_font_size)

        % Annotation of the p_value
        if round(TgtHold_CursorPos_p_values(1, jj), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(TgtHold_CursorPos_p_values(1, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(TgtHold_CursorPos_p_values(1, jj), 3), 0)
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
        if ~isequal(round(TgtHold_CursorPos_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(TgtHold_CursorPos_perc_changes(1, jj), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(TgtHold_CursorPos_perc_changes(1, jj), 3), 0)
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

        % Setting the axis limits
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        xlim([0, idx_length])

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

        % Setting the axis limits
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        xlim([0, idx_length])

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes, 'FontName', font_name);

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End of the Plot if-statement

end % End of target loop
