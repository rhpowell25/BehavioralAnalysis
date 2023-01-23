function [TgtHold_CursorPos_p_values, TgtHold_CursorPos_perc_changes, ...
    Morn_TgtHold_CursorPos, Err_Morn_TgtHold_CursorPos, Noon_TgtHold_CursorPos, Err_Noon_TgtHold_CursorPos] = ...
    TgtHold_CursorPos_Stats(xds_morn, xds_noon, Cursor_Norm_Factor, Plot_Figs, Save_Figs)

%% Display the function being used
disp('TgtHold Cursor Position Statistics:');

%% Basic Settings, some variable extractions, & definitions

% Font specifications
legend_font_size = 12;
title_font_size = 15;
font_name = 'Arial';

% Save Counter
ss = 1;
if ~isequal(Save_Figs, 0)
    close all
    save_title = strings;
end

% Define the window for the baseline phase
tgt_hold_time = xds_morn.meta.TgtHold;

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

    %% Cursor position and time aligned to trial end
    % Find the rewarded times in the whole trial time frame
    rewarded_end_idx_morn = zeros(height(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        rewarded_end_idx_morn(ii) = find(xds_morn.time_frame == rewarded_end_time_morn(ii));
    end

    rewarded_end_idx_noon = zeros(height(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        rewarded_end_idx_noon(ii) = find(xds_noon.time_frame == rewarded_end_time_noon(ii));
    end

    aligned_end_cursor_p_morn = struct([]); % Cursor position at the end
    aligned_end_timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_morn)
        aligned_end_cursor_p_morn{ii, 1} = xds_morn.curs_p((rewarded_end_idx_morn(ii) - ...
            (tgt_hold_time / xds_morn.bin_width) : rewarded_end_idx_morn(ii)), :);
        aligned_end_timing_morn{ii, 1} = xds_morn.time_frame((rewarded_end_idx_morn(ii) - ... 
            (tgt_hold_time / xds_morn.bin_width) : rewarded_end_idx_morn(ii)));
    end
    
    aligned_end_cursor_p_noon = struct([]); % Cursor position at the end
    aligned_end_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_noon)
        aligned_end_cursor_p_noon{ii, 1} = xds_noon.curs_p((rewarded_end_idx_noon(ii) - ...
            (tgt_hold_time / xds_noon.bin_width) : rewarded_end_idx_noon(ii)), :);
        aligned_end_timing_noon{ii, 1} = xds_noon.time_frame((rewarded_end_idx_noon(ii) - ... 
            (tgt_hold_time / xds_noon.bin_width) : rewarded_end_idx_noon(ii)));
    end

    % Finding the absolute timing
    absolute_end_cursor_p_timing_morn = aligned_end_timing_morn{1,1} - rewarded_end_time_morn(1);
    absolute_end_cursor_p_timing_noon = aligned_end_timing_noon{1,1} - rewarded_end_time_noon(1);

    %% Recompose the cursor positions
    z_cursor_p_morn = struct([]);
    for ii = 1:length(rewarded_end_time_morn)
        z_cursor_p_morn{ii,1} = sqrt(aligned_end_cursor_p_morn{ii,1}(:, 2).^2 + aligned_end_cursor_p_morn{ii, 1}(:, 1).^2);
    end

    z_cursor_p_noon = struct([]);
    for ii = 1:length(rewarded_end_time_noon)
        z_cursor_p_noon{ii,1} = sqrt(aligned_end_cursor_p_noon{ii,1}(:, 2).^2 + aligned_end_cursor_p_noon{ii, 1}(:, 1).^2);
    end
        
    %% Putting all succesful trials in one array
    all_trials_end_cursor_p_morn = zeros(length(z_cursor_p_morn{1,1}), length(rewarded_end_time_morn));
    for ii = 1:length(rewarded_end_time_morn)
        all_trials_end_cursor_p_morn(:,ii) = z_cursor_p_morn{ii, 1};
    end

    all_trials_end_cursor_p_noon = zeros(length(z_cursor_p_noon{1,1}), length(rewarded_end_time_noon));
    for ii = 1:length(rewarded_end_time_noon)
        all_trials_end_cursor_p_noon(:,ii) = z_cursor_p_noon{ii, 1};
    end

    %% Normalizing the average cursor position
    all_trials_end_cursor_p_morn = all_trials_end_cursor_p_morn / Cursor_Norm_Factor*100;
    all_trials_end_cursor_p_noon = all_trials_end_cursor_p_noon / Cursor_Norm_Factor*100;

    %% Calculating average cursor position (Average per trial)
    per_trial_avg_end_cursor_p_morn = zeros(length(aligned_end_cursor_p_morn), 1);
    per_trial_avg_end_cursor_p_noon = zeros(length(aligned_end_cursor_p_noon), 1);
    for ii = 1:length(aligned_end_cursor_p_morn)
        per_trial_avg_end_cursor_p_morn(ii, 1) = mean(all_trials_end_cursor_p_morn(:, ii));
    end
    for ii = 1:length(aligned_end_cursor_p_noon)
        per_trial_avg_end_cursor_p_noon(ii, 1) = mean(all_trials_end_cursor_p_noon(:, ii));
    end

    %% Calculate the output variables
    [~, TgtHold_CursorPos_p_values(1, ss)] = ...
        ttest2(per_trial_avg_end_cursor_p_morn, per_trial_avg_end_cursor_p_noon);
    % Average Cursor Position
    Morn_TgtHold_CursorPos(1, ss) = mean(per_trial_avg_end_cursor_p_morn);
    Noon_TgtHold_CursorPos(1, ss) = mean(per_trial_avg_end_cursor_p_noon);
    % Standard deviation
    End_cursor_p_morn_std = std(per_trial_avg_end_cursor_p_morn);
    End_cursor_p_noon_std = std(per_trial_avg_end_cursor_p_noon);
    % Standard error
    Err_Morn_TgtHold_CursorPos(1, ss) = End_cursor_p_morn_std / sqrt(length(per_trial_avg_end_cursor_p_morn));
    Err_Noon_TgtHold_CursorPos(1, ss) = End_cursor_p_noon_std / sqrt(length(per_trial_avg_end_cursor_p_noon));
    % Cursor position percent change
    TgtHold_CursorPos_perc_changes(1, ss) = (Noon_TgtHold_CursorPos - Morn_TgtHold_CursorPos) / abs(Morn_TgtHold_CursorPos);

    %% Plot the target hold cursor position

    if isequal(Plot_Figs, 1)

        %% Plot the EMG box and whisker plot

        % Combining the morning and afternoon into one matrix
        avg_TgtHold_cursor_p =  [per_trial_avg_end_cursor_p_morn; per_trial_avg_end_cursor_p_noon];
        % Put in the categories for the x-axis
        categ_morn = cell(length(per_trial_avg_end_cursor_p_morn),1);
        categ_noon = cell(length(per_trial_avg_end_cursor_p_noon),1);
        for cc = 1:length(categ_morn)
            categ_morn{cc} = 'Morning';
        end
        for cc = 1:length(categ_noon)
            categ_noon{cc} = 'Afternoon';
        end
        categ = [categ_morn; categ_noon];

        End_cursor_p_mean(1, 1) = Morn_TgtHold_CursorPos(1, ss);
        End_cursor_p_mean(1, 2) = Noon_TgtHold_CursorPos(1, ss);
        End_cursor_p_std(1, 1) = End_cursor_p_morn_std;
        End_cursor_p_std(1, 2) = End_cursor_p_noon_std;

        % Finding the min and max for the y-axis limit
        End_y_min = End_cursor_p_mean - End_cursor_p_std;
        End_y_max = End_cursor_p_mean + End_cursor_p_std;

        figure
        subplot(211) % Top Plot
        hold on
        boxplot(avg_TgtHold_cursor_p, categ)
    
        % Setting the y-axis limits
        y_max = max(End_y_max);
        y_min = min(End_y_min);
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        title(sprintf('Trial End Wrist Position, %i°, TgtCenter at %0.1f', ...
            target_dirs_noon(jj), target_centers_noon(jj)), 'FontSize', title_font_size)

        % Get the top subplot title for saving
        if ~isequal(Save_Figs, 0)
            fig_info = get(gca,'title');
            save_title(ss) = get(fig_info, 'string');
        end

        % Annotation of the p_value
        if round(TgtHold_CursorPos_p_values(1, ss), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(TgtHold_CursorPos_p_values(1, ss), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

        if isequal(round(TgtHold_CursorPos_p_values(1, ss), 3), 0)
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
        if ~isequal(round(TgtHold_CursorPos_perc_changes(1, ss), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('PC =', {' '}, mat2str(round(TgtHold_CursorPos_perc_changes(1, ss), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(TgtHold_CursorPos_perc_changes(1, ss), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('PC ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

        per_trial_ymax_morn = max(all_trials_end_cursor_p_morn, [], 'All');
        per_trial_ymax_noon = max(all_trials_end_cursor_p_noon, [], 'All');
        per_trial_ymin_morn = min(all_trials_end_cursor_p_morn, [], 'All');
        per_trial_ymin_noon = min(all_trials_end_cursor_p_noon, [], 'All');
        y_max = max(per_trial_ymax_morn, per_trial_ymax_noon);
        y_min = min(per_trial_ymin_morn, per_trial_ymin_noon);

        subplot(223) % Bottom Left Plot
        hold on
        % Set the title
        title('Morning', 'FontSize', title_font_size);
        for pp = 1:width(all_trials_end_cursor_p_morn)
            plot(absolute_end_cursor_p_timing_morn, all_trials_end_cursor_p_morn(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.175 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_end_cursor_p_morn)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;

        % Setting the y-axis limits
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        xlim([-tgt_hold_time, 0])

        subplot(224) % Bottom Right Plot
        hold on
        % Set the title
        title('Afternoon', 'FontSize', title_font_size);
        for pp = 1:width(all_trials_end_cursor_p_noon)
            plot(absolute_end_cursor_p_timing_noon, all_trials_end_cursor_p_noon(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.625 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_end_cursor_p_noon)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;

        % Setting the y-axis limits
        ylim([y_min - 1.5*abs(y_min), y_max + 1.5*abs(y_max)])
        xlim([-tgt_hold_time, 0])

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes, 'FontName', font_name);

    end % End of the Plot if-statement

    % Add to the loop counter
    ss = ss + 1;

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
