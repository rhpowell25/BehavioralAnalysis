function [Ramp_CursorPos_p_values, Ramp_CursorPos_perc_changes, ...
    Morn_Ramp_CursorPos, Err_Morn_Ramp_CursorPos, Noon_Ramp_CursorPos, Err_Noon_Ramp_CursorPos] = ...
    Ramp_CursorPos_Stats(xds_morn, xds_noon, norm_cursor, Plot_Figs, Save_File)

%% File Description:

% This function finds changes in cursor position during the ramp phase 
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
disp('Ramp Cursor Position Statistics:');

%% Basic Settings, some variable extractions, & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

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
    [rewarded_onset_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'task_onset');
    [rewarded_end_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_end');
    rewarded_end_time_morn = rewarded_end_time_morn - TgtHold_time;
    % Round to match the neural bin size
    rewarded_end_time_morn = round(rewarded_end_time_morn, abs(floor(log10(xds_morn.bin_width))));

    [rewarded_onset_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'task_onset');
    [rewarded_end_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_end');
    rewarded_end_time_noon = rewarded_end_time_noon - TgtHold_time;
    % Round to match the neural bin size
    rewarded_end_time_noon = round(rewarded_end_time_noon, abs(floor(log10(xds_noon.bin_width))));
    
    %% Define the output variables
    if jj == 1
        Morn_Ramp_CursorPos = zeros(1, num_dirs);
        STD_Morn_TgtHold_cursor_p = zeros(1, num_dirs);
        Err_Morn_Ramp_CursorPos = zeros(1, num_dirs);
        Noon_Ramp_CursorPos = zeros(1, num_dirs);
        STD_Noon_TgtHold_cursor_p = zeros(1, num_dirs);
        Err_Noon_Ramp_CursorPos = zeros(1, num_dirs);
        Ramp_CursorPos_p_values = zeros(1, num_dirs);
        Ramp_CursorPos_perc_changes = zeros(1, num_dirs);
    end

    %% Cursor position aligned to specified event
    % Find the rewarded times in the whole trial time frame
    rewarded_start_idx_morn = zeros(height(rewarded_onset_time_morn),1);
    for ii = 1:length(rewarded_onset_time_morn)
        rewarded_start_idx_morn(ii) = find(xds_morn.time_frame == rewarded_onset_time_morn(ii));
    end

    rewarded_end_idx_morn = zeros(height(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        rewarded_end_idx_morn(ii) = find(xds_morn.time_frame == rewarded_end_time_morn(ii));
    end

    rewarded_start_idx_noon = zeros(height(rewarded_onset_time_noon),1);
    for ii = 1:length(rewarded_onset_time_noon)
        rewarded_start_idx_noon(ii) = find(xds_noon.time_frame == rewarded_onset_time_noon(ii));
    end

    rewarded_end_idx_noon = zeros(height(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        rewarded_end_idx_noon(ii) = find(xds_noon.time_frame == rewarded_end_time_noon(ii));
    end
    
    % Cursor position during each succesful trial
    cursor_p_morn = struct([]);
    for ii = 1:length(rewarded_end_time_morn)
        cursor_p_morn{ii, 1} = xds_morn.curs_p(rewarded_start_idx_morn(ii) : rewarded_end_idx_morn(ii), :);
    end

    % Cursor position during each succesful trial
    cursor_p_noon = struct([]);
    for ii = 1:length(rewarded_end_time_noon)
        cursor_p_noon{ii, 1} = xds_noon.curs_p(rewarded_start_idx_noon(ii) : rewarded_end_idx_noon(ii), :);
    end

    %% Find the vector sum of the cursor position
    z_cursor_p_morn = struct([]);
    for ii = 1:length(rewarded_onset_time_morn)
        z_cursor_p_morn{ii,1} = sqrt(cursor_p_morn{ii,1}(:, 2).^2 + cursor_p_morn{ii, 1}(:, 1).^2);
    end

    z_cursor_p_noon = struct([]);
    for ii = 1:length(rewarded_onset_time_noon)
        z_cursor_p_noon{ii,1} = sqrt(cursor_p_noon{ii,1}(:, 2).^2 + cursor_p_noon{ii, 1}(:, 1).^2);
    end

    %% Normalizing the average cursor position

    Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_cursor);
    for ii = 1:length(z_cursor_p_morn)
        z_cursor_p_morn{ii,1} = (z_cursor_p_morn{ii,1} / Cursor_Norm_Factor)*100;
    end
    for ii = 1:length(z_cursor_p_noon)
        z_cursor_p_noon{ii,1} = (z_cursor_p_noon{ii,1} / Cursor_Norm_Factor)*100;
    end

    %% Calculating average cursor position (Average per trial)
    per_trial_avg_z_cursor_p_morn = zeros(length(cursor_p_morn), 1);
    per_trial_avg_z_cursor_p_noon = zeros(length(cursor_p_noon), 1);
    for ii = 1:length(z_cursor_p_morn)
        per_trial_avg_z_cursor_p_morn(ii,1) = mean(z_cursor_p_morn{ii});
    end
    for ii = 1:length(z_cursor_p_noon)
        per_trial_avg_z_cursor_p_noon(ii,1) = mean(z_cursor_p_noon{ii});
    end

    %% Do the statistics on the cursor position
    [~, Ramp_CursorPos_p_values(1, jj)] = ttest2(per_trial_avg_z_cursor_p_morn, per_trial_avg_z_cursor_p_noon);
    % Average cursor position
    Morn_Ramp_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_morn, 'omitnan');
    Noon_Ramp_CursorPos(1, jj) = mean(per_trial_avg_z_cursor_p_noon, 'omitnan');
    % Standard deviation
    STD_Morn_TgtHold_cursor_p(1, jj) = std(per_trial_avg_z_cursor_p_morn);
    STD_Noon_TgtHold_cursor_p(1, jj) = std(per_trial_avg_z_cursor_p_noon);
    % Standard error
    Err_Morn_Ramp_CursorPos(1, jj) = STD_Morn_TgtHold_cursor_p(1, jj) / sqrt(length(per_trial_avg_z_cursor_p_morn));
    Err_Noon_Ramp_CursorPos(1, jj) = STD_Noon_TgtHold_cursor_p(1, jj) / sqrt(length(per_trial_avg_z_cursor_p_noon));
    % Percent change
    Ramp_CursorPos_perc_changes(1, jj) = (Noon_Ramp_CursorPos - Morn_Ramp_CursorPos) / abs(Morn_Ramp_CursorPos);

    %% Plot the cursor position box & whisker plot
    if isequal(Plot_Figs, 1)

        % Top Plot
        % Combining the morning & afternoon into one matrix
        avg_cursor_p = [per_trial_avg_z_cursor_p_morn; per_trial_avg_z_cursor_p_noon];
        ramp_cursor_p_labels = zeros(length(avg_cursor_p),1);
        ramp_cursor_p_labels(length(per_trial_avg_z_cursor_p_morn) + 1:end) = 1;
        
        End_cursor_p_mean(1,1) = Morn_Ramp_CursorPos(1, jj);
        End_cursor_p_mean(1,2) = Noon_Ramp_CursorPos(1, jj);
        End_cursor_p_STD(1,1) = STD_Morn_TgtHold_cursor_p(1, jj);
        End_cursor_p_STD(1,2) = STD_Noon_TgtHold_cursor_p(1, jj);

        figure
        subplot(211) % Top Plot
        hold on
        boxplot(avg_cursor_p, ramp_cursor_p_labels)

        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        Fig_Title = sprintf('Ramp Wrist Position, %i°, TgtCenter at %0.1f', ... 
            target_dirs_noon(jj), target_centers_noon(jj));
        title(Fig_Title, 'FontSize', Plot_Params.title_font_size)
        
        % Annotation of the p_value
        if round(Ramp_CursorPos_p_values(1, jj), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(Ramp_CursorPos_p_values(1, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(Ramp_CursorPos_p_values(1, jj), 3), 0)
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end

        % Annotation of the percent change
        if ~isequal(round(Ramp_CursorPos_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(Ramp_CursorPos_perc_changes(1, jj), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(Ramp_CursorPos_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'xticklabel',{[]})
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes,'FontName', Plot_Params.font_name);

        subplot(223) % Bottom Left Plot
        hold on
        % Set the title
        title('Morning', 'FontSize', Plot_Params.title_font_size);
        for pp = 1:length(z_cursor_p_morn)
            plot(z_cursor_p_morn{pp,1})
        end

        % Annotation of number of trials
        legend_dims = [0.175 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(length(z_cursor_p_morn)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes,'FontName', Plot_Params.font_name);

        subplot(224) % Bottom Right Plot
        hold on
        % Set the title
        title('Afternoon', 'FontSize', Plot_Params.title_font_size);
        for pp = 1:length(z_cursor_p_noon)
            plot(z_cursor_p_noon{pp,1})
        end

        % Annotation of number of trials
        legend_dims = [0.625 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(length(z_cursor_p_noon)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;

        % Set ticks to outside
        figure_axes = gca;
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set The Font
        set(figure_axes,'FontName', Plot_Params.font_name);

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End of the Plot if-statement

end % End of target loop
