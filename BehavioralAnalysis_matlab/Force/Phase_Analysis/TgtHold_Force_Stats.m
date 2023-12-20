function [TgtHold_Force_p_values, TgtHold_Force_perc_changes, ...
    Morn_TgtHold_Force, Err_Morn_TgtHold_Force, Noon_TgtHold_Force, Err_Noon_TgtHold_Force] = ...
    TgtHold_Force_Stats(xds_morn, xds_noon, norm_force, Plot_Figs, Save_File)

%% Display the function being used
disp('Target Hold Force Statistics:');

%% Ending the function if there is no force

if strcmp(xds_morn.meta.task, 'WS')
    disp('Event cannot be force related for this task');
    return
end

%% Basic Settings, some variable extractions, & definitions

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

% Save Counter
if ~isequal(Save_File, 0)
    close all
end

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

%% Begin the loop through all directions
for jj = 1:num_dirs

    %% Times for rewarded trials
    [rewarded_end_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_end');
    [rewarded_end_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_end');
    
    %% Define the output variables
    if jj == 1
        Morn_TgtHold_Force = zeros(1, num_dirs);
        STD_Morn_TgtHold_Force = zeros(1, num_dirs);
        Err_Morn_TgtHold_Force = zeros(1, num_dirs);
        Noon_TgtHold_Force = zeros(1, num_dirs);
        STD_Noon_TgtHold_Force = zeros(1, num_dirs);
        Err_Noon_TgtHold_Force = zeros(1, num_dirs);
        TgtHold_Force_p_values = zeros(1, num_dirs);
        TgtHold_Force_perc_changes = zeros(1, num_dirs);
    end

    %% Force and time aligned to specified event
    % Find the rewarded times in the whole trial time frame
    rewarded_end_idx_morn = zeros(height(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        rewarded_end_idx_morn(ii) = find(xds_morn.time_frame == rewarded_end_time_morn(ii));
    end

    rewarded_end_idx_noon = zeros(height(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        rewarded_end_idx_noon(ii) = find(xds_noon.time_frame == rewarded_end_time_noon(ii));
    end
    
    Force_morn = struct([]); % Force during each successful trial
    timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_morn)
        Force_morn{ii, 1} = xds_morn.force((rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width) : ...
            rewarded_end_idx_morn(ii)), :);
        timing_morn{ii, 1} = xds_morn.time_frame((rewarded_end_idx_morn(ii) - ... 
            (TgtHold_time / xds_morn.bin_width) : rewarded_end_idx_morn(ii)));
    end

    Force_noon = struct([]); % Force during each successful trial
    Force_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_noon)
        Force_noon{ii, 1} = xds_noon.force((rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width) : ...
            rewarded_end_idx_noon(ii)), :);
        Force_timing_noon{ii, 1} = xds_noon.time_frame((rewarded_end_idx_noon(ii) - ... 
            (TgtHold_time / xds_noon.bin_width) : rewarded_end_idx_noon(ii)));
    end
    
    % Finding the absolute timing
    absolute_end_Force_timing_morn = timing_morn{1,1} - rewarded_end_time_morn(1);
    absolute_end_Force_timing_noon = Force_timing_noon{1,1} - rewarded_end_time_noon(1);

    %% Recompose the force
    [Sigma_Force_morn] = Sum_Force(xds_morn.meta.task, Force_morn);
    [Sigma_Force_noon] = Sum_Force(xds_noon.meta.task, Force_noon);

    %% Putting all succesful trials in one array
    
    all_trials_end_Force_morn = zeros(length(Force_morn{1,1}), length(rewarded_end_time_morn));
    for ii = 1:length(rewarded_end_time_morn)
        all_trials_end_Force_morn(:,ii) = Sigma_Force_morn{ii, 1};
    end

    all_trials_end_Force_noon = zeros(length(Force_noon{1,1}), length(rewarded_end_time_noon));
    for ii = 1:length(rewarded_end_time_noon)
        all_trials_end_Force_noon(:,ii) = Sigma_Force_noon{ii, 1};
    end

    %% Normalizing the average force

    if isequal(norm_force, 1)
        Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force);
        all_trials_end_Force_morn = (all_trials_end_Force_morn / Force_Norm_Factor)*100;
        all_trials_end_Force_noon = (all_trials_end_Force_noon / Force_Norm_Factor)*100;
    elseif strcmp(norm_force, 'Convert')
        all_trials_end_Force_morn = all_trials_end_Force_morn / 1000*5; % Millivolt conversion * gain
        all_trials_end_Force_noon = all_trials_end_Force_noon / 1000*5; % Millivolt conversion * gain
    end

    %% Calculating average Force (Average per trial)
    per_trial_avg_end_Force_morn = zeros(length(Force_morn), 1);
    per_trial_avg_end_Force_noon = zeros(length(Force_noon), 1);
    for ii = 1:length(Force_morn)
        per_trial_avg_end_Force_morn(ii,1) = mean(all_trials_end_Force_morn(:,ii));
    end
    for ii = 1:length(Force_noon)
        per_trial_avg_end_Force_noon(ii,1) = mean(all_trials_end_Force_noon(:,ii));
    end
    
    %% Calculating average Force (Average Across trials)
    cross_trial_avg_end_Force_morn = zeros(length(Force_morn), 1);
    cross_trial_avg_end_Force_noon = zeros(length(Force_noon), 1);
    for ii = 1:length(Force_morn{1,1})
        cross_trial_avg_end_Force_morn(ii,1) = mean(all_trials_end_Force_morn(ii,:));
    end
    for ii = 1:length(Force_noon{1,1})
        cross_trial_avg_end_Force_noon(ii,1) = mean(all_trials_end_Force_noon(ii,:));
    end

    %% Calculate the output variables
    [~, TgtHold_Force_p_values(1, jj)] = ttest2(per_trial_avg_end_Force_morn, per_trial_avg_end_Force_noon);
    % Average Force
    Morn_TgtHold_Force(1, jj) = mean(per_trial_avg_end_Force_morn);
    Noon_TgtHold_Force(1, jj) = mean(per_trial_avg_end_Force_noon);
    % Standard deviation
    STD_Morn_TgtHold_Force(1, jj) = std(per_trial_avg_end_Force_morn);
    STD_Noon_TgtHold_Force(1, jj) = std(per_trial_avg_end_Force_noon);
    % Standard error
    Err_Morn_TgtHold_Force(1, jj) = STD_Morn_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_end_Force_morn));
    Err_Noon_TgtHold_Force(1, jj) = STD_Noon_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_end_Force_noon));
    % Percent change
    TgtHold_Force_perc_changes(1, jj) = (Noon_TgtHold_Force - Morn_TgtHold_Force) / abs(Morn_TgtHold_Force);

    %% Plot the target hold Force

    if isequal(Plot_Figs, 1)

        %% Plot the Force box and whisker plot

        % Combining the morning and afternoon into one matrix
        avg_end_Force(:,1) = cross_trial_avg_end_Force_morn;
        avg_end_Force(:,2) = cross_trial_avg_end_Force_noon;
        End_Force_mean(1,1) = Morn_TgtHold_Force(1, jj);
        End_Force_mean(1,2) = Noon_TgtHold_Force(1, jj);
        End_Force_STD(1,1) = STD_Morn_TgtHold_Force(1, jj);
        End_Force_STD(1,2) = STD_Noon_TgtHold_Force(1, jj);

        % Finding the min and max for the y-axis limit
        End_y_min = End_Force_mean - End_Force_STD;
        End_y_max = End_Force_mean + End_Force_STD;

        figure
        subplot(211) % Top Plot
        hold on
        boxplot(avg_end_Force)

        % Setting the y-axis limits
        y_max = max(End_y_max);
        y_min = min(End_y_min);
        ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])
        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        Fig_Title = sprintf('TgtHold Force, %i°, TgtCenter at %0.1f', ... 
            target_dirs_noon(jj), target_centers_noon(jj));
        title(Fig_Title, 'FontSize', Plot_Params.title_font_size)

        % Annotation of the p_value
        if round(TgtHold_Force_p_values(1, jj), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(TgtHold_Force_p_values(1, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(TgtHold_Force_p_values(1, jj), 3), 0)
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
        if ~isequal(round(TgtHold_Force_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(TgtHold_Force_perc_changes(1, jj), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(TgtHold_Force_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end

        per_trial_ymax_morn = max(all_trials_end_Force_morn, [], 'All');
        per_trial_ymax_noon = max(all_trials_end_Force_noon, [], 'All');
        per_trial_ymin_morn = min(all_trials_end_Force_morn, [], 'All');
        per_trial_ymin_noon = min(all_trials_end_Force_noon, [], 'All');
        y_max = max(per_trial_ymax_morn, per_trial_ymax_noon);
        y_min = min(per_trial_ymin_morn, per_trial_ymin_noon);

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
        for pp = 1:width(all_trials_end_Force_morn)
            plot(absolute_end_Force_timing_morn, all_trials_end_Force_morn(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.175 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_end_Force_morn)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;

        % Setting the axis limits
        ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])
        xlim([-TgtHold_time, 0])

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
        for pp = 1:width(all_trials_end_Force_noon)
            plot(absolute_end_Force_timing_noon, all_trials_end_Force_noon(:,pp))
        end

        % Annotation of number of trials
        legend_dims = [0.625 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_end_Force_noon)));
        legend_string = {char(trial_num_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;

        % Setting the axis limits
        ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])
        xlim([-TgtHold_time, 0])

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

