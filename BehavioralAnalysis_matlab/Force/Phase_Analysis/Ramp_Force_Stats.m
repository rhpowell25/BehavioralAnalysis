function [Ramp_Force_p_values, Ramp_Force_perc_changes, ...
    Morn_Ramp_Force, Err_Morn_Ramp_Force, Noon_Ramp_Force, Err_Noon_Ramp_Force] = ...
    Ramp_Force_Stats(xds_morn, xds_noon, norm_force, Plot_Figs, Save_File)

%% Display the function being used
disp('Ramp Force Statistics:');

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
        Morn_Ramp_Force = zeros(1, num_dirs);
        STD_Morn_TgtHold_Force = zeros(1, num_dirs);
        Err_Morn_Ramp_Force = zeros(1, num_dirs);
        Noon_Ramp_Force = zeros(1, num_dirs);
        STD_Noon_TgtHold_Force = zeros(1, num_dirs);
        Err_Noon_Ramp_Force = zeros(1, num_dirs);
        Ramp_Force_p_values = zeros(1, num_dirs);
        Ramp_Force_perc_changes = zeros(1, num_dirs);
    end

    %% Force and time aligned to specified event
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
    
    Force_morn = struct([]); % Force during each successful trial
    Force_timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_morn)
        Force_morn{ii, 1} = xds_morn.force(rewarded_start_idx_morn(ii) : rewarded_end_idx_morn(ii), :);
        Force_timing_morn{ii, 1} = xds_morn.time_frame(rewarded_start_idx_morn(ii) : rewarded_end_idx_morn(ii));
    end

    Force_noon = struct([]); % Force during each successful trial
    Force_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_noon)
        Force_noon{ii, 1} = xds_noon.force(rewarded_start_idx_noon(ii) : rewarded_end_idx_noon(ii), :);
        Force_timing_noon{ii, 1} = xds_noon.time_frame(rewarded_start_idx_noon(ii) : rewarded_end_idx_noon(ii));
    end

    %% Recompose the force
    [Sigma_Force_morn] = Sum_Force(xds_morn.meta.task, Force_morn);
    [Sigma_Force_noon] = Sum_Force(xds_noon.meta.task, Force_noon);

    %% Normalizing the average force

    if isequal(norm_force, 1)
        Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force);
        for ii = 1:length(Sigma_Force_morn)
            Sigma_Force_morn{ii,1} = (Sigma_Force_morn{ii,1} / Force_Norm_Factor)*100;
        end
        for ii = 1:length(Sigma_Force_noon)
            Sigma_Force_noon{ii,1} = (Sigma_Force_noon{ii,1} / Force_Norm_Factor)*100;
        end
    elseif strcmp(norm_force, 'Convert')
        for ii = 1:length(Sigma_Force_morn)
            Sigma_Force_morn{ii,1} = Sigma_Force_morn{ii,1} / 1000*5; % Millivolt conversion * gain
        end
        for ii = 1:length(Sigma_Force_noon)
            Sigma_Force_noon{ii,1} = Sigma_Force_noon{ii,1} / 1000*5; % Millivolt conversion * gain
        end
    end

    %% Calculating average Force (Average per trial)
    per_trial_avg_end_Force_morn = zeros(length(Force_morn), 1);
    per_trial_avg_end_Force_noon = zeros(length(Force_noon), 1);
    for ii = 1:length(Sigma_Force_morn)
        per_trial_avg_end_Force_morn(ii,1) = mean(Sigma_Force_morn{ii});
    end
    for ii = 1:length(Sigma_Force_noon)
        per_trial_avg_end_Force_noon(ii,1) = mean(Sigma_Force_noon{ii});
    end

    %% Calculate the output variables
    [~, Ramp_Force_p_values(1, jj)] = ttest2(per_trial_avg_end_Force_morn, per_trial_avg_end_Force_noon);
    % Average Force
    Morn_Ramp_Force(1, jj) = mean(per_trial_avg_end_Force_morn);
    Noon_Ramp_Force(1, jj) = mean(per_trial_avg_end_Force_noon);
    % Standard deviation
    STD_Morn_TgtHold_Force(1, jj) = std(per_trial_avg_end_Force_morn);
    STD_Noon_TgtHold_Force(1, jj) = std(per_trial_avg_end_Force_noon);
    % Standard error
    Err_Morn_Ramp_Force(1, jj) = STD_Morn_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_end_Force_morn));
    Err_Noon_Ramp_Force(1, jj) = STD_Noon_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_end_Force_noon));
    % Percent change
    Ramp_Force_perc_changes(1, jj) = (Noon_Ramp_Force - Morn_Ramp_Force) / abs(Morn_Ramp_Force);

    %% Plot the ramp phase Force

    if isequal(Plot_Figs, 1)

        %% Plot the Force box and whisker plot

        % Combining the morning and afternoon into one matrix
        avg_ramp_Force = cat(1, per_trial_avg_end_Force_morn, per_trial_avg_end_Force_noon);
        ramp_force_labels = zeros(length(avg_ramp_Force),1);
        ramp_force_labels(length(per_trial_avg_end_Force_morn) + 1:end) = 1;
        
        End_Force_mean(1,1) = Morn_Ramp_Force(1, jj);
        End_Force_mean(1,2) = Noon_Ramp_Force(1, jj);
        End_Force_STD(1,1) = STD_Morn_TgtHold_Force(1, jj);
        End_Force_STD(1,2) = STD_Noon_TgtHold_Force(1, jj);

        figure
        subplot(211) % Top Plot
        hold on
        boxplot(avg_ramp_Force, ramp_force_labels)

        % Setting the x-axis limits
        xlim([0.5, 2.5]);

        % Titling the plot
        Fig_Title = sprintf('Ramp Phase Force, %i°, TgtCenter at %0.1f', ... 
            target_dirs_noon(jj), target_centers_noon(jj));
        title(Fig_Title, 'FontSize', Plot_Params.title_font_size)

        % Annotation of the p_value
        if round(Ramp_Force_p_values(1, jj), 3) > 0
            legend_dims = [0 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(Ramp_Force_p_values(1, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(Ramp_Force_p_values(1, jj), 3), 0)
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
        if ~isequal(round(Ramp_Force_perc_changes(1, jj), 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(Ramp_Force_perc_changes(1, jj), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(Ramp_Force_perc_changes(1, jj), 3), 0)
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
        for pp = 1:length(Sigma_Force_morn)
            plot(Sigma_Force_morn{pp,1})
        end

        % Annotation of number of trials
        legend_dims = [0.175 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(length(Sigma_Force_morn)));
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
        for pp = 1:length(Sigma_Force_noon)
            plot(Sigma_Force_noon{pp,1})
        end

        % Annotation of number of trials
        legend_dims = [0.625 0.02 0.44 0.44];
        trial_num_string = strcat('n =', {' '}, mat2str(length(Sigma_Force_noon)));
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


