function [rxn_time] = Task_Metric_ViolinPlot(xds_morn, xds_noon, Task_Metric, Plot_Figs, Save_File)

%% File Description:

% This function plots a violin plot of reaction time (defined as the time after the
% go-cue when the EMG of interest exceeds 2 std of the baseline EMG), or the 
% trial length
% The EMG of interest is chosen based on the task / target.
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% Task_Metric: 'Rxn_Time', 'Trial_Length'
% Plot_Figs: 1 or 0
% Save_Figs: 'pdf', 'png', 'fig', 'All', or 0

%% Display the function being used
disp('Task Metric Violin Plot Function:')

%% Basic settings, some variable extractions, & definitions

% Do you want a vioin plot for each target / direction combo? (1 = Yes, 0 = No)
per_dir_plot = 0;

% How much do you want to expand the y-axis
axis_expansion = 0.25;

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

% Extract the target directions & centers
[target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
[target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

% Close all previously open figures if you're saving 
if ~isequal(Save_File, 0)
    close all
end

%% Check to see if both sessions use a consistent number of targets

% Find matching targets between the two sessions
[Matching_Idxs_Morn, Matching_Idxs_Noon] = ...
    Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon);

% Only use the info of target centers conserved between morn & noon
if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
    disp('Uneven Targets Between Morning & Afternoon');
    target_dirs_morn = target_dirs_morn(Matching_Idxs_Morn);
    target_dirs_noon = target_dirs_noon(Matching_Idxs_Noon);
    target_centers_morn = target_centers_morn(Matching_Idxs_Morn);
    target_centers_noon = target_centers_noon(Matching_Idxs_Noon);
end

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dir = length(target_dirs_morn);

%% Begin the loop through all directions
for jj = 1:num_dir

    %% Times for rewarded trials
    if strcmp(Task_Metric, 'Trial_Length')
        event = 'trial_end';
    elseif strcmp(Task_Metric, 'Rxn_Time')
        event = 'task_onset';
    end
    if isequal(per_dir_plot, 0)
        [rewarded_gocue_time_morn] = TrialAlignmentTimes(xds_morn, NaN, 'Max', 'trial_goCue');
        [rewarded_gocue_time_noon] = TrialAlignmentTimes(xds_noon, NaN, 'Max', 'trial_goCue');
        [Alignment_Times_morn] = EventAlignmentTimes(xds_morn, NaN, 'Max', event);
        [Alignment_Times_noon] = EventAlignmentTimes(xds_noon, NaN, 'Max', event);
    else
        [rewarded_gocue_time_morn] = ...
            TrialAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_goCue');
        [rewarded_gocue_time_noon] = ...
            TrialAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_goCue');
        [Alignment_Times_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), event);
        [Alignment_Times_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), event);
    end

    %% Define the output variables
    if jj == 1 && ~isequal(per_dir_plot, 0)
        rxn_time = zeros(num_dir, 1);
    else
        rxn_time = zeros(1, 1);
    end

    %% Find the difference between times
    time_length_morn = Alignment_Times_morn - rewarded_gocue_time_morn;
    time_length_noon = Alignment_Times_noon - rewarded_gocue_time_noon;

    %% Plot the violin plot

    if isequal(Plot_Figs, 1)

        violin_fig = figure;
        violin_fig.Position = [200 50 Plot_Params.fig_size Plot_Params.fig_size];
        
        % Title
        if ~isequal(per_dir_plot, 0)
            Fig_Title = strcat(Date, {' '}, Task, ',', {' '}, Drug, ':', {' '}, ...
                num2str(target_dirs_noon(jj)), '°,', {' '}, 'TgtCenter at', {' '}, ...
                num2str(target_centers_morn(jj)));
        else
            Fig_Title = strcat(Date, {' '}, Task, ',', {' '}, Drug);
        end
        sgtitle(Fig_Title, 'FontSize', Plot_Params.title_font_size);
    
        % Morning violin plot
        subplot(1,2,1);
        hold on
        violin_positions = (1:length(time_length_morn));
        Violin_Plot({time_length_morn}, violin_positions, 'ViolinColor', [0.9290, 0.6940, 0.1250]);
    
        % Annotation of the n-count
        morn_succ_trials = strcat('n =', {' '}, mat2str(length(time_length_morn)));
        morn_legend = annotation('textbox', [0.35 0.1 0.1 0.1], 'String', ...
            morn_succ_trials, 'FitBoxToText', 'on', 'EdgeColor','none', ...
            'VerticalAlignment', 'top', 'horizontalalignment', 'right');
        morn_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;
    
        % Find the axis limits for plotting
        y_min = min(cat(1, time_length_morn, time_length_noon)) - axis_expansion;
        y_max = max(cat(1, time_length_morn, time_length_noon)) + axis_expansion;
    
        % Axis limits
        ylim([y_min, y_max])
        xlim([0.5, 1.5]);
    
        % Labels
        xlabel('Morning', 'FontSize', Plot_Params.label_font_size)
        if strcmp(Task_Metric, 'Rxn_Time')
            ylabel('Reaction Time (Sec.)', 'FontSize', Plot_Params.label_font_size);
        end
        if strcmp(Task_Metric, 'Trial_Length')
            ylabel('Trial Length (Sec.)', 'FontSize', Plot_Params.label_font_size);
        end
    
        % Axis Editing
        figure_axes = gca;
        % Set ticks to outside
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set the tick label font size
        figure_axes.FontSize = Plot_Params.label_font_size;
    
        % Only label every other tick
        x_labels = string(figure_axes.XAxis.TickLabels);
        y_labels = string(figure_axes.YAxis.TickLabels);
        x_labels(1:end) = NaN;
        y_labels(2:2:end) = NaN;
        figure_axes.XAxis.TickLabels = x_labels;
        figure_axes.YAxis.TickLabels = y_labels;
    
        % Afternoon violin plot
        subplot(1,2,2);
        hold on
        violin_positions = (1:length(time_length_noon));
        Violin_Plot({time_length_noon}, violin_positions, 'ViolinColor', [.5 0 .5]);
    
        % Annotation of the n-count
        noon_succ_trials = strcat('n =', {' '}, mat2str(length(time_length_noon)));
        noon_legend = annotation('textbox', [0.6 0.1 0.1 0.1], 'String', ...
            noon_succ_trials, 'FitBoxToText', 'on', 'EdgeColor','none', ...
            'VerticalAlignment', 'top', 'horizontalalignment', 'right');
        noon_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;
    
        % Axis limits
        ylim([y_min, y_max])
        xlim([0.5, 1.5]);
    
        % Labels
        xlabel('Afternoon', 'FontSize', Plot_Params.label_font_size)
    
        % Axis Editing
        figure_axes = gca;
        % Set ticks to outside
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set the tick label font size
        figure_axes.FontSize = Plot_Params.label_font_size;
    
        % Only label every other tick
        x_labels = string(figure_axes.XAxis.TickLabels);
        y_labels = string(figure_axes.YAxis.TickLabels);
        x_labels(1:end) = NaN;
        y_labels(1:end) = NaN;
        figure_axes.XAxis.TickLabels = x_labels;
        figure_axes.YAxis.TickLabels = y_labels;

        % Do the statistics
        normality_morn = kstest(time_length_morn);
        normality_noon = kstest(time_length_noon);
        if normality_morn == 0 && normality_noon == 0
            disp('Distrubutions are normal:')
            [~, violin_plot_p_val] = ttest2(time_length_morn, time_length_noon);
        else
            disp('Distrubutions are not normal:')
            [violin_plot_p_val, ~] = ranksum(time_length_morn, time_length_noon);
        end
    
        % Annotation of the p-value
        if round(violin_plot_p_val, 3) > 0
            legend_dims = [0.015 0.45 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(violin_plot_p_val, 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
        end
        if isequal(round(violin_plot_p_val, 3), 0)
            legend_dims = [0.015 0.45 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
        end

        % Find the percent change
        avg_violin_plot_morn = mean(time_length_morn);
        avg_violin_plot_noon = mean(time_length_noon);
    
        violin_plot_perc_change = avg_violin_plot_noon - avg_violin_plot_morn / ...
                abs(avg_violin_plot_morn);
    
        % Annotation of the percent change
        if ~isequal(round(violin_plot_perc_change, 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(violin_plot_perc_change, 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end
        if isequal(round(violin_plot_perc_change, 3), 0)
            legend_dims = [0.55 0.45 0.44 0.44];
            perc_change_string = strcat('Δ% ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
                'FitBoxToText', 'on', 'verticalalignment', 'top', ...
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End of Plot_Figs statement

    %% Defining the output variables

    % Mean reaction time
    rxn_time(jj,1) = (mean(time_length_morn) + mean(time_length_noon)) / 2;

    

    %% End the function after one loop if using all targets
    if isequal(per_dir_plot, 0)
        break
    end

end




