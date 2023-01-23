function [EMG_Names, Baseline_EMG_p_values, Baseline_EMG_perc_changes, ...
    Morn_Baseline_EMG, Err_Morn_Baseline_EMG, Noon_Baseline_EMG, Err_Noon_Baseline_EMG] = ...
    Baseline_EMG_Stats(xds_morn, xds_noon, muscle_groups, EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, Save_Figs)

%% Display the function being used
disp('Baseline EMG Statistics:');

%% Find the EMG index

[M] = EMG_Index(xds_morn, muscle_groups);

EMG_Names = strings;
for ii = 1:length(M)
    EMG_Names(ii,1) = xds_morn.EMG_names(M(ii));
end

%% Basic Settings, some variable extractions, & definitions

% Font specifications
legend_font_size = 12;
title_font_size = 15;
font_name = 'Arial';

if ~isequal(Save_Figs, 0)
    save_title = strings;
end

% Define the window for the baseline phase
time_before_gocue = 0.4;

% Do you want a baseline EMG for each target / direction combo? (1 = Yes, 0 = No)
per_dir_EMG = 0;

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
    if ~isequal(per_dir_EMG, 0)
        [rewarded_gocue_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_gocue');
        [rewarded_gocue_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_gocue');
    else
        [rewarded_gocue_time_morn] = EventAlignmentTimes(xds_morn, NaN, NaN, 'trial_gocue');
        [rewarded_gocue_time_noon] = EventAlignmentTimes(xds_noon, NaN, NaN, 'trial_gocue');
    end

    %% Define the output variables
    if jj == 1 && ~isequal(per_dir_EMG, 0)
        Morn_Baseline_EMG = zeros(length(M), num_dirs);
        STD_Morn_Baseline_EMG = zeros(length(M), num_dirs);
        Err_Morn_Baseline_EMG = zeros(length(M), num_dirs);
        Noon_Baseline_EMG = zeros(length(M), num_dirs);
        STD_Noon_Baseline_EMG = zeros(length(M), num_dirs);
        Err_Noon_Baseline_EMG = zeros(length(M), num_dirs);
        Baseline_EMG_p_values = zeros(length(M), num_dirs);
        Baseline_EMG_perc_changes = zeros(length(M), num_dirs);
    elseif jj == 1 && ~isequal(per_dir_EMG, 1)
        Morn_Baseline_EMG = zeros(length(M), 1);
        STD_Morn_Baseline_EMG = zeros(length(M), 1);
        Err_Morn_Baseline_EMG = zeros(length(M), 1);
        Noon_Baseline_EMG = zeros(length(M), 1);
        STD_Noon_Baseline_EMG = zeros(length(M), 1);
        Err_Noon_Baseline_EMG = zeros(length(M), 1);
        Baseline_EMG_p_values = zeros(length(M), 1);
        Baseline_EMG_perc_changes = zeros(length(M), 1);
    end

    %% EMG and time aligned to specified event
    % Find the rewarded times in the whole trial time frame
    rewarded_gocue_idx_morn = zeros(height(rewarded_gocue_time_morn),1);
    for ii = 1:length(rewarded_gocue_time_morn)
        rewarded_gocue_idx_morn(ii) = find(xds_morn.time_frame == rewarded_gocue_time_morn(ii));
    end

    rewarded_gocue_idx_noon = zeros(height(rewarded_gocue_time_noon),1);
    for ii = 1:length(rewarded_gocue_time_noon)
        rewarded_gocue_idx_noon(ii) = find(xds_noon.time_frame == rewarded_gocue_time_noon(ii));
    end
    
    aligned_baseline_EMG_morn = struct([]); % EMG during each successful trial
    aligned_baseline_EMG_timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_gocue_time_morn)
        aligned_baseline_EMG_morn{ii, 1} = xds_morn.EMG((rewarded_gocue_idx_morn(ii) - (time_before_gocue / xds_morn.bin_width) : ...
            rewarded_gocue_idx_morn(ii)), :);
        aligned_baseline_EMG_timing_morn{ii, 1} = xds_morn.time_frame((rewarded_gocue_idx_morn(ii) - ... 
            (time_before_gocue / xds_morn.bin_width) : rewarded_gocue_idx_morn(ii)));
    end

    aligned_baseline_EMG_noon = struct([]); % EMG during each successful trial
    aligned_baseline_EMG_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_gocue_time_noon)
        aligned_baseline_EMG_noon{ii, 1} = xds_noon.EMG((rewarded_gocue_idx_noon(ii) - (time_before_gocue / xds_noon.bin_width) : ...
            rewarded_gocue_idx_noon(ii)), :);
        aligned_baseline_EMG_timing_noon{ii, 1} = xds_noon.time_frame((rewarded_gocue_idx_noon(ii) - ... 
            (time_before_gocue / xds_noon.bin_width) : rewarded_gocue_idx_noon(ii)));
    end
    
    % Finding the absolute timing
    absolute_baseline_EMG_timing_morn = aligned_baseline_EMG_timing_morn{1,1} - rewarded_gocue_time_morn(1);
    absolute_baseline_EMG_timing_noon = aligned_baseline_EMG_timing_noon{1,1} - rewarded_gocue_time_noon(1);

    %% Putting all succesful trials in one array
    all_trials_baseline_EMG_morn = struct([]);
    for ii = 1:length(M)
        all_trials_baseline_EMG_morn{ii,1} = zeros(length(aligned_baseline_EMG_morn{1,1}),length(rewarded_gocue_time_morn));
        for mm = 1:length(rewarded_gocue_time_morn)
            all_trials_baseline_EMG_morn{ii,1}(:,mm) = aligned_baseline_EMG_morn{mm, 1}(:, M(ii));
        end
    end

    all_trials_baseline_EMG_noon = struct([]);
    for ii = 1:length(M)
        all_trials_baseline_EMG_noon{ii,1} = zeros(length(aligned_baseline_EMG_noon{1,1}),length(rewarded_gocue_time_noon));
        for mm = 1:length(rewarded_gocue_time_noon)
            all_trials_baseline_EMG_noon{ii,1}(:,mm) = aligned_baseline_EMG_noon{mm, 1}(:, M(ii));
        end
    end

    %% Zeroing the EMG
    for ii = 1:length(M)
        all_trials_baseline_EMG_morn{ii,1} = all_trials_baseline_EMG_morn{ii,1} - EMG_Zero_Factor(ii);
        all_trials_baseline_EMG_noon{ii,1} = all_trials_baseline_EMG_noon{ii,1} - EMG_Zero_Factor(ii);
    end

    %% Normalizing the average EMG's
    for ii = 1:length(M)
        all_trials_baseline_EMG_morn{ii,1} = (all_trials_baseline_EMG_morn{ii,1} / EMG_Norm_Factor(ii))*100;
        all_trials_baseline_EMG_noon{ii,1} = (all_trials_baseline_EMG_noon{ii,1} / EMG_Norm_Factor(ii))*100;
    end

    %% Calculating average EMG (Average per trial)
    per_trial_avg_baseline_EMG_morn = zeros(length(aligned_baseline_EMG_morn),length(M));
    per_trial_avg_baseline_EMG_noon = zeros(length(aligned_baseline_EMG_noon),length(M));
    for ii = 1:length(M)
        for mm = 1:length(aligned_baseline_EMG_morn)
            per_trial_avg_baseline_EMG_morn(mm,ii) = mean(all_trials_baseline_EMG_morn{ii,1}(:,mm));
        end
        for mm = 1:length(aligned_baseline_EMG_noon)
            per_trial_avg_baseline_EMG_noon(mm,ii) = mean(all_trials_baseline_EMG_noon{ii,1}(:,mm));
        end
    end
    
    %% Calculating average EMG (Average Across trials)
    cross_trial_avg_baseline_EMG_morn = struct([]);
    cross_trial_avg_baseline_EMG_noon = struct([]);
    for ii = 1:length(M)
        cross_trial_avg_baseline_EMG_morn{ii,1} = zeros(length(aligned_baseline_EMG_morn{1,1}),1);
        cross_trial_avg_baseline_EMG_noon{ii,1} = zeros(length(aligned_baseline_EMG_noon{1,1}),1);
        for mm = 1:length(aligned_baseline_EMG_morn{1,1})
            cross_trial_avg_baseline_EMG_morn{ii,1}(mm) = mean(all_trials_baseline_EMG_morn{ii,1}(mm,:));
        end
        for mm = 1:length(aligned_baseline_EMG_noon{1,1})
            cross_trial_avg_baseline_EMG_noon{ii,1}(mm) = mean(all_trials_baseline_EMG_noon{ii,1}(mm,:));
        end
    end

    %% Calculate the output variables
    for ii = 1:length(M)
        [~, Baseline_EMG_p_values(ii, jj)] = ...
            ttest2(per_trial_avg_baseline_EMG_morn(:,ii), per_trial_avg_baseline_EMG_noon(:,ii));
        Morn_Baseline_EMG(ii, jj) = mean(per_trial_avg_baseline_EMG_morn(:,ii));
        STD_Morn_Baseline_EMG(ii, jj) = std(per_trial_avg_baseline_EMG_morn(:,ii));
        Err_Morn_Baseline_EMG(ii, jj) =  STD_Morn_Baseline_EMG(ii, jj) / sqrt(length(per_trial_avg_baseline_EMG_morn(:,ii)));
        Noon_Baseline_EMG(ii, jj) = mean(per_trial_avg_baseline_EMG_noon(:,ii));
        STD_Noon_Baseline_EMG(ii, jj) = std(per_trial_avg_baseline_EMG_noon(:,ii));
        Err_Noon_Baseline_EMG(ii, jj) =  STD_Noon_Baseline_EMG(ii, jj) / sqrt(length(per_trial_avg_baseline_EMG_noon(:,ii)));
        % EMG percent change
        if all(EMG_Norm_Factor + 1)
            Baseline_EMG_perc_changes(ii, jj) = (Noon_Baseline_EMG(ii, jj) - Morn_Baseline_EMG(ii, jj)) / 100;
        else
            Baseline_EMG_perc_changes(ii, jj) = ...
                (Noon_Baseline_EMG(ii, jj) - Morn_Baseline_EMG(ii, jj)) / abs(Morn_Baseline_EMG(ii, jj));
        end
    end

    %% Combining the morning and afternoon into one matrix

    avg_baseline_EMG = struct([]);
    Baseline_EMG_mean = zeros(length(M),2);
    Baseline_EMG_STD = zeros(length(M),2);
    for ii = 1:length(M)
        avg_baseline_EMG{ii,1}(:,1) = cross_trial_avg_baseline_EMG_morn{ii};
        avg_baseline_EMG{ii,1}(:,2) = cross_trial_avg_baseline_EMG_noon{ii};
        Baseline_EMG_mean(ii,1) = Morn_Baseline_EMG(ii);
        Baseline_EMG_mean(ii,2) = Noon_Baseline_EMG(ii);
        Baseline_EMG_STD(ii,1) = STD_Morn_Baseline_EMG(ii);
        Baseline_EMG_STD(ii,2) = STD_Noon_Baseline_EMG(ii);
    end

    % Finding the min and max for the y-axis limit
    Baseline_y_min = Baseline_EMG_mean - Baseline_EMG_STD;
    Baseline_y_max = Baseline_EMG_mean + Baseline_EMG_STD;

    if isequal(Plot_Figs, 1)

        %% Plot the EMG box and whisker plot

        for ii = 1:length(M)

            figure
            subplot(211) % Top Plot
            hold on
            boxplot(avg_baseline_EMG{ii})

            % Setting the y-axis limits
            y_max = max(Baseline_y_max(ii,:));
            y_min = min(Baseline_y_min(ii,:));
            ylim([y_min - abs(y_min/8), y_max + abs(y_max/8)])
            % Setting the x-axis limits
            xlim([0.5, 2.5]);

            % Titling the plot
            EMG_title = strrep(string(xds_morn.EMG_names(M(ii))),'EMG_','');
            if ~isequal(per_dir_EMG, 0)
                title(sprintf('Baseline EMG, %i°, TgtCenter at %0.1f: %s', ...
                    target_dir_noon(jj), unique_targets_noon(kk), EMG_title), 'FontSize', title_font_size)
            else
                title(sprintf('Baseline EMG, All Trials: %s', ...
                    EMG_title), 'FontSize', title_font_size)
            end      

            % Get the top subplot title for saving
            if ~isequal(Save_Figs, 0)
                fig_info = get(gca,'title');
                save_title{ii} = get(fig_info, 'string');
            end

            % Annotation of the p_value
            if round(Baseline_EMG_p_values(ii, jj), 3) > 0
                legend_dims = [0 0.45 0.44 0.44];
                p_value_string = strcat('p =', {' '}, mat2str(round(Baseline_EMG_p_values(ii, jj), 3)));
                legend_string = {char(p_value_string)};
                ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                    'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                    'EdgeColor','none', 'horizontalalignment', 'center');
                ann_legend.FontSize = legend_font_size;
                ann_legend.FontName = font_name;
            end

            if isequal(round(Baseline_EMG_p_values(ii, jj), 3), 0)
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
            if ~isequal(round(Baseline_EMG_perc_changes(ii, jj), 3), 0)
                legend_dims = [0.55 0.45 0.44 0.44];
                perc_change_string = strcat('Δ% =', {' '}, mat2str(round(Baseline_EMG_perc_changes(ii, jj), 3)));
                legend_string = {char(perc_change_string)};
                ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                    'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                    'EdgeColor','none', 'horizontalalignment', 'center');
                ann_legend.FontSize = legend_font_size;
                ann_legend.FontName = font_name;
            end
            if isequal(round(Baseline_EMG_perc_changes(ii, jj), 3), 0)
                legend_dims = [0.55 0.45 0.44 0.44];
                perc_change_string = strcat('Δ% ≈', {' '}, '0');
                legend_string = {char(perc_change_string)};
                ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                    'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                    'EdgeColor','none', 'horizontalalignment', 'center');
                ann_legend.FontSize = legend_font_size;
                ann_legend.FontName = font_name;
            end

            per_trial_ymax_morn = max(all_trials_baseline_EMG_morn{ii,1}, [], 'All');
            per_trial_ymax_noon = max(all_trials_baseline_EMG_noon{ii,1}, [], 'All');
            per_trial_ymin_morn = min(all_trials_baseline_EMG_morn{ii,1}, [], 'All');
            per_trial_ymin_noon = min(all_trials_baseline_EMG_noon{ii,1}, [], 'All');
            y_max = max(per_trial_ymax_morn, per_trial_ymax_noon);
            y_min = min(per_trial_ymin_morn, per_trial_ymin_noon);

            subplot(223) % Bottom Left Plot
            hold on
            % Set the title
            title('Morning', 'FontSize', title_font_size);
            for pp = 1:width(all_trials_baseline_EMG_morn{ii})
                plot(absolute_baseline_EMG_timing_morn, all_trials_baseline_EMG_morn{ii}(:,pp))
            end

            % Annotation of number of trials
            legend_dims = [0.175 0.02 0.44 0.44];
            trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_baseline_EMG_morn{ii})));
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
            for pp = 1:width(all_trials_baseline_EMG_noon{ii})
                plot(absolute_baseline_EMG_timing_noon, all_trials_baseline_EMG_noon{ii}(:,pp))
            end

            % Annotation of number of trials
            legend_dims = [0.625 0.02 0.44 0.44];
            trial_num_string = strcat('n =', {' '}, mat2str(width(all_trials_baseline_EMG_noon{ii})));
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

        end % End of the EMG loop

    end % End of the Plot if-statement

    % End the function if you only want one baseline EMG
    if ~isequal(per_dir_EMG, 1) && isequal(jj, 1)
        return
    end

end % End of target loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        save_title{ii} = strrep(save_title{ii}, ':', '');
        save_title{ii} = strrep(save_title{ii}, 'vs.', 'vs');
        save_title{ii} = strrep(save_title{ii}, 'mg.', 'mg');
        save_title{ii} = strrep(save_title{ii}, 'kg.', 'kg');
        save_title{ii} = strrep(save_title{ii}, '.', '_');
        save_title{ii} = strrep(save_title{ii}, '/', '_');
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
