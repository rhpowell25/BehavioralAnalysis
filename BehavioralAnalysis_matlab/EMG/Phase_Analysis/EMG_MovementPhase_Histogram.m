function [EMG_Names, TgtHold_EMG_p_values, TgtHold_EMG_perc_changes] = ...
    EMG_MovementPhase_Histogram(xds_morn, xds_noon, muscle_groups, EMG_Zero_Factor, EMG_Norm_Factor, Save_Figs)

%% Display the function being used
disp('Movement Phase EMG Histogram:');

%% Find the EMG index

[M] = EMG_Index(xds_morn, muscle_groups);

EMG_Names = strings;
for ii = 1:length(M)
    EMG_Names(ii,1) = xds_morn.EMG_names(M(ii));
end

%% Basic Settings, some variable extractions, & definitions

% Do you want to plot the morning / afternoon legend? (1 = Yes, 0 = No)
plot_legend = 1;

% Font specifications
label_font_size = 30;
legend_size = 30;
mean_line_width = 5;
if isequal(plot_legend, 1)
    p_value_dims = [0.51 0.3 0.44 0.44];
else
    p_value_dims = [0.51 0.45 0.44 0.44];
end
if isequal(plot_legend, 1)
    perc_change_dims = [0.49 0.225 0.44 0.44];
else
    perc_change_dims = [0.49 0.375 0.44 0.44];
end
title_font_size = 15;
font_name = 'Arial';

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

    if xds_morn.has_force == 1
        [movement_start_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'force_onset');
        [movement_start_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'force_onset');
    else
        [movement_start_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'cursor_onset');
        [movement_start_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'cursor_onset');
    end
    %% Define the output variables
    if jj == 1
        Morn_TgtHold_EMG = zeros(length(M), num_dirs);
        STD_Morn_TgtHold_EMG = zeros(length(M), num_dirs);
        Err_Morn_TgtHold_EMG = zeros(length(M), num_dirs);
        Noon_TgtHold_EMG = zeros(length(M), num_dirs);
        STD_Noon_TgtHold_EMG = zeros(length(M), num_dirs);
        Err_Noon_TgtHold_EMG = zeros(length(M), num_dirs);
        TgtHold_EMG_p_values = zeros(length(M), num_dirs);
        TgtHold_EMG_perc_changes = zeros(length(M), num_dirs);
    end
        
    %% Find the rewarded indices
    rewarded_tgt_hold_idx_morn = zeros(height(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        rewarded_tgt_hold_idx_morn(ii) = find(xds_morn.time_frame == ...
            round(rewarded_end_time_morn(ii) - tgt_hold_time, abs(floor(log10(xds_morn.bin_width)))));
    end

    rewarded_movement_start_idx_morn = zeros(height(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        rewarded_movement_start_idx_morn(ii) = find(xds_morn.time_frame == movement_start_time_morn(ii));
    end

    rewarded_tgt_hold_idx_noon = zeros(height(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        rewarded_tgt_hold_idx_noon(ii) = find(xds_noon.time_frame == ...
            round(rewarded_end_time_noon(ii) - tgt_hold_time, abs(floor(log10(xds_noon.bin_width)))));
    end

    rewarded_movement_start_idx_noon = zeros(height(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        rewarded_movement_start_idx_noon(ii) = find(xds_noon.time_frame == movement_start_time_noon(ii));
    end
    
    %% EMG and time aligned to specified event
    aligned_end_EMG_morn = struct([]); % EMG during each successful trial
    aligned_end_EMG_timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_morn)
        aligned_end_EMG_morn{ii, 1} = xds_morn.EMG((rewarded_movement_start_idx_morn(ii) : ...
            rewarded_tgt_hold_idx_morn(ii)), :);
        aligned_end_EMG_timing_morn{ii, 1} = xds_morn.time_frame((rewarded_movement_start_idx_morn(ii) : ...
            rewarded_tgt_hold_idx_morn(ii)));
    end

    aligned_end_EMG_noon = struct([]); % EMG during each successful trial
    aligned_end_EMG_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_noon)
        aligned_end_EMG_noon{ii, 1} = xds_noon.EMG((rewarded_movement_start_idx_noon(ii) : ...
            rewarded_tgt_hold_idx_noon(ii)), :);
        aligned_end_EMG_timing_noon{ii, 1} = xds_noon.time_frame((rewarded_movement_start_idx_noon(ii) : ...
            rewarded_tgt_hold_idx_noon(ii)));
    end

    %% Putting all succesful trials in one array
    all_trials_end_EMG_morn = struct([]);
    for ii = 1:length(M)
        all_trials_end_EMG_morn{ii,1} = zeros(length(aligned_end_EMG_morn{1,1}),length(rewarded_end_time_morn));
        for mm = 1:length(rewarded_end_time_morn)
            all_trials_end_EMG_morn{ii,1}(:,mm) = aligned_end_EMG_morn{mm, 1}(:, M(ii));
        end
    end

    all_trials_end_EMG_noon = struct([]);
    for ii = 1:length(M)
        all_trials_end_EMG_noon{ii,1} = zeros(length(aligned_end_EMG_noon{1,1}),length(rewarded_end_time_noon));
        for mm = 1:length(rewarded_end_time_noon)
            all_trials_end_EMG_noon{ii,1}(:,mm) = aligned_end_EMG_noon{mm, 1}(:, M(ii));
        end
    end

    %% Zeroing the EMG
    for ii = 1:length(M)
        all_trials_end_EMG_morn{ii,1} = all_trials_end_EMG_morn{ii,1} - EMG_Zero_Factor(ii);
        all_trials_end_EMG_noon{ii,1} = all_trials_end_EMG_noon{ii,1} - EMG_Zero_Factor(ii);
    end

    %% Normalizing the average EMG's
    for ii = 1:length(M)
        all_trials_end_EMG_morn{ii,1} = (all_trials_end_EMG_morn{ii,1} / EMG_Norm_Factor(ii))*100;
        all_trials_end_EMG_noon{ii,1} = (all_trials_end_EMG_noon{ii,1} / EMG_Norm_Factor(ii))*100;
    end

    %% Calculating average EMG (Average per trial)
    per_trial_avg_end_EMG_morn = zeros(length(aligned_end_EMG_morn),length(M));
    per_trial_avg_end_EMG_noon = zeros(length(aligned_end_EMG_noon),length(M));
    for ii = 1:length(M)
        for mm = 1:length(aligned_end_EMG_morn)
            per_trial_avg_end_EMG_morn(mm,ii) = mean(all_trials_end_EMG_morn{ii,1}(:,mm));
        end
        for mm = 1:length(aligned_end_EMG_noon)
            per_trial_avg_end_EMG_noon(mm,ii) = mean(all_trials_end_EMG_noon{ii,1}(:,mm));
        end
    end
    
    %% Calculating average EMG (Average Across trials)
    cross_trial_avg_end_EMG_morn = struct([]);
    cross_trial_avg_end_EMG_noon = struct([]);
    for ii = 1:length(M)
        cross_trial_avg_end_EMG_morn{ii,1} = zeros(length(aligned_end_EMG_morn{1,1}),1);
        cross_trial_avg_end_EMG_noon{ii,1} = zeros(length(aligned_end_EMG_noon{1,1}),1);
        for mm = 1:length(aligned_end_EMG_morn{1,1})
            cross_trial_avg_end_EMG_morn{ii,1}(mm) = mean(all_trials_end_EMG_morn{ii,1}(mm,:));
        end
        for mm = 1:length(aligned_end_EMG_noon{1,1})
            cross_trial_avg_end_EMG_noon{ii,1}(mm) = mean(all_trials_end_EMG_noon{ii,1}(mm,:));
        end
    end

    %% Calculate the output variables
    for ii = 1:length(M)
        [~, TgtHold_EMG_p_values(ii, jj)] = ...
            ttest2(per_trial_avg_end_EMG_morn(:,ii), per_trial_avg_end_EMG_noon(:,ii));
        % EMG means
        Morn_TgtHold_EMG(ii, jj) = mean(per_trial_avg_end_EMG_morn(:,ii));
        Noon_TgtHold_EMG(ii, jj) = mean(per_trial_avg_end_EMG_noon(:,ii));
        % EMG standard devs
        STD_Morn_TgtHold_EMG(ii, jj) = std(per_trial_avg_end_EMG_morn(:,ii));
        Err_Morn_TgtHold_EMG(ii, jj) =  STD_Morn_TgtHold_EMG(ii, jj) / sqrt(length(per_trial_avg_end_EMG_morn(:,ii)));
        STD_Noon_TgtHold_EMG(ii, jj) = std(per_trial_avg_end_EMG_noon(:,ii));
        Err_Noon_TgtHold_EMG(ii, jj) =  STD_Noon_TgtHold_EMG(ii, jj) / sqrt(length(per_trial_avg_end_EMG_noon(:,ii)));
        % EMG percent change
        TgtHold_EMG_perc_changes(ii, jj) = ...
            (Noon_TgtHold_EMG(ii, jj) - Morn_TgtHold_EMG(ii, jj)) / abs(Morn_TgtHold_EMG(ii, jj));
    end

    %% Plot the target hold EMG

    % Combining the morning and afternoon into one matrix
    avg_end_EMG = struct([]);
    End_EMG_mean = zeros(length(M),2);
    End_EMG_STD = zeros(length(M),2);
    for ii = 1:length(M)
        avg_end_EMG{ii,1}(:,1) = cross_trial_avg_end_EMG_morn{ii};
        avg_end_EMG{ii,1}(:,2) = cross_trial_avg_end_EMG_noon{ii};
        End_EMG_mean(ii,1) = Morn_TgtHold_EMG(ii, jj);
        End_EMG_mean(ii,2) = Noon_TgtHold_EMG(ii, jj);
        End_EMG_STD(ii,1) = STD_Morn_TgtHold_EMG(ii, jj);
        End_EMG_STD(ii,2) = STD_Noon_TgtHold_EMG(ii, jj);
    end

    %% Plot the EMG histogram

    for ii = 1:length(M)

        EMG_Hist_fig = figure;
        EMG_Hist_fig.Position = [200 50 800 600];
        hold on

        histogram(all_trials_end_EMG_morn{ii}, 15, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
        histogram(all_trials_end_EMG_noon{ii}, 15, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])

        fprintf('The EMG percent change is %0.2f percent \n', TgtHold_EMG_perc_changes(ii, jj))

        % Set the axis
        x_limits = xlim;
        y_limits = ylim;
        xlim([x_limits(1) - 10, x_limits(2) + 10])
        ylim([y_limits(1), y_limits(2)])
        
        % Plot dummy points for the legend
        if isequal(plot_legend, 1)
            dummy_morn = plot(-10,-10, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], 'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
            dummy_noon = plot(-20,-10, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
                    'MarkerEdgeColor',[.5 0 .5], 'MarkerFaceColor',[.5 0 .5]);
        end

        % Plot the means
        line([mean(avg_end_EMG{ii}(:,1)) mean(avg_end_EMG{ii}(:,1))], [y_limits(1) y_limits(2) + 0.25], ... 
            'LineStyle','--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', mean_line_width)
        line([mean(avg_end_EMG{ii}(:,2)) mean(avg_end_EMG{ii}(:,2))], [y_limits(1) y_limits(2) + 0.25], ... 
            'LineStyle','--', 'Color', [.5 0 .5], 'LineWidth', mean_line_width)

        % Set the title
        EMG_title = strrep(string(xds_morn.EMG_names(M(ii))),'EMG_','');
        title(sprintf('Movement Phase EMG, %i°, TgtCenter at %0.1f: %s', ...
            target_dirs_noon(jj), target_centers_noon(jj), EMG_title), 'FontSize', title_font_size)

        % Get the top subplot title for saving
        if ~isequal(Save_Figs, 0)
            fig_info = get(gca,'title');
            save_title{ii} = get(fig_info, 'string');
        end

        % Annotation of the p_value
        if round(TgtHold_EMG_p_values(ii, jj), 3) > 0
            p_value_string = strcat('p =', {' '}, mat2str(round(TgtHold_EMG_p_values(ii, jj), 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', p_value_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_size;
            ann_legend.FontName = font_name;
        end

        if isequal(round(TgtHold_EMG_p_values(ii, jj), 3), 0)
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', p_value_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_size;
            ann_legend.FontName = font_name;
        end

        % Annotation of the percent change
        if ~isequal(round(TgtHold_EMG_perc_changes(ii, jj), 2), 0)
            perc_change_string = strcat('Δ% =', {' '}, mat2str(round(TgtHold_EMG_perc_changes(ii, jj), 3)));
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', perc_change_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_size;
            ann_legend.FontName = font_name;
        end
        if isequal(round(TgtHold_EMG_perc_changes(ii, jj), 2), 0)
            perc_change_string = strcat('Δ% ≈', {' '}, '0');
            legend_string = {char(perc_change_string)};
            ann_legend = annotation('textbox', perc_change_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_size;
            ann_legend.FontName = font_name;
        end
        
        % Axis Editing
        figure_axes = gca;
        % Set ticks to outside
        set(figure_axes,'TickDir','out');
        % Remove the top and right tick marks
        set(figure_axes,'box','off')
        % Set the tick label font size
        figure_axes.FontSize = label_font_size - 15;

        % Label the axis
        ylabel('# of EMG Samples', 'FontSize', label_font_size)
        xlabel('EMG', 'FontSize', label_font_size)

        if isequal(plot_legend, 1)
            legend([dummy_morn, dummy_noon], ... 
                {'Morning', 'Afternoon'}, ... 
                'FontSize', legend_size, 'Location', 'NorthEast')
            legend boxoff
        end

        % Only label every other tick
        %x_labels = string(figure_axes.XAxis.TickLabels);
        %y_labels = string(figure_axes.YAxis.TickLabels);
        %x_labels(2:2:end) = NaN;
        %y_labels(2:2:end) = NaN;
        %figure_axes.XAxis.TickLabels = x_labels;
        %figure_axes.YAxis.TickLabels = y_labels;

    end % End of the EMG loop

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
            saveas(gcf, fullfile(save_dir, char(save_title{ii})), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title{ii})), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title{ii})), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title{ii})), 'fig')
        end
        close gcf
    end
end
