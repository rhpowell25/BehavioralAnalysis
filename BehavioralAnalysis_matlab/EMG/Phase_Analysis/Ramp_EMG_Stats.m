function [EMG_Names, Ramp_EMG_p_values, Ramp_EMG_perc_changes, ...
    Morn_Ramp_EMG, Err_Morn_Ramp_EMG, Noon_Ramp_EMG, Err_Noon_Ramp_EMG] = ...
    Ramp_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
    EMG_Zero_Factor, EMG_Norm_Factor, Plot_Figs, Save_File)

%% Display the function being used
disp('Ramp EMG Statistics:');

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

% Save Counter
ss = 1;
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
        Morn_Ramp_EMG = zeros(length(M), num_dirs);
        STD_Morn_Ramp_EMG = zeros(length(M), num_dirs);
        Err_Morn_Ramp_EMG = zeros(length(M), num_dirs);
        Noon_Ramp_EMG = zeros(length(M), num_dirs);
        STD_Noon_Ramp_EMG = zeros(length(M), num_dirs);
        Err_Noon_Ramp_EMG = zeros(length(M), num_dirs);
        Ramp_EMG_p_values = zeros(length(M), num_dirs);
        Ramp_EMG_perc_changes = zeros(length(M), num_dirs);
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
    
    aligned_end_EMG_morn = struct([]); % EMG during each successful trial
    aligned_end_EMG_timing_morn = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_morn)
        aligned_end_EMG_morn{ii, 1} = xds_morn.EMG(rewarded_start_idx_morn(ii) : rewarded_end_idx_morn(ii), M);
        aligned_end_EMG_timing_morn{ii, 1} = xds_morn.time_frame(rewarded_start_idx_morn(ii) : rewarded_end_idx_morn(ii));
    end

    aligned_end_EMG_noon = struct([]); % EMG during each successful trial
    aligned_end_EMG_timing_noon = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_end_time_noon)
        aligned_end_EMG_noon{ii, 1} = xds_noon.EMG(rewarded_start_idx_noon(ii) : rewarded_end_idx_noon(ii), M);
        aligned_end_EMG_timing_noon{ii, 1} = xds_noon.time_frame(rewarded_start_idx_noon(ii) : rewarded_end_idx_noon(ii));
    end

    %% Zeroing the EMG

    for ii = 1:length(aligned_end_EMG_morn)
        for mm = 1:length(M)
            aligned_end_EMG_morn{ii,1}(:,mm) = aligned_end_EMG_morn{ii,1}(:,mm) - EMG_Zero_Factor(mm);
        end
    end
    for ii = 1:length(aligned_end_EMG_noon)
        for mm = 1:length(M)
            aligned_end_EMG_noon{ii,1}(:,mm) = aligned_end_EMG_noon{ii,1}(:,mm) - EMG_Zero_Factor(mm);
        end
    end

    %% Normalizing the EMG

    for ii = 1:length(aligned_end_EMG_morn)
        for mm = 1:length(M)
            aligned_end_EMG_morn{ii,1}(:,mm) = (aligned_end_EMG_morn{ii,1}(:,mm) / EMG_Norm_Factor(mm))*100;
        end
    end
    for ii = 1:length(aligned_end_EMG_noon)
        for mm = 1:length(M)
            aligned_end_EMG_noon{ii,1}(:,mm) = (aligned_end_EMG_noon{ii,1}(:,mm) / EMG_Norm_Factor(mm))*100;
        end
    end

    %% Calculating average EMG (Average per trial)
    per_trial_avg_end_EMG_morn = zeros(length(aligned_end_EMG_morn), length(M));
    per_trial_avg_end_EMG_noon = zeros(length(aligned_end_EMG_noon), length(M));
    for ii = 1:length(aligned_end_EMG_morn)
        for mm = 1:length(M)
            per_trial_avg_end_EMG_morn(ii,mm) = mean(aligned_end_EMG_morn{ii}(:,mm));
        end
    end
    for ii = 1:length(aligned_end_EMG_noon)
        for mm = 1:length(M)
            per_trial_avg_end_EMG_noon(ii,mm) = mean(aligned_end_EMG_noon{ii}(:,mm));
        end
    end

    %% Calculate the output variables
    for ii = 1:length(M)
        [~, Ramp_EMG_p_values(ii, jj)] = ...
            ttest2(per_trial_avg_end_EMG_morn(:,ii), per_trial_avg_end_EMG_noon(:,ii));
        % EMG means
        Morn_Ramp_EMG(ii, jj) = mean(per_trial_avg_end_EMG_morn(:,ii));
        Noon_Ramp_EMG(ii, jj) = mean(per_trial_avg_end_EMG_noon(:,ii));
        % EMG standard devs
        STD_Morn_Ramp_EMG(ii, jj) = std(per_trial_avg_end_EMG_morn(:,ii));
        Err_Morn_Ramp_EMG(ii, jj) =  STD_Morn_Ramp_EMG(ii, jj) / sqrt(length(per_trial_avg_end_EMG_morn(:,ii)));
        STD_Noon_Ramp_EMG(ii, jj) = std(per_trial_avg_end_EMG_noon(:,ii));
        Err_Noon_Ramp_EMG(ii, jj) =  STD_Noon_Ramp_EMG(ii, jj) / sqrt(length(per_trial_avg_end_EMG_noon(:,ii)));
        % EMG percent change
        Ramp_EMG_perc_changes(ii, jj) = ...
            (Noon_Ramp_EMG(ii, jj) - Morn_Ramp_EMG(ii, jj)) / abs(Morn_Ramp_EMG(ii, jj));
    end

    %% Plot the ramp phase EMG

    if isequal(Plot_Figs, 1)

        %% Plot the EMG box and whisker plot
        for ii = 1:length(M)

            % Combining the morning and afternoon into one matrix
            avg_ramp_EMG = cat(1, per_trial_avg_end_EMG_morn(:,ii), per_trial_avg_end_EMG_noon(:,ii));
            ramp_EMG_labels = zeros(length(avg_ramp_EMG),1);
            ramp_EMG_labels(length(per_trial_avg_end_EMG_morn) + 1:end) = 1;
            
            End_EMG_mean(1,1) = Morn_Ramp_EMG(ii, jj);
            End_EMG_mean(1,2) = Noon_Ramp_EMG(ii, jj);
            End_EMG_STD(1,1) = STD_Morn_Ramp_EMG(ii, jj);
            End_EMG_STD(1,2) = STD_Noon_Ramp_EMG(ii, jj);
    
            figure
            subplot(211) % Top Plot
            hold on
            boxplot(avg_ramp_EMG, ramp_EMG_labels)
    
            % Setting the x-axis limits
            xlim([0.5, 2.5]);
    
            % Titling the plot
            Fig_Title = strrep(string(xds_morn.EMG_names(M(ii))),'EMG_','');
            title(sprintf('Ramp phase EMG, %i°, TgtCenter at %0.1f: %s', ...
                target_dirs_noon(jj), target_centers_noon(jj), Fig_Title), 'FontSize', title_font_size)
    
            % Annotation of the p_value
            if round(Ramp_EMG_p_values(ii, jj), 3) > 0
                legend_dims = [0 0.45 0.44 0.44];
                p_value_string = strcat('p =', {' '}, mat2str(round(Ramp_EMG_p_values(ii, jj), 3)));
                legend_string = {char(p_value_string)};
                ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                    'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                    'EdgeColor','none', 'horizontalalignment', 'center');
                ann_legend.FontSize = legend_font_size;
                ann_legend.FontName = font_name;
            end
            if isequal(round(Ramp_EMG_p_values(ii, jj), 3), 0)
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
            if ~isequal(round(Ramp_EMG_perc_changes(ii, jj), 3), 0)
                legend_dims = [0.55 0.45 0.44 0.44];
                perc_change_string = strcat('Δ% =', {' '}, mat2str(round(Ramp_EMG_perc_changes(ii, jj), 3)));
                legend_string = {char(perc_change_string)};
                ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                    'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                    'EdgeColor','none', 'horizontalalignment', 'center');
                ann_legend.FontSize = legend_font_size;
                ann_legend.FontName = font_name;
            end
            if isequal(round(Ramp_EMG_perc_changes(ii, jj), 3), 0)
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
    
            subplot(223) % Bottom Left Plot
            hold on
            % Set the title
            title('Morning', 'FontSize', title_font_size);
            for pp = 1:length(aligned_end_EMG_morn)
                plot(aligned_end_EMG_morn{pp,1}(:,ii))
            end
    
            % Annotation of number of trials
            legend_dims = [0.175 0.02 0.44 0.44];
            trial_num_string = strcat('n =', {' '}, mat2str(length(aligned_end_EMG_morn)));
            legend_string = {char(trial_num_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
    
            % Set ticks to outside
            figure_axes = gca;
            set(figure_axes,'TickDir','out');
            % Remove the top and right tick marks
            set(figure_axes,'box','off')
            % Set The Font
            set(figure_axes,'FontName', font_name);
    
            subplot(224) % Bottom Right Plot
            hold on
            % Set the title
            title('Afternoon', 'FontSize', title_font_size);
            for pp = 1:length(aligned_end_EMG_noon)
                plot(aligned_end_EMG_noon{pp,1}(:,ii))
            end
    
            % Annotation of number of trials
            legend_dims = [0.625 0.02 0.44 0.44];
            trial_num_string = strcat('n =', {' '}, mat2str(length(aligned_end_EMG_noon)));
            legend_string = {char(trial_num_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
    
            % Set ticks to outside
            figure_axes = gca;
            set(figure_axes,'TickDir','out');
            % Remove the top and right tick marks
            set(figure_axes,'box','off')
            % Set The Font
            set(figure_axes,'FontName', font_name);

            ss = ss + 1;

        end % End of the EMG loop

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End of the Plot if-statement

end % End of target loop