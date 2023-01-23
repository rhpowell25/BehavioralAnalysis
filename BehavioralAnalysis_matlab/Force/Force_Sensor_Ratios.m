function Force_Sensor_Ratios(xds_morn, xds_noon, norm_force, Save_Figs)

%% Display the function being used

disp('Force Sensor Ratios:');

%% Basic Settings, some variable extractions, & definitions

% Font specifications
label_font_size = 17;
title_font_size = 14;
mean_line_width = 5;
p_value_dims = [0.575 0.375 0.44 0.44];
legend_size = 15;
font_name = 'Arial';

% Save Counter
if ~isequal(Save_Figs, 0)
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
    
    First_Force_morn = struct([]); % Force during each successful trial
    Second_Force_morn = struct([]);
    for ii = 1:length(rewarded_end_time_morn)
        First_Force_morn{ii, 1} = xds_morn.force((rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width) : ...
            rewarded_end_idx_morn(ii)), 1);
        Second_Force_morn{ii, 1} = xds_morn.force((rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width) : ...
            rewarded_end_idx_morn(ii)), 2);
    end
    
    First_Force_noon = struct([]); % Force during each successful trial
    Second_Force_noon = struct([]);
    for ii = 1:length(rewarded_end_time_noon)
        First_Force_noon{ii, 1} = xds_noon.force((rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width) : ...
            rewarded_end_idx_noon(ii)), 1);
        Second_Force_noon{ii, 1} = xds_noon.force((rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width) : ...
            rewarded_end_idx_noon(ii)), 2);
    end
    
    %% Putting all succesful trials in one array
    
    all_trials_First_Force_morn = zeros(length(First_Force_morn{1,1}), length(rewarded_end_time_morn));
    all_trials_Second_Force_morn = zeros(length(Second_Force_morn{1,1}), length(rewarded_end_time_morn));
    for ii = 1:length(rewarded_end_time_morn)
        all_trials_First_Force_morn(:,ii) = First_Force_morn{ii, 1};
        all_trials_Second_Force_morn(:,ii) = Second_Force_morn{ii, 1};
    end
    
    all_trials_First_Force_noon = zeros(length(First_Force_noon{1,1}), length(rewarded_end_time_noon));
    all_trials_Second_Force_noon = zeros(length(Second_Force_noon{1,1}), length(rewarded_end_time_noon));
    for ii = 1:length(rewarded_end_time_noon)
        all_trials_First_Force_noon(:,ii) = First_Force_noon{ii, 1};
        all_trials_Second_Force_noon(:,ii) = Second_Force_noon{ii, 1};
    end
    
    %% Normalizing the average force
    if isequal(norm_force, 1)
        Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force);
        all_trials_First_Force_morn = all_trials_First_Force_morn / Force_Norm_Factor*100;
        all_trials_Second_Force_morn = all_trials_Second_Force_morn / Force_Norm_Factor*100;
        all_trials_First_Force_noon = all_trials_First_Force_noon / Force_Norm_Factor*100;
        all_trials_Second_Force_noon = all_trials_Second_Force_noon / Force_Norm_Factor*100;
    elseif strcmp(norm_force, 'Convert')
        all_trials_First_Force_morn = all_trials_First_Force_morn / 1000*5; % Millivolt conversion * gain
        all_trials_Second_Force_morn = all_trials_Second_Force_morn / 1000*5; % Millivolt conversion * gain
        all_trials_First_Force_noon = all_trials_First_Force_noon / 1000*5; % Millivolt conversion * gain
        all_trials_Second_Force_noon = all_trials_Second_Force_noon / 1000*5; % Millivolt conversion * gain
    end
    
    %% Calculating average Force (Average per trial)
    per_trial_avg_Force_1_morn = zeros(length(First_Force_morn), 1);
    per_trial_avg_Force_2_morn = zeros(length(Second_Force_morn), 1);
    per_trial_avg_Force_1_noon = zeros(length(First_Force_noon), 1);
    per_trial_avg_Force_2_noon = zeros(length(Second_Force_noon), 1);
    for ii = 1:length(First_Force_morn)
        per_trial_avg_Force_1_morn(ii,1) = mean(all_trials_First_Force_morn(:,ii));
        per_trial_avg_Force_2_morn(ii,1) = mean(all_trials_Second_Force_morn(:,ii));
    end
    for ii = 1:length(First_Force_noon)
        per_trial_avg_Force_1_noon(ii,1) = mean(all_trials_First_Force_noon(:,ii));
        per_trial_avg_Force_2_noon(ii,1) = mean(all_trials_Second_Force_noon(:,ii));
    end
    
    %% Calculate the ratios between the two force sensors
    % Per trial ratios
    per_trial_Force_ratio_morn = zeros(length(rewarded_end_time_morn),1);
    for ii = 1:length(rewarded_end_time_morn)
        per_trial_Force_ratio_morn(ii) = per_trial_avg_Force_1_morn(ii) / per_trial_avg_Force_2_morn(ii);
    end
    per_trial_Force_ratio_noon = zeros(length(rewarded_end_time_noon),1);
    for ii = 1:length(rewarded_end_time_noon)
        per_trial_Force_ratio_noon(ii) = per_trial_avg_Force_1_noon(ii) / per_trial_avg_Force_2_noon(ii);
    end
    
    %% Calculate the average force ratios
    Force_ratio_morn = mean(per_trial_Force_ratio_morn);
    Force_ratio_noon = mean(per_trial_Force_ratio_noon);
    
    %% Plot the force sensor ratio histograms
    
    force_fig = figure;
    force_fig.Position = [200 50 800 600];
    hold on
    
    % Label the axis
    xlabel('Force Sensor Ratios', 'FontSize', label_font_size);
    ylabel('Trial Numbers', 'FontSize', label_font_size);
    
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
    
    % Set the title
    title_string = strcat(Date, {' '}, Task, ',', {' '}, Drug, ': Force Sensor Ratios');
    title(strcat(title_string, sprintf(': %i°, TgtCenter at %0.1f', ... 
        target_dirs_morn(jj), target_centers_morn(jj))), 'FontSize', title_font_size)
    
    histogram(per_trial_Force_ratio_morn, 15, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(per_trial_Force_ratio_noon, 15, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])
    
    % Set the axis
    y_limits = ylim;
    ylim([y_limits(1), y_limits(2) + 0.25])
    
    % Plot the means
    line([Force_ratio_morn Force_ratio_morn], [y_limits(1) y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', mean_line_width)
    line([Force_ratio_noon Force_ratio_noon], [y_limits(1) y_limits(2) + 0.25], ... 
        'LineStyle','--', 'Color', [.5 0 .5], 'LineWidth', mean_line_width)
    
    % Plot dummy points for the legend
    dummy_morn = plot(-1,-1, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
            'MarkerEdgeColor',[0.9290, 0.6940, 0.1250], 'MarkerFaceColor',[0.9290, 0.6940, 0.1250]);
    dummy_noon = plot(-2,-1, 's', 'MarkerSize',20, 'LineWidth', 1.5, ...
            'MarkerEdgeColor',[.5 0 .5], 'MarkerFaceColor',[.5 0 .5]);
    
    % Plot the legend
    legend([dummy_morn, dummy_noon], ... 
        {'Morning', 'Afternoon'}, ... 
        'FontSize', legend_size, 'Location', 'NorthEast')
    legend boxoff
    
    % Annotation of the p_value
    [~, force_ratio_p_val] = ttest2(per_trial_Force_ratio_morn, per_trial_Force_ratio_noon);
    if round(force_ratio_p_val, 3) > 0
        p_value_string = strcat('p =', {' '}, mat2str(round(force_ratio_p_val, 3)));
        p_value_string = {char(p_value_string)};
        ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_p_value.FontSize = legend_size;
        ann_p_value.FontName = font_name;
    end
    
    if isequal(round(force_ratio_p_val, 3), 0)
        p_value_string = strcat('p <', {' '}, '0.001');
        p_value_string = {char(p_value_string)};
        ann_p_value = annotation('textbox', p_value_dims, 'String', p_value_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_p_value.FontSize = legend_size;
        ann_p_value.FontName = font_name;
    end

end % End of target loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, 'vs.', 'vs');
        save_title = strrep(save_title, 'mg.', 'mg');
        save_title = strrep(save_title, 'kg.', 'kg');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        if ~strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), Save_Figs)
        end
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(save_title)), 'fig')
        end
        close gcf
    end
end


