function Target_Force_Traces(xds_morn, xds_noon, norm_force, Save_Figs)

%% Display the function being used

disp('Target Force Traces:');

%% Basic Settings, some variable extractions, & definitions

% Find the force y-limits
force_YLims = ForceYLimit(xds_morn, xds_noon, norm_force);

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

% Font specifications
label_font_size = 17;
legend_font_size = 13;
title_font_size = 14;
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
    [rewarded_gocue_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_gocue');
    [rewarded_end_time_morn] = EventAlignmentTimes(xds_morn, target_dirs_morn(jj), target_centers_morn(jj), 'trial_end');
    [rewarded_gocue_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_gocue');
    [rewarded_end_time_noon] = EventAlignmentTimes(xds_noon, target_dirs_noon(jj), target_centers_noon(jj), 'trial_end');
    
    %% Define the output variables
    if jj == 1
    
        per_trial_avg_First_Force_morn = struct([]);
        per_trial_avg_Second_Force_morn = struct([]);
        per_trial_avg_First_Force_noon = struct([]);
        per_trial_avg_Second_Force_noon = struct([]);
    
        all_trials_First_Force_morn = struct([]);
        all_trials_Second_Force_morn = struct([]);
        all_trials_First_Force_noon = struct([]);
        all_trials_Second_Force_noon = struct([]);
    
        First_Trial_Force_morn = struct([]);
        Second_Trial_Force_morn = struct([]);
        First_Trial_Force_noon = struct([]);
        Second_Trial_Force_noon = struct([]);
    
        TgtHold_First_Force_p_values = zeros(1, num_dirs);
        TgtHold_Second_Force_p_values = zeros(1, num_dirs);

    end

    %% Force and time aligned to specified event
    % Find the rewarded times in the whole trial time frame
    rewarded_gocue_idx_morn = zeros(height(rewarded_gocue_time_morn),1);
    for ii = 1:length(rewarded_gocue_time_morn)
        rewarded_gocue_idx_morn(ii) = find(xds_morn.time_frame == rewarded_gocue_time_morn(ii));
    end

    rewarded_end_idx_morn = zeros(height(rewarded_gocue_time_morn),1);
    for ii = 1:length(rewarded_gocue_time_morn)
        rewarded_end_idx_morn(ii) = find(xds_morn.time_frame == rewarded_end_time_morn(ii));
    end

    rewarded_gocue_idx_noon = zeros(height(rewarded_gocue_time_noon),1);
    for ii = 1:length(rewarded_gocue_time_noon)
        rewarded_gocue_idx_noon(ii) = find(xds_noon.time_frame == rewarded_gocue_time_noon(ii));
    end

    rewarded_end_idx_noon = zeros(height(rewarded_gocue_time_noon),1);
    for ii = 1:length(rewarded_gocue_time_noon)
        rewarded_end_idx_noon(ii) = find(xds_noon.time_frame == rewarded_end_time_noon(ii));
    end

    % Force during each successful trial from gocue to target hold
    for ii = 1:length(rewarded_gocue_time_morn)
        First_Trial_Force_morn{jj,1}{ii, 1} = xds_morn.force(rewarded_gocue_idx_morn(ii) : ...
            (rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width)), 1);
        Second_Trial_Force_morn{jj,1}{ii, 1} = xds_morn.force(rewarded_gocue_idx_morn(ii) : ...
            (rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width)), 2);
    end
    
    % Force during each successful trial from target hold to reward
    First_TgtHold_Force_morn = struct([]);
    Second_TgtHold_Force_morn = struct([]);
    for ii = 1:length(rewarded_gocue_time_morn)
        First_TgtHold_Force_morn{ii, 1} = xds_morn.force((rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width) : ...
            rewarded_end_idx_morn(ii)), 1);
        Second_TgtHold_Force_morn{ii, 1} = xds_morn.force((rewarded_end_idx_morn(ii) - (TgtHold_time / xds_morn.bin_width) : ...
            rewarded_end_idx_morn(ii)), 2);
    end

    % Force during each successful trial from gocue to target hold
    for ii = 1:length(rewarded_gocue_time_noon)
        First_Trial_Force_noon{jj,1}{ii, 1} = xds_noon.force(rewarded_gocue_idx_noon(ii) : ...
            (rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width)), 1);
        Second_Trial_Force_noon{jj,1}{ii, 1} = xds_noon.force(rewarded_gocue_idx_noon(ii) : ...
            (rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width)), 2);
    end

    % Force during each successful trial from target hold to reward
    First_TgtHold_Force_noon = struct([]);
    Second_TgtHold_Force_noon = struct([]);
    for ii = 1:length(rewarded_gocue_time_noon)
        First_TgtHold_Force_noon{ii, 1} = xds_noon.force((rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width) : ...
            rewarded_end_idx_noon(ii)), 1);
        Second_TgtHold_Force_noon{ii, 1} = xds_noon.force((rewarded_end_idx_noon(ii) - (TgtHold_time / xds_noon.bin_width) : ...
            rewarded_end_idx_noon(ii)), 2);
    end

    %% Putting all succesful trials in one array
    
    all_trials_First_Force_morn{jj,1} = zeros(length(First_TgtHold_Force_morn{1,1}), length(rewarded_gocue_time_morn));
    all_trials_Second_Force_morn{jj,1} = zeros(length(Second_TgtHold_Force_morn{1,1}), length(rewarded_gocue_time_morn));
    for ii = 1:length(rewarded_gocue_time_morn)
        all_trials_First_Force_morn{jj,1}(:,ii) = First_TgtHold_Force_morn{ii, 1};
        all_trials_Second_Force_morn{jj,1}(:,ii) = Second_TgtHold_Force_morn{ii, 1};
    end

    all_trials_First_Force_noon{jj,1} = zeros(length(First_TgtHold_Force_noon{1,1}), length(rewarded_gocue_time_noon));
    all_trials_Second_Force_noon{jj,1} = zeros(length(Second_TgtHold_Force_noon{1,1}), length(rewarded_gocue_time_noon));
    for ii = 1:length(rewarded_gocue_time_noon)
        all_trials_First_Force_noon{jj,1}(:,ii) = First_TgtHold_Force_noon{ii, 1};
        all_trials_Second_Force_noon{jj,1}(:,ii) = Second_TgtHold_Force_noon{ii, 1};
    end

    %% Normalizing the average force
    if isequal(norm_force, 1)
        Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force);
        for pp = 1:length(First_Trial_Force_morn{jj,1})
            First_Trial_Force_morn{jj,1}{pp,1} = First_Trial_Force_morn{jj,1}{pp,1} / Force_Norm_Factor*100;
            Second_Trial_Force_morn{jj,1}{pp,1} = Second_Trial_Force_morn{jj,1}{pp,1} / Force_Norm_Factor*100;
        end
        for pp = 1:length(First_Trial_Force_noon{jj,1})
            First_Trial_Force_noon{jj,1}{pp,1} = First_Trial_Force_noon{jj,1}{pp,1} / Force_Norm_Factor*100;
            Second_Trial_Force_noon{jj,1}{pp,1} = Second_Trial_Force_noon{jj,1}{pp,1} / Force_Norm_Factor*100;
        end
        all_trials_First_Force_morn{jj,1} = all_trials_First_Force_morn{jj,1} / Force_Norm_Factor*100;
        all_trials_Second_Force_morn{jj,1} = all_trials_Second_Force_morn{jj,1} / Force_Norm_Factor*100;
        all_trials_First_Force_noon{jj,1} = all_trials_First_Force_noon{jj,1} / Force_Norm_Factor*100;
        all_trials_Second_Force_noon{jj,1} = all_trials_Second_Force_noon{jj,1} / Force_Norm_Factor*100;
    elseif strcmp(norm_force, 'Convert')
        for pp = 1:length(First_Trial_Force_morn{jj,1})
            First_Trial_Force_morn{jj,1}{pp,1} = First_Trial_Force_morn{jj,1}{pp,1} / 1000*5; % Millivolt conversion * gain
            Second_Trial_Force_morn{jj,1}{pp,1} = Second_Trial_Force_morn{jj,1}{pp,1} / 1000*5; % Millivolt conversion * gain
        end
        for pp = 1:length(First_Trial_Force_noon{jj,1})
            First_Trial_Force_noon{jj,1}{pp,1} = First_Trial_Force_noon{jj,1}{pp,1} / 1000*5; % Millivolt conversion * gain
            Second_Trial_Force_noon{jj,1}{pp,1} = Second_Trial_Force_noon{jj,1}{pp,1} / 1000*5; % Millivolt conversion * gain
        end
        all_trials_First_Force_morn{jj,1} = all_trials_First_Force_morn{jj,1} / 1000*5; % Millivolt conversion * gain
        all_trials_Second_Force_morn{jj,1} = all_trials_Second_Force_morn{jj,1} / 1000*5; % Millivolt conversion * gain
        all_trials_First_Force_noon{jj,1} = all_trials_First_Force_noon{jj,1} / 1000*5; % Millivolt conversion * gain
        all_trials_Second_Force_noon{jj,1} = all_trials_Second_Force_noon{jj,1} / 1000*5; % Millivolt conversion * gain
    end

    %% Calculating average Force (Average per trial)
    per_trial_avg_First_Force_morn{jj,1} = zeros(length(First_TgtHold_Force_morn), 1);
    per_trial_avg_Second_Force_morn{jj,1} = zeros(length(Second_TgtHold_Force_morn), 1);
    per_trial_avg_First_Force_noon{jj,1} = zeros(length(First_TgtHold_Force_noon), 1);
    per_trial_avg_Second_Force_noon{jj,1} = zeros(length(Second_TgtHold_Force_noon), 1);
    for ii = 1:length(First_TgtHold_Force_morn)
        per_trial_avg_First_Force_morn{jj,1}(ii,1) = mean(all_trials_First_Force_morn{jj,1}(:,ii));
        per_trial_avg_Second_Force_morn{jj,1}(ii,1) = mean(all_trials_Second_Force_morn{jj,1}(:,ii));
    end
    for ii = 1:length(First_TgtHold_Force_noon)
        per_trial_avg_First_Force_noon{jj,1}(ii,1) = mean(all_trials_First_Force_noon{jj,1}(:,ii));
        per_trial_avg_Second_Force_noon{jj,1}(ii,1) = mean(all_trials_Second_Force_noon{jj,1}(:,ii));
    end
    
    %% Calculating average Force (Average Across trials)
    cross_trial_avg_First_Force_morn = zeros(length(First_TgtHold_Force_morn), 1);
    cross_trial_avg_Second_Force_morn = zeros(length(Second_TgtHold_Force_morn), 1);
    cross_trial_avg_First_Force_noon = zeros(length(First_TgtHold_Force_noon), 1);
    cross_trial_avg_Second_Force_noon = zeros(length(Second_TgtHold_Force_noon), 1);
    for ii = 1:length(First_TgtHold_Force_morn{1,1})
        cross_trial_avg_First_Force_morn(ii,1) = mean(all_trials_First_Force_morn{jj,1}(ii,:));
        cross_trial_avg_Second_Force_morn(ii,1) = mean(all_trials_Second_Force_morn{jj,1}(ii,:));
    end
    for ii = 1:length(First_TgtHold_Force_noon{1,1})
        cross_trial_avg_First_Force_noon(ii,1) = mean(all_trials_First_Force_noon{jj,1}(ii,:));
        cross_trial_avg_Second_Force_noon(ii,1) = mean(all_trials_Second_Force_noon{jj,1}(ii,:));
    end

    %% Calculate the output variables
    [~, TgtHold_First_Force_p_values(1, jj)] = ttest2(per_trial_avg_First_Force_morn{jj,1}, per_trial_avg_First_Force_noon{jj,1});
    [~, TgtHold_Second_Force_p_values(1, jj)] = ttest2(per_trial_avg_Second_Force_morn{jj,1}, per_trial_avg_Second_Force_noon{jj,1});
    
end % End of target loop

%% Plot the force traces for each target

for kk = 1:length(target_centers_morn)

    figure
    hold on
    
    % Morning Force before the Target Zone
    for jj = 1:width(all_trials_First_Force_morn{kk,1})
        plot(First_Trial_Force_morn{kk,1}{jj,1} - Second_Trial_Force_morn{kk,1}{jj,1}, ...
            First_Trial_Force_morn{kk,1}{jj,1} + Second_Trial_Force_morn{kk,1}{jj,1}, ...
            'LineWidth', 0.5, 'LineStyle', '--', 'Color', [0.9290, 0.6940, 0.1250])
    end
    % Afternoon Force before the Target Zone
    for jj = 1:width(all_trials_First_Force_noon{kk,1})
        plot(First_Trial_Force_noon{kk,1}{jj,1} - Second_Trial_Force_noon{kk,1}{jj,1}, ...
            First_Trial_Force_noon{kk,1}{jj,1} + Second_Trial_Force_noon{kk,1}{jj,1}, ...
            'LineWidth', 0.5, 'LineStyle', '--', 'Color', [0.5 0 0.5])
    end
    % Morning Force in the Target Zone
    for jj = 1:width(all_trials_First_Force_morn{kk,1})
        plot(all_trials_First_Force_morn{kk,1}(:,jj) - all_trials_Second_Force_morn{kk,1}(:,jj), ...
            all_trials_First_Force_morn{kk,1}(:,jj) + all_trials_Second_Force_morn{kk,1}(:,jj), ...
            'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
    end
    % Afternoon Force in the Target Zone
    for jj = 1:width(all_trials_First_Force_noon{kk,1})
        plot(all_trials_First_Force_noon{kk,1}(:,jj) - all_trials_Second_Force_noon{kk,1}(:,jj), ...
            all_trials_First_Force_noon{kk,1}(:,jj) + all_trials_Second_Force_noon{kk,1}(:,jj), ...
            'LineWidth', 2, 'Color', [0.5 0 0.5])
    end

    % Set the title
    title_string = strcat(Date, {' '}, Task, ',', {' '}, Drug, ': TgtCenter at', {' '}, num2str(target_centers_morn(kk)));
    title(title_string, 'FontSize', title_font_size)

    % Label the axis
    xlabel('FS1 - FS2', 'FontSize', label_font_size);
    ylabel('FS1 + FS2', 'FontSize', label_font_size);

    % Set the axis
    %xlim([-90, 90])
    ylim([force_YLims(2), force_YLims(1)]);

    % Plot dummy points for the bottom right legend
    dummy_yellow = plot(-40,-40, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [0.9290, 0.6940, 0.1250], 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
    dummy_purple = plot(-45,-45, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor', [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5], 'LineWidth', 1.5);

    % Plot the bottom bottom right legend
    right_legend = legend([dummy_yellow, dummy_purple], {'Morning', 'Afternoon'}, ...
        'FontSize', legend_font_size, 'Location', 'SouthEast');
    right_legend.ItemTokenSize(1) = 15;
    legend boxoff

     % Annotation of the FS1 p-value
     if round(TgtHold_First_Force_p_values(kk), 3) > 0
        legend_dims = [0.05 0 0.4 0.275];
        p_value_string = strcat('FS1: p =', {' '}, mat2str(round(TgtHold_First_Force_p_values(kk), 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end
    if isequal(round(TgtHold_First_Force_p_values(kk), 3), 0)
        legend_dims = [0.05 0 0.4 0.275];
        p_value_string = strcat('FS1: p <', {' '}, '0.001');
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end
    % Annotation of the FS2 p-value
     if round(TgtHold_Second_Force_p_values(kk), 3) > 0
        legend_dims = [0.05 0 0.4 0.225];
        p_value_string = strcat('FS2: p =', {' '}, mat2str(round(TgtHold_Second_Force_p_values(kk), 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end
    if isequal(round(TgtHold_Second_Force_p_values(kk), 3), 0)
        legend_dims = [0.05 0 0.4 0.225];
        p_value_string = strcat('FS2: p <', {' '}, '0.001');
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:numel(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
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


