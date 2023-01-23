
function Decomposed_Force_Scatter(xds_morn, xds_noon, norm_force, Save_Figs)

%% Display the function being used

disp('Decomposed Force Scatter Plot:');

%% Basic Settings, some variable extractions, & definitions

% Plot the mean recomposed force (1 = Yes, 0 = No)
recom_Force = 1;
% Plot the mean decomposed force (1 = Yes, 0 = No)
decom_Force = 1;

% How much to enlarge the axis (to make room for the legend)
axis_expansion = 1;

% Font specifications
label_font_size = 17;
legend_font_size = 13;
title_font_size = 14;

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
    
    %% Define the output variables
    if jj == 1

        per_trial_avg_First_Force_morn = struct([]);
        per_trial_avg_Second_Force_morn = struct([]);
        per_trial_avg_First_Force_noon = struct([]);
        per_trial_avg_Second_Force_noon = struct([]);

        First_Morn_TgtHold_Force = zeros(1, num_dirs);
        Second_Morn_TgtHold_Force = zeros(1, num_dirs);
        Morn_TgtHold_Force = zeros(1, num_dirs);

        STD_First_Morn_TgtHold_Force = zeros(1, num_dirs);
        STD_Second_Morn_TgtHold_Force = zeros(1, num_dirs);
        STD_Morn_TgtHold_Force = zeros(1, num_dirs);

        Err_First_Morn_TgtHold_Force = zeros(1, num_dirs);
        Err_Second_Morn_TgtHold_Force = zeros(1, num_dirs);
        Err_Morn_TgtHold_Force = zeros(1, num_dirs);

        First_Noon_TgtHold_Force = zeros(1, num_dirs);
        Second_Noon_TgtHold_Force = zeros(1, num_dirs);
        Noon_TgtHold_Force = zeros(1, num_dirs);

        STD_First_Noon_TgtHold_Force = zeros(1, num_dirs);
        STD_Second_Noon_TgtHold_Force = zeros(1, num_dirs);
        STD_Noon_TgtHold_Force = zeros(1, num_dirs);

        Err_First_Noon_TgtHold_Force = zeros(1, num_dirs);
        Err_Second_Noon_TgtHold_Force = zeros(1, num_dirs);
        Err_Noon_TgtHold_Force = zeros(1, num_dirs);

        TgtHold_First_Force_p_values = zeros(1, num_dirs);
        TgtHold_Second_Force_p_values = zeros(1, num_dirs);
        TgtHold_Force_p_values = zeros(1, num_dirs);

        First_TgtHold_Force_perc_changes = zeros(1, num_dirs);
        Second_TgtHold_Force_perc_changes = zeros(1, num_dirs);
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
    per_trial_avg_First_Force_morn{jj,1} = zeros(length(First_Force_morn), 1);
    per_trial_avg_Second_Force_morn{jj,1} = zeros(length(Second_Force_morn), 1);
    per_trial_avg_First_Force_noon{jj,1} = zeros(length(First_Force_noon), 1);
    per_trial_avg_Second_Force_noon{jj,1} = zeros(length(Second_Force_noon), 1);
    for ii = 1:length(First_Force_morn)
        per_trial_avg_First_Force_morn{jj,1}(ii,1) = mean(all_trials_First_Force_morn(:,ii));
        per_trial_avg_Second_Force_morn{jj,1}(ii,1) = mean(all_trials_Second_Force_morn(:,ii));
    end
    for ii = 1:length(First_Force_noon)
        per_trial_avg_First_Force_noon{jj,1}(ii,1) = mean(all_trials_First_Force_noon(:,ii));
        per_trial_avg_Second_Force_noon{jj,1}(ii,1) = mean(all_trials_Second_Force_noon(:,ii));
    end
    
    %% Calculating average Force (Average Across trials)
    cross_trial_avg_First_Force_morn = zeros(length(First_Force_morn), 1);
    cross_trial_avg_Second_Force_morn = zeros(length(Second_Force_morn), 1);
    cross_trial_avg_First_Force_noon = zeros(length(First_Force_noon), 1);
    cross_trial_avg_Second_Force_noon = zeros(length(Second_Force_noon), 1);
    for ii = 1:length(First_Force_morn{1,1})
        cross_trial_avg_First_Force_morn(ii,1) = mean(all_trials_First_Force_morn(ii,:));
        cross_trial_avg_Second_Force_morn(ii,1) = mean(all_trials_Second_Force_morn(ii,:));
    end
    for ii = 1:length(First_Force_noon{1,1})
        cross_trial_avg_First_Force_noon(ii,1) = mean(all_trials_First_Force_noon(ii,:));
        cross_trial_avg_Second_Force_noon(ii,1) = mean(all_trials_Second_Force_noon(ii,:));
    end

    %% Calculate the output variables
    [~, TgtHold_First_Force_p_values(1, jj)] = ttest2(per_trial_avg_First_Force_morn{jj,1}, per_trial_avg_First_Force_noon{jj,1});
    [~, TgtHold_Second_Force_p_values(1, jj)] = ttest2(per_trial_avg_Second_Force_morn{jj,1}, per_trial_avg_Second_Force_noon{jj,1});

    per_trial_avg_Force_morn = per_trial_avg_First_Force_morn{jj,1} + per_trial_avg_Second_Force_morn{jj,1};
    per_trial_avg_Force_noon = per_trial_avg_First_Force_noon{jj,1} + per_trial_avg_Second_Force_noon{jj,1};
    [~, TgtHold_Force_p_values(1, jj)] = ttest2(per_trial_avg_Force_morn, per_trial_avg_Force_noon);

    % Average Force
    First_Morn_TgtHold_Force(1, jj) = mean(per_trial_avg_First_Force_morn{jj,1});
    Second_Morn_TgtHold_Force(1, jj) = mean(per_trial_avg_Second_Force_morn{jj,1});
    Morn_TgtHold_Force(1, jj) = First_Morn_TgtHold_Force(1, jj) + Second_Morn_TgtHold_Force(1, jj);

    First_Noon_TgtHold_Force(1, jj) = mean(per_trial_avg_First_Force_noon{jj,1});
    Second_Noon_TgtHold_Force(1, jj) = mean(per_trial_avg_Second_Force_noon{jj,1});
    Noon_TgtHold_Force(1, jj) = First_Noon_TgtHold_Force(1, jj) + Second_Noon_TgtHold_Force(1, jj);

    % Standard deviation
    STD_First_Morn_TgtHold_Force(1, jj) = std(per_trial_avg_First_Force_morn{jj,1});
    STD_Second_Morn_TgtHold_Force(1, jj) = std(per_trial_avg_Second_Force_morn{jj,1});
    STD_Morn_TgtHold_Force(1, jj) = std(per_trial_avg_Force_morn);

    STD_First_Noon_TgtHold_Force(1, jj) = std(per_trial_avg_First_Force_noon{jj,1});
    STD_Second_Noon_TgtHold_Force(1, jj) = std(per_trial_avg_Second_Force_noon{jj,1});
    STD_Noon_TgtHold_Force(1, jj) = std(per_trial_avg_Force_noon);

    % Standard error
    Err_First_Morn_TgtHold_Force(1, jj) = STD_First_Morn_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_First_Force_morn{jj,1}));
    Err_Second_Morn_TgtHold_Force(1, jj) = STD_Second_Morn_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_Second_Force_morn{jj,1}));
    Err_Morn_TgtHold_Force(1, jj) = STD_Morn_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_Force_morn));

    Err_First_Noon_TgtHold_Force(1, jj) = STD_First_Noon_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_First_Force_noon{jj,1}));
    Err_Second_Noon_TgtHold_Force(1, jj) = STD_Second_Noon_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_Second_Force_noon{jj,1}));
    Err_Noon_TgtHold_Force(1, jj) = STD_Noon_TgtHold_Force(1, jj) / sqrt(length(per_trial_avg_Force_noon));
    
    % Percent change
    First_TgtHold_Force_perc_changes(1, jj) = (First_Noon_TgtHold_Force - First_Morn_TgtHold_Force) / abs(First_Morn_TgtHold_Force);
    Second_TgtHold_Force_perc_changes(1, jj) = (Second_Noon_TgtHold_Force - Second_Morn_TgtHold_Force) / abs(Second_Morn_TgtHold_Force);

end % End of target loop

%% Plot the force scatter

force_fig = figure;
force_fig.Position = [200 50 800 600];
hold on

% Label the axis
xlabel('Force: Sensor 1', 'FontSize', label_font_size);
ylabel('Force: Sensor 2', 'FontSize', label_font_size);

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
title_string = strcat(Date, {' '}, Task, ',', {' '}, Drug, ': Decomposed Force');
title(title_string, 'FontSize', title_font_size)

% Force marker shape & marker size
decom_marker_metric ='.';
decom_sz = 100;
if ~isequal(decom_Force, 0)
    decom_mean_marker_metric ='s';
    decom_mean_sz = 50;
end
if ~isequal(recom_Force, 0)
    recom_marker_metric ='*';
    recom_sz = 100;
end

for jj = 1:length(per_trial_avg_First_Force_morn)

    % Plot the normalized TgtHold Force
    scatter(per_trial_avg_First_Force_morn{jj}, per_trial_avg_Second_Force_morn{jj}, decom_sz, decom_marker_metric, 'MarkerEdgeColor', ...
        [0.9290, 0.6940, 0.1250], 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
    scatter(per_trial_avg_First_Force_noon{jj}, per_trial_avg_Second_Force_noon{jj}, decom_sz, decom_marker_metric, 'MarkerEdgeColor', ...
        [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5], 'LineWidth', 1.5);

end

for jj = 1:width(TgtHold_First_Force_p_values)

    % If the TgtHold Force change is significant
    if TgtHold_First_Force_p_values(jj) <= 0.05
        First_color_metric = 1;
    else
        First_color_metric = 0;
    end
    if TgtHold_Second_Force_p_values(jj) <= 0.05
        Second_color_metric = 1;
    else
        Second_color_metric = 0;
    end
    if TgtHold_Force_p_values(jj) <= 0.05
        color_metric = 1;
    else
        color_metric = 0;
    end

    % TgtHold Force standard error
    if ~isequal(decom_Force, 0)
        err_First_morn = errorbar(First_Morn_TgtHold_Force(jj), Second_Morn_TgtHold_Force(jj), ...
            Err_First_Morn_TgtHold_Force(jj), 'horizontal');
        err_First_morn.Color = [First_color_metric, 0, 0];
        err_First_morn.LineWidth = 1;
        err_Second_morn = errorbar(First_Morn_TgtHold_Force(jj), Second_Morn_TgtHold_Force(jj), ...
            Err_Second_Morn_TgtHold_Force(jj), 'vertical');
        err_Second_morn.Color = [Second_color_metric, 0, 0];
        err_Second_morn.LineWidth = 1;

        err_First_noon = errorbar(First_Noon_TgtHold_Force(jj), Second_Noon_TgtHold_Force(jj), ...
            Err_First_Noon_TgtHold_Force(jj), 'horizontal');
        err_First_noon.Color = [First_color_metric, 0, 0];
        err_First_noon.LineWidth = 1;
        err_Second_noon = errorbar(First_Noon_TgtHold_Force(jj), Second_Noon_TgtHold_Force(jj), ...
            Err_Second_Noon_TgtHold_Force(jj), 'vertical');
        err_Second_noon.Color = [Second_color_metric, 0, 0];
        err_Second_noon.LineWidth = 1;
    end

    if ~isequal(recom_Force, 0)
        err_morn = errorbar(Morn_TgtHold_Force(jj) / 2, Noon_TgtHold_Force(jj) / 2, ...
            Err_Morn_TgtHold_Force(jj) / 2, 'horizontal');
        err_morn.Color = [color_metric, 0, 0];
        err_morn.LineWidth = 1;
        err_noon = errorbar(Morn_TgtHold_Force(jj) / 2, Noon_TgtHold_Force(jj) / 2, ...
            Err_Noon_TgtHold_Force(jj) / 2, 'vertical');
        err_noon.Color = [color_metric, 0, 0];
        err_noon.LineWidth = 1;
    end

    % Plot the normalized TgtHold Force
    if ~isequal(decom_Force, 0)
        scatter(First_Morn_TgtHold_Force(jj), Second_Morn_TgtHold_Force(jj), decom_mean_sz, decom_mean_marker_metric, 'MarkerEdgeColor', ...
            [0.9290, 0.6940, 0.1250], 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
        scatter(First_Noon_TgtHold_Force(jj), Second_Noon_TgtHold_Force(jj), decom_mean_sz, decom_mean_marker_metric, 'MarkerEdgeColor', ...
            [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5], 'LineWidth', 1.5);
    end
    if ~isequal(recom_Force, 0)
        scatter(Morn_TgtHold_Force(jj) / 2, Noon_TgtHold_Force(jj) / 2, recom_sz, recom_marker_metric, 'MarkerEdgeColor', ...
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);
    end
end

% Calculate the axis limits
curr_axis = gca;
min_x = curr_axis.XLim(1);
min_y = curr_axis.YLim(1);
axis_min = round(min(min_x, min_y)/5)*5 - axis_expansion;
max_x = curr_axis.XLim(2);
max_y = curr_axis.YLim(2);
axis_max = round(max(max_x, max_y)/5)*5 + axis_expansion;

% Draw the identity line 
line([axis_min, axis_max],[axis_min, axis_max], ...
    'Color', 'k', 'Linewidth', 1, 'Linestyle','--')

% Plot dummy points for the top left legend
dummy_circle = scatter(-47, -47, 75, 'o', 'MarkerEdgeColor',[0 0 0]);
dummy_square = scatter(-50, -50, 75, 's', 'MarkerEdgeColor',[0 0 0]);
dummy_star = scatter(-55, -55, 75, '*', 'MarkerEdgeColor',[0 0 0]);
% Plot dummy points for the bottom right legend
dummy_yellow = plot(-40,-40, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor', [0.9290, 0.6940, 0.1250], 'MarkerFaceColor', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
dummy_purple = plot(-45,-45, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor', [0.5 0 0.5], 'MarkerFaceColor', [0.5 0 0.5], 'LineWidth', 1.5);
if ~isequal(recom_Force, 0) && ~isequal(decom_Force, 0)
    dummy_red = plot(-40,-40, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor',[1, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
    dummy_black = plot(-45,-45, '.', 'MarkerSize',20, ...
        'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);
end

% Set the axis
xlim([axis_min, axis_max])
ylim([axis_min, axis_max])

% Record the title
fig_info = get(gca,'title');
save_title = get(fig_info, 'string');

% Plot the bottom right legend
if ~isequal(recom_Force, 0) || ~isequal(decom_Force, 0)
    right_legend = legend([dummy_yellow, dummy_purple, dummy_red, dummy_black], {'Morning', 'Afternoon', 'p <= 0.05', 'p > 0.05'}, ...
        'FontSize', legend_font_size, 'Location', 'SouthEast');
else
    right_legend = legend([dummy_yellow, dummy_purple], {'Morning', 'Afternoon'}, ...
        'FontSize', legend_font_size, 'Location', 'SouthEast');
end
right_legend.ItemTokenSize(1) = 15;
legend boxoff

% Plot the top top left legend
if ~isequal(recom_Force, 0) || ~isequal(decom_Force, 0)
    legend_axes = axes('position', get(gca,'position'), 'visible', 'off');
end
if ~isequal(recom_Force, 0) && ~isequal(decom_Force, 0)
    left_legend = legend(legend_axes, [dummy_circle, dummy_square, dummy_star], ...
        {'Single Trials', 'Mean Force', 'Recomposed Mean Force'}, ...
        'FontSize', legend_font_size, 'Location', 'northwest');
elseif isequal(recom_Force, 0) && ~isequal(decom_Force, 0)
    left_legend = legend(legend_axes, [dummy_circle, dummy_square], ...
        {'Single Trials', 'Mean Force'}, 'FontSize', legend_font_size, 'Location', 'northwest');
elseif ~isequal(recom_Force, 0) && isequal(decom_Force, 0)
    left_legend = legend(legend_axes, [dummy_circle, dummy_star], ...
        {'Single Trials', 'Recomposed Mean Force'}, ...
        'FontSize', legend_font_size, 'Location', 'northwest');
end
left_legend.ItemTokenSize(1) = 15;
legend boxoff

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


