%% Loading the morning and afternoon files
%clear
clearvars -except xds & xds_morn & xds_noon & unit_name & event & Save_Figs & ... 
    Plot_Figs & muscle_groups & zero_EMG & norm_EMG
clc

% Select The Date & Task To Analyze
Date = '20211001';
Task = 'WS';
% Do You Want To Process The XDS File? (1 = Yes; 0 = No)
Process_XDS = 1;

[xds_morn, xds_noon, xds_excel] = Load_XDS(Date, Task, Process_XDS);

% Which targets do you want the mnovement phase firing rate calculated from? ('Max', 'Min', 'All')
tgt_mpfr = 'Max';

%% Basic Settings, some variable extractions, & definitions

% Save the figures to desktop? ('pdf', 'png', 'fig', 0 = No)
Save_Figs = 0;
if ~isequal(Save_Figs, 0)
    close all
end

% Normalization percentile
norm_perctile = 95;
% Normalize Force? (1 = Yes, 0 = No)
norm_force = 1;
% Normalize Cursor? (1 = Yes, 0 = No)
norm_cursor = 1;
% Normalize EMG? (1 = Yes, 0 = No)
norm_EMG = 1;

% How much to enlarge the axis (to make room for the legend)
axis_expansion = 0;

% Font specifications
label_font_size = 17;
title_font_size = 14;
legend_font_size = 13;
font_name = 'Arial';

% Save Counter
ss = 1;

%% Calculate the number of target directions 

% Find the number of targets in the morning and afternoon
[dir_idx_morn, unique_tgts_morn, ~, ~, ~, ~] = ...
    WindowTrialGoCueFiringRate(xds_morn, 1, tgt_mpfr);
[dir_idx_noon, unique_tgts_noon, ~, ~, ~, ~] = ...
    WindowTrialGoCueFiringRate(xds_noon, 1, tgt_mpfr);

% Check to see if both sessions use a consistent number of directions
if ~isequal(unique(dir_idx_morn), unique(dir_idx_noon))
    disp('Uneven Target Directions Between Morning & Afternoon');
    % Only use the info of target directions conserved between morn & noon
    shared_dir_idx = ismember(dir_idx_morn, dir_idx_noon);
    dir_idx_morn = dir_idx_morn(shared_dir_idx);
    unique_tgts_morn = unique_tgts_morn(shared_dir_idx);
    dir_idx_noon = dir_idx_noon(shared_dir_idx);
    unique_tgts_noon = unique_tgts_noon(shared_dir_idx);
end

% Check to see if both sessions use a consistent number of targets
if ~isequal(unique(unique_tgts_morn), unique(unique_tgts_noon))
    disp('Uneven Target Centers Between Morning & Afternoon');
    % Only use the info of target centers conserved between morn & noon
    shared_target_centers_idx = ismember(unique_tgts_morn, unique_tgts_noon);
    dir_idx_morn = dir_idx_morn(shared_target_centers_idx);
    unique_tgts_morn = unique_tgts_morn(shared_target_centers_idx);
    dir_idx_noon = dir_idx_noon(shared_target_centers_idx);
    unique_tgts_noon = unique_tgts_noon(shared_target_centers_idx);
end

%% EMG

if strcmp(Task, 'PG')
    % All, Flex, Exten, Uln_Dev, Rad_Dev, Grasp, or Custom?
    muscle_groups = 'Grasp';
elseif strcmp(Task, 'WS')
    if length(unique_tgts_morn) > 3
        muscle_groups = 'Both';
    end
end

% Zero? (1 = Yes, 0 = No)
zero_EMG = 1;
zero_method = 'Percentile';
EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);

EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_perctile, norm_EMG);

% EMG Statistics
[EMG_Names, Baseline_EMG_p_values, Baseline_EMG_perc_changes, ...
    Morn_Baseline_EMG, Err_Morn_Baseline_EMG, Noon_Baseline_EMG, Err_Noon_Baseline_EMG] = ...
    Baseline_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
    EMG_Zero_Factor, EMG_Norm_Factor, 0, 0);
[~, TgtHold_EMG_p_values, TgtHold_EMG_perc_changes, ...
    Morn_TgtHold_EMG, Err_Morn_TgtHold_EMG, Noon_TgtHold_EMG, Err_Noon_TgtHold_EMG] = ...
    TgtHold_EMG_Stats(xds_morn, xds_noon, muscle_groups, ... 
    EMG_Zero_Factor, EMG_Norm_Factor, 0, 0);

%% Force
if strcmp(Task, 'PG')
    Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_perctile, norm_force);

    % Force Statistics
    [TgtHold_p_values, TgtHold_perc_changes, ...
        Morn_TgtHold, Err_Morn_TgtHold, Noon_TgtHold, Err_Noon_TgtHold] = ...
        TgtHold_Force_Stats(xds_morn, xds_noon, Force_Norm_Factor, 0, 0);
end

%% Cursor
if strcmp(Task, 'WS')
    Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_perctile, norm_cursor);

    % Cursor Position Statistics
    [Baseline_Cursor_Pos_p_values, Baseline_Cursor_Pos_perc_changes, ...
    Morn_Baseline_CursorPos, Err_Morn_Baseline_CursorPos, Noon_Baseline_CursorPos, Err_Noon_Baseline_CursorPos] = ...
        Baseline_CursorPos_Stats(xds_morn, xds_noon, Cursor_Norm_Factor, 0, 0);
    [TgtHold_p_values, TgtHold_perc_changes, ...
        Morn_TgtHold, Err_Morn_TgtHold, Noon_TgtHold, Err_Noon_TgtHold] = ...
        TgtHold_CursorPos_Stats(xds_morn, xds_noon, Cursor_Norm_Factor, 0, 0);
end

%% Plot the Baseline Behavioral Stats Scatter

figure
hold on

% Label the axis
xlabel('Morning Normalized Baseline EMG', 'FontSize', label_font_size);
ylabel('Afternoon Normalized Baseline EMG', 'FontSize', label_font_size);

% Set the title
title_string = strcat(Date, {' '}, Task, ': Baseline Behavioral Metrics');
title(title_string, 'FontSize', title_font_size)

% Baseline EMG marker shape & marker size
marker_metric ='.';
sz = 500;

for jj = 1:length(EMG_Names)
    for kk = 1:width(Baseline_EMG_p_values)

        % If the baseline EMG change is significant
        if Baseline_EMG_p_values(jj, kk) <= 0.05
            color_metric = 1;
        else
            color_metric = 0;
        end

        % Plot the normalized baseline EMG
        scatter(Morn_Baseline_EMG(jj, kk), Noon_Baseline_EMG(jj, kk), sz, marker_metric, 'MarkerEdgeColor', ...
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

        % Baseline EMG standard deviation
        err_morn = errorbar(Morn_Baseline_EMG(jj, kk), Noon_Baseline_EMG(jj, kk), ...
            Err_Morn_Baseline_EMG(jj, kk), 'horizontal');
        err_morn.Color = [color_metric, 0, 0];
        err_morn.LineWidth = 1;
        err_noon = errorbar(Morn_Baseline_EMG(jj, kk), Noon_Baseline_EMG(jj, kk), ...
            Err_Noon_Baseline_EMG(jj, kk), 'vertical');
        err_noon.Color = [color_metric, 0, 0];
        err_noon.LineWidth = 1;

    end
end

if strcmp(Task, 'WS')
    % Baseline Cursor marker shape & marker size
    marker_metric ='*';
    sz = 100;

    for jj = 1:width(Baseline_Cursor_Pos_p_values)

        % If the Baseline Cursor change is significant
        if Baseline_Cursor_Pos_p_values(jj) <= 0.05
            color_metric = 1;
        else
            color_metric = 0;
        end

        % Plot the normalized Baseline Cursor
        scatter(Morn_Baseline_CursorPos(jj), Noon_Baseline_CursorPos(jj), sz, marker_metric, 'MarkerEdgeColor', ...
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

        % Baseline Cursor standard error
        err_morn = errorbar(Morn_Baseline_CursorPos(jj), Noon_Baseline_CursorPos(jj), ...
            Err_Morn_Baseline_CursorPos(jj), 'horizontal');
        err_morn.Color = [color_metric, 0, 0];
        err_morn.LineWidth = 1;
        err_noon = errorbar(Morn_Baseline_CursorPos(jj), Noon_Baseline_CursorPos(jj), ...
            Err_Noon_Baseline_CursorPos(jj), 'vertical');
        err_noon.Color = [color_metric, 0, 0];
        err_noon.LineWidth = 1;
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
dummy_star = scatter(-55, -55, 75, '*', 'MarkerEdgeColor',[0 0 0]);
% Plot dummy points for the bottom right legend
dummy_red = plot(-40,-40, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[1, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
dummy_black = plot(-45,-45, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);

% Set the axis
xlim([axis_min, axis_max])
ylim([axis_min, axis_max])

% Set ticks to outside
figure_axes = gca;
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')
% Set The Font
set(figure_axes,'FontName', font_name);

% Record the title
fig_info = get(gca,'title');
save_title(ss) = get(fig_info, 'string');
% Add to the counter
ss = ss + 1;

% Plot the bottom right legend
right_legend = legend([dummy_red, dummy_black], {'p <= 0.05', 'p > 0.05'}, ...
    'FontSize', legend_font_size, 'Location', 'southeast');
right_legend.ItemTokenSize(1) = 15;
legend boxoff
% Plot the top top left legend
legend_axes = axes('position', get(gca,'position'), 'visible', 'off');

if strcmp(Task, 'PG')
    left_legend = legend(legend_axes, (dummy_circle), {'Baseline EMG'}, ...
        'FontSize', legend_font_size, 'Location', 'northwest');
    elseif strcmp(Task, 'WS')
        left_legend = legend(legend_axes, [dummy_circle, dummy_star], ...
            {'Baseline EMG', 'Baseline Wrist Position'}, ...
            'FontSize', legend_font_size, 'Location', 'northwest');
end

left_legend.ItemTokenSize(1) = 15;
legend boxoff

%% Plot the TgtHold Behavioral Stats Scatter

figure
hold on

% Label the axis
xlabel('Morning Normalized TgtHold EMG', 'FontSize', label_font_size);
ylabel('Afternoon Normalized TgtHold EMG', 'FontSize', label_font_size);

% Set the title
title_string = strcat(Date, {' '}, Task, ': TgtHold Behavioral Metrics');
title(title_string, 'FontSize', title_font_size)

% TgtHold EMG marker shape & marker size
marker_metric ='.';
sz = 500;

for jj = 1:length(EMG_Names)
    for kk = 1:width(TgtHold_EMG_p_values)

        % If the TgtHold EMG change is significant
        if TgtHold_EMG_p_values(jj, kk) <= 0.05
            color_metric = 1;
        else
            color_metric = 0;
        end

        % Plot the normalized TgtHold EMG
        scatter(Morn_TgtHold_EMG(jj, kk), Noon_TgtHold_EMG(jj, kk), sz, marker_metric, 'MarkerEdgeColor', ...
            [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

        % TgtHold EMG standard deviation
        err_morn = errorbar(Morn_TgtHold_EMG(jj, kk), Noon_TgtHold_EMG(jj, kk), ...
            Err_Morn_TgtHold_EMG(jj, kk), 'horizontal');
        err_morn.Color = [color_metric, 0, 0];
        err_morn.LineWidth = 1;
        err_noon = errorbar(Morn_TgtHold_EMG(jj, kk), Noon_TgtHold_EMG(jj, kk), ...
            Err_Noon_TgtHold_EMG(jj, kk), 'vertical');
        err_noon.Color = [color_metric, 0, 0];
        err_noon.LineWidth = 1;

    end
end

% Force / Cursror marker shape & marker size
marker_metric ='*';
sz = 100;

for jj = 1:width(TgtHold_p_values)

    % If the TgtHold Force / Cursror change is significant
    if TgtHold_p_values(jj) <= 0.05
        color_metric = 1;
    else
        color_metric = 0;
    end

    % Plot the normalized TgtHold Force / Cursror
    scatter(Morn_TgtHold(jj), Noon_TgtHold(jj), sz, marker_metric, 'MarkerEdgeColor', ...
        [color_metric 0 0], 'MarkerFaceColor', [color_metric 0 0], 'LineWidth', 1.5);

    % TgtHold Force / Cursror standard error
    err_morn = errorbar(Morn_TgtHold(jj), Noon_TgtHold(jj), ...
        Err_Morn_TgtHold(jj), 'horizontal');
    err_morn.Color = [color_metric, 0, 0];
    err_morn.LineWidth = 1;
    err_noon = errorbar(Morn_TgtHold(jj), Noon_TgtHold(jj), ...
        Err_Noon_TgtHold(jj), 'vertical');
    err_noon.Color = [color_metric, 0, 0];
    err_noon.LineWidth = 1;

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
dummy_star = scatter(-55, -55, 75, '*', 'MarkerEdgeColor',[0 0 0]);
% Plot dummy points for the bottom right legend
dummy_red = plot(-40,-40, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[1, 0, 0], 'MarkerFaceColor',[1, 0, 0], 'LineWidth', 1.5);
dummy_black = plot(-45,-45, '.', 'MarkerSize',20, ...
    'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0, 0, 0], 'LineWidth', 1.5);

% Set the axis
xlim([axis_min, axis_max])
ylim([axis_min, axis_max])

% Set ticks to outside
figure_axes = gca;
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')
% Set The Font
set(figure_axes,'FontName', font_name);

% Record the title
fig_info = get(gca,'title');
save_title(ss) = get(fig_info, 'string');
% Add to the counter
ss = ss + 1;

if strcmp(Task, 'WS')
    force_curs = 'TgtHold Wrist Position';
elseif strcmp(Task, 'PG')
    force_curs = 'TgtHold Force';
end

% Plot the bottom right legend
right_legend = legend([dummy_red, dummy_black], {'p <= 0.05', 'p > 0.05'}, ...
    'FontSize', legend_font_size, 'Location', 'southeast');
right_legend.ItemTokenSize(1) = 15;
legend boxoff
% Plot the top top left legend
legend_axes = axes('position', get(gca,'position'), 'visible', 'off');

left_legend = legend(legend_axes, [dummy_circle, dummy_star], ...
    {'TgtHold EMG', force_curs}, ...
    'FontSize', legend_font_size, 'Location', 'northwest');

left_legend.ItemTokenSize(1) = 15;
legend boxoff

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = numel(findobj('type','figure')):-1:1
        save_title(ii) = strrep(save_title(ii), ':', '');
        save_title(ii) = strrep(save_title(ii), '.0', '');
        save_title(ii) = strrep(save_title(ii), 'vs.', 'vs');
        save_title(ii) = strrep(save_title(ii), 'mg.', 'mg');
        save_title(ii) = strrep(save_title(ii), 'kg.', 'kg');
        save_title(ii) = strrep(save_title(ii), '.', '_');
        save_title(ii) = strrep(save_title(ii), '/', '_');
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


