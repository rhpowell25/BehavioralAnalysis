function Rxn_Time_BoxPlot(xds_morn, xds_noon, Save_Figs)

%% File Description:

% This function plots a box plot of the reaction time (defined as the time after the
% go-cue when the EMG of interest exceeds 2 std of the baseline EMG).
% The EMG of interest is chosen based on the task / target.
% This function can also plot a box plot of the trial lengths for each
% succesful trial in each unique target direction and target center combination
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% Save_Figs: 'pdf', 'png', 'fig', 'All', or 0

%% Display the function being used
disp('Reaction Time Box Plot Function:')

%% Basic settings, some variable extractions, & definitions

% Which time metric? 'Rxn_Time' or 'Trial_Length'
Time_Metric = 'Rxn_Time';

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

% Close all previously open figures if you're saving 
if ~isequal(Save_Figs, 0)
    close all
end

%% Get the reaction times & trial variables
disp('Morning:')
[target_dirs_morn, target_centers_morn, rxn_time_morn, trial_length_morn] = Trial_Length_And_Rxn_Time(xds_morn);
disp('Afternoon:')
[target_dirs_noon, target_centers_noon, rxn_time_noon, trial_length_noon] = Trial_Length_And_Rxn_Time(xds_noon);

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
    rxn_time_morn = rxn_time_morn(Matching_Idxs_Morn);
    rxn_time_noon = rxn_time_noon(Matching_Idxs_Noon);
    trial_length_morn = trial_length_morn(Matching_Idxs_Morn);
    trial_length_noon = trial_length_noon(Matching_Idxs_Noon);
end

%% Put the trial times in the same struct of the plots

box_plot = struct([]);
for pp = 1:length(target_dirs_morn)
    % Reaction time
    if strcmp(Time_Metric, 'Rxn_Time')
        box_plot{pp,1} = [rxn_time_morn{pp,1}; rxn_time_noon{pp,1}];
    % Trial length
    elseif strcmp(Time_Metric, 'Trial_Length')
        box_plot{pp,1} = [trial_length_morn{pp,1}; trial_length_noon{pp,1}];
    end
end

%% Do the statistics

box_plot_p_val = zeros(length(target_dirs_morn),1);
for pp = 1:length(target_dirs_morn)
    % Reaction time
    if strcmp(Time_Metric, 'Rxn_Time')
        [~, box_plot_p_val(pp,1)] = ttest2(rxn_time_morn{pp,1}, rxn_time_noon{pp,1});
    % Trial length
    elseif strcmp(Time_Metric, 'Trial_Length')
        [~, box_plot_p_val(pp,1)] = ttest2(trial_length_morn{pp,1}, trial_length_noon{pp,1});
    end
end

%% Find the percent change

avg_box_plot_morn = zeros(length(target_dirs_morn),1);
avg_box_plot_noon = zeros(length(target_dirs_morn),1);
box_plot_perc_change = zeros(length(target_dirs_morn),1);
for pp = 1:length(target_dirs_morn)
    % Reaction Tim
    if strcmp(Time_Metric, 'Rxn_Time')
        avg_box_plot_morn(pp,1) = mean(rxn_time_morn{pp,1});
        avg_box_plot_noon(pp,1) = mean(rxn_time_noon{pp,1});
        box_plot_perc_change(pp,1) = avg_box_plot_noon(pp,1) - avg_box_plot_morn(pp,1) / ...
            abs(avg_box_plot_morn(pp,1));
    % Trial Length
    elseif strcmp(Time_Metric, 'Trial_Length')
        avg_box_plot_morn(pp,1) = mean(trial_length_morn{pp,1});
        avg_box_plot_noon(pp,1) = mean(trial_length_noon{pp,1});
        box_plot_perc_change(pp,1) = avg_box_plot_noon(pp,1) - avg_box_plot_morn(pp,1) / ...
            abs(avg_box_plot_morn(pp,1));
    end
end

%% Find the y-axis limits
max_box_plot = zeros(length(box_plot),1);
for pp = 1:length(max_box_plot)
    max_box_plot(pp) = max(box_plot{pp});
end

y_max = round(max(max_box_plot)) + 0.5;
y_min = -0.3;

%% Plot the box plot data 

for pp = 1:length(target_dirs_morn)

    % Put in the categories for the x-axis
    % Reaction Tim
    if strcmp(Time_Metric, 'Rxn_Time')
        categ_morn = cell(length(rxn_time_morn{pp,1}),1);
        categ_noon = cell(length(rxn_time_noon{pp,1}),1);
    % Trial Length
    elseif strcmp(Time_Metric, 'Trial_Length')
        categ_morn = cell(length(trial_length_morn{pp,1}),1);
        categ_noon = cell(length(trial_length_noon{pp,1}),1);
    end
    for kk = 1:length(categ_morn)
        categ_morn{kk} = 'Morning';
    end
    for kk = 1:length(categ_noon)
        categ_noon{kk} = 'Afternoon';
    end
    categ = [categ_morn; categ_noon];
    
    % Boxplot
    figure
    hold on
    boxplot(box_plot{pp,1}, categ);

    % Increase the axes font
    categ_axes = gca;
    categ_axes.FontSize = label_font_size;

    % Annotation of the n-count
    morn_succ_trials = strcat('n =', {' '}, mat2str(length(categ_morn)));
    noon_succ_trials = strcat('n =', {' '}, mat2str(length(categ_noon)));

    morn_legend = annotation('textbox', [0.35 0.1 0.1 0.1], 'String', ...
        morn_succ_trials, 'FitBoxToText', 'on', 'EdgeColor','none', ...
        'VerticalAlignment', 'top', 'horizontalalignment', 'right');
    morn_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
    noon_legend = annotation('textbox', [0.6 0.1 0.1 0.1], 'String', ...
        noon_succ_trials, 'FitBoxToText', 'on', 'EdgeColor','none', ...
        'VerticalAlignment', 'top', 'horizontalalignment', 'right');
    noon_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;

    % Annotation of the p-value
    if round(box_plot_p_val(pp,1), 3) > 0
        legend_dims = [0.025 0.45 0.44 0.44];
        p_value_string = strcat('p =', {' '}, mat2str(round(box_plot_p_val(pp,1), 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
    end
    if isequal(round(box_plot_p_val(pp,1), 3), 0)
        legend_dims = [0.025 0.45 0.44 0.44];
        p_value_string = strcat('p <', {' '}, '0.001');
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
    end

    % Annotation of the percent change
    if ~isequal(round(box_plot_perc_change(pp,1), 3), 0)
        legend_dims = [0.55 0.45 0.44 0.44];
        perc_change_string = strcat('Δ% =', {' '}, mat2str(round(box_plot_perc_change(pp,1), 3)));
        legend_string = {char(perc_change_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end
    if isequal(round(box_plot_perc_change(pp,1), 3), 0)
        legend_dims = [0.55 0.45 0.44 0.44];
        perc_change_string = strcat('Δ% ≈', {' '}, '0');
        legend_string = {char(perc_change_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

    % Title
    if ~strcmp(Task, 'PG')
        title_string = strcat(Date, {' '}, Task, ',', {' '}, Drug, ':', {' '}, ...
            num2str(target_dirs_noon(pp)), '°,', {' '}, 'TgtCenter at', {' '}, ...
            num2str(target_centers_morn(pp)));
    else
        title_string = strcat(Date, {' '}, Task, ',', {' '}, Drug, ':', {' '}, ...
            'TgtCenter at', {' '}, num2str(target_centers_noon(pp)));
    end
    title(title_string, 'FontSize', title_font_size);
    
    % Labels
    if strcmp(Time_Metric, 'Rxn_Time')
        ylabel('Reaction Time (Sec.)', 'FontSize', label_font_size);
    end
    if strcmp(Time_Metric, 'Trial_Length')
        ylabel('Trial Length (Sec.)', 'FontSize', label_font_size);
    end

    % Axis limits
    ylim([y_min, y_max])
    xlim([0,3]);

end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        save_title = get(fig_info, 'string');
        save_title = strrep(save_title, ':', '');
        save_title = strrep(save_title, 'vs.', 'vs');
        save_title = strrep(save_title, 'mg.', 'mg');
        save_title = strrep(save_title, 'kg.', 'kg');
        save_title = strrep(save_title, '.', '_');
        save_title = strrep(save_title, '/', '_');
        save_title = strrep(save_title, '{ }', ' ');
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








