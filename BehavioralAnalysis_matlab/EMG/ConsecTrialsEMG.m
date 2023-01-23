function ConsecTrialsEMG(xds, EMG_Zero_Factor, EMG_Norm_Factor, trial_num, muscle_groups, Save_Figs)

%% Find the EMG index

[M] = EMG_Index(xds, muscle_groups);

%% Some variable extraction & definitions

% Font specifications
label_font_size = 17;
legend_font_size = 10;
title_font_size = 14;
font_name = 'Arial';

% Do you want to plot the trial lines (Yes = 1; No = 0)
trial_lines = 1;

% Do you want to manually set the y-axis?
man_y_axis = 'No';
%man_y_axis = [-10, 150];

% Times for rewarded trials
[rewarded_start_time] = TrialStartAlignmentTimes(xds, NaN, NaN);
[rewarded_gocue_time] = GoCueAlignmentTimes(xds, NaN, NaN);
[rewarded_end_time] = TrialEndAlignmentTimes(xds, NaN, NaN);

%% If a trial number was selected

if ~strcmp(trial_num, 'All')
    if gt(trial_num, length(rewarded_gocue_time))
        fprintf('trial_num cannot exceed %0.f \n', length(rewarded_gocue_time));
        return
    end
end

%% Find the relative times of the succesful trials (in seconds)
start_to_gocue = zeros(length(rewarded_gocue_time),1);
for ii = 1:length(rewarded_gocue_time)
    start_to_gocue(ii) = rewarded_gocue_time(ii) - rewarded_start_time(ii);
end

start_to_end = zeros(length(rewarded_gocue_time),1);
for ii = 1:length(rewarded_gocue_time)
    start_to_end(ii) = rewarded_end_time(ii) - rewarded_start_time(ii);
end

% Find when the trial started
relative_start_time = zeros(length(rewarded_gocue_time),1);
for ii = 1:length(rewarded_gocue_time)
    if ii == 1
        relative_start_time(ii) = 0;
    else
        relative_start_time(ii) = relative_start_time(ii-1) + start_to_end(ii-1) + xds.bin_width;
    end
end

% Find when the go cue took place in each trial
relative_gocue_time = zeros(length(rewarded_gocue_time),1);
for ii = 1:length(rewarded_gocue_time)
    relative_gocue_time(ii) = relative_start_time(ii) + start_to_gocue(ii);
end

% Find when the end of the trial took place
relative_end_time = zeros(length(rewarded_gocue_time),1);
for ii = 1:length(rewarded_gocue_time)
    relative_end_time(ii) = relative_start_time(ii) + start_to_end(ii);
end

%% Getting the timestamps based on the behavior timings above

trial_timing = struct([]);
for ii = 1:height(rewarded_start_time)
    trial_timing{ii,1} = xds.time_frame((xds.time_frame > rewarded_start_time(ii)) & ... 
        (xds.time_frame < rewarded_end_time(ii)));
end

% Subtracting differences between trial timings to make them consecutive
relative_timing = struct([]);
for ii = 1:height(rewarded_start_time)
    if ii == 1
        relative_timing{ii,1} = trial_timing{ii,1} - rewarded_start_time(ii);
        continue
    end
    relative_timing{ii,1} = trial_timing{ii,1} - rewarded_start_time(ii) + ...
        relative_end_time(ii-1) + xds.bin_width;
end

%% Pull the EMG corresponding to the extracted time frames

% Pull all the EMG
EMG = struct([]);
for ii = 1:height(rewarded_start_time)
    EMG{ii,1} = xds.EMG((xds.time_frame > rewarded_start_time(ii)) & ... 
        (xds.time_frame < rewarded_end_time(ii)),:);
end

% Pull only the selected EMG
selected_EMG = struct([]);
for ii = 1:height(rewarded_start_time)
    for jj = 1:length(M)
        selected_EMG{ii,1}(:,jj) = EMG{ii,1}(:,M(jj));
    end
end

%% Zeroing the EMG

for ii = 1:length(selected_EMG)
    for jj = 1:length(M)
        selected_EMG{ii,1}(:,jj) = selected_EMG{ii,1}(:,jj) - EMG_Zero_Factor(jj);
    end
end

%% Normalizing the average EMG's

norm_EMG = selected_EMG;
for ii = 1:height(norm_EMG)
    for jj = 1:length(M)
        norm_EMG{ii,1}(:,jj) = (norm_EMG{ii,1}(:,jj) /  (EMG_Norm_Factor(jj))*100);
    end
end

%% Concatenating the EMG & relative timing

cat_relative_timing = [];
for ii = 1:height(relative_timing)
    cat_relative_timing = cat(1,cat_relative_timing, relative_timing{ii,1});
end

cat_norm_EMG = [];
for ii = 1:height(norm_EMG)
    cat_norm_EMG = cat(1,cat_norm_EMG, norm_EMG{ii,:});
end

%% Define the number of trials that will be plotted
if ~strcmp(trial_num, 'All')
    line_number = trial_num;
else
    line_number = length(rewarded_start_time);
end

%% Plot the first trials (Top Plot)

% Find the length of the first number of trials
if ~strcmp(trial_num, 'All')
    first_trial_length = relative_timing{trial_num,1}(end,1);
else
    first_trial_length = relative_timing{end,1}(end,1);
end
first_trial_idx = find(cat_relative_timing == first_trial_length);

consec_EMG = figure;
consec_EMG.Position = [300 300 750 450];

hold on
subplot(211)
plot(cat_relative_timing(1:first_trial_idx,1), cat_norm_EMG(1:first_trial_idx,:), ...
    'linewidth', 2)

% Titling the top plot
if ~strcmp(trial_num, 'All')
    title(sprintf('First %i Succesful Trials: EMG', trial_num), 'FontSize', title_font_size)
else
    title('All Succesful Trials: EMG', 'FontSize', title_font_size)
end

% Axis Labels
ylabel('EMG', 'FontSize', label_font_size)
xlabel('Time (Sec.)', 'FontSize', label_font_size)

%% Line indicating go cue and rewards (Top Plot)

% Setting the y-axis limits
y_max = max(max(cat_norm_EMG(1:first_trial_idx,:))) + 30;
y_min = min(min(cat_norm_EMG(1:first_trial_idx,:))) - 1;
if isnan(y_max)
    y_max = 1;
end
if isnan(y_min)
    y_min = 0;
end

if ischar(man_y_axis)
    ylim([y_min, y_max])
else
    ylim(man_y_axis)
end
ylims = ylim;
% Setting the x-axis limits
xlim([0, relative_end_time(line_number)]);

if isequal(trial_lines, 1)
    for ii = 1:line_number
        % Dotted green line indicating beginning of measured window
        line([relative_gocue_time(ii) - 0.4, relative_gocue_time(ii) - 0.4], [ylims(1), ylims(2)], ...
            'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([relative_gocue_time(ii), relative_gocue_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', [0 0.5 0]);
        % Dotted red line indicating beginning of measured window
        line([relative_end_time(ii) - xds.meta.TgtHold, relative_end_time(ii) - xds.meta.TgtHold], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color', 'r','linestyle','--');
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', 'r');
    end
end

% Querry the axes
top_figure_axes = gca;

%% Plot the last number of trials (Bottom Plot)

% Find the length of the last number of trials trials
last_five_start = relative_timing{length(rewarded_start_time) + 1 - line_number,1}(1,1);
last_five_end = relative_timing{length(rewarded_start_time),1}(end,1);
last_five_start_idx = find(cat_relative_timing == last_five_start);
last_five_end_idx = find(cat_relative_timing == last_five_end);

subplot(212)
hold on
plot(cat_relative_timing(last_five_start_idx:last_five_end_idx,1), ... 
    cat_norm_EMG(last_five_start_idx:last_five_end_idx,:), ...
    'linewidth', 2)

% Titling the bottom plot
if ~strcmp(trial_num, 'All')
    title(sprintf('Last %i Succesful Trials: EMG', trial_num), 'FontSize', title_font_size)
else
    title('All Trials: EMG', 'FontSize', title_font_size)
end

% Axes Labels
ylabel('EMG', 'FontSize', label_font_size)
xlabel('Time (Sec.)', 'FontSize', label_font_size)

%% Line indicating go cue and rewards (Bottom Plot)

% Setting the y-axis limits
y_max = max(max(cat_norm_EMG(last_five_start_idx:last_five_end_idx,:))) + 30;
y_min = min(min(cat_norm_EMG(last_five_start_idx:last_five_end_idx,:))) - 1;
if isnan(y_max)
    y_max = 1;
end
if isnan(y_min)
    y_min = 0;
end

if ischar(man_y_axis)
    ylim([y_min, y_max])
else
    ylim(man_y_axis)
end
ylims = ylim;
% Setting the x-axis limits
xlim([relative_start_time(length(relative_end_time) - line_number + 1), relative_end_time(end)]);

if isequal(trial_lines, 1)
    for ii = length(rewarded_start_time) + 1 - line_number:length(rewarded_start_time)
        % Dotted green line indicating beginning of measured window
        line([relative_gocue_time(ii) - 0.4, relative_gocue_time(ii) - 0.4], [ylims(1), ylims(2)], ...
            'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([relative_gocue_time(ii), relative_gocue_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', [0 0.5 0]);
        % Dotted red line indicating beginning of measured window
        line([relative_end_time(ii) - xds.meta.TgtHold, relative_end_time(ii) - xds.meta.TgtHold], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color', 'r','linestyle','--');
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', 'r');
    end
end

% Querry the axes
bottom_figure_axes = gca;

%% Legend containing EMG names
if length(M) == 12
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(4))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(5))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(6))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(7))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(8))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(9))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(10))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(11))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(12))),'EMG_',' ')), ... 
        'NumColumns', 6, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

if length(M) == 6
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(4))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(5))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(6))),'EMG_',' ')), ... 
        'NumColumns', 6, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

if length(M) == 4
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(4))),'EMG_',' ')), ... 
        'NumColumns', 4, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

if length(M) == 3
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        'NumColumns', 3, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

if length(M) == 2
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ... 
        'NumColumns', 2, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

if length(M) == 1
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        'NumColumns', 2, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthEast');
end

legend boxoff

%% Only label every other tick

% Set The Font
set(top_figure_axes,'fontname', font_name);
set(bottom_figure_axes,'fontname', font_name);
% Set ticks to outside
set(top_figure_axes,'TickDir','out');
set(bottom_figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(top_figure_axes,'box','off')
set(bottom_figure_axes,'box','off')

% X Labels
top_x_labels = string(top_figure_axes.XAxis.TickLabels);
bottom_x_labels = string(bottom_figure_axes.XAxis.TickLabels);
top_x_labels(2:2:end) = NaN;
bottom_x_labels(1:1:end) = NaN;
top_figure_axes.XAxis.TickLabels = top_x_labels;
bottom_figure_axes.XAxis.TickLabels = bottom_x_labels;

% Bottom Y Labels
bottom_y_labels = string(bottom_figure_axes.YAxis.TickLabels);
bottom_y_labels(1:1:end) = NaN;
bottom_figure_axes.YAxis.TickLabels = bottom_y_labels;

% Top Y Labels
%top_y_labels = string(top_figure_axes.YAxis.TickLabels);
%top_y_labels(2:2:end) = NaN;
%top_figure_axes.YAxis.TickLabels = top_y_labels;

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:numel(findobj('type','figure'))
        save_title = strcat(num2str(trial_num), {' '}, 'Consecutive Succesful Trials, EMG', {' '}, muscle_groups);
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








    
    
    
    
    
    