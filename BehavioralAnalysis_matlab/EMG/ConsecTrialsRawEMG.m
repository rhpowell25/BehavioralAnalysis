function ConsecTrialsRawEMG(xds, trial_num, muscle_groups, Save_File)

%% Find the EMG index

[M] = EMG_Index(xds, muscle_groups);

%% Some variable extraction & definitions

% Font specifications
label_font_size = 12;
legend_font_size = 12;
title_font_size = 12;
font_name = 'Arial';

% Do you want to plot the trial lines (Yes = 1; No = 0)
trial_lines = 1;

% Extract the raw EMG
xds_raw_EMG = xds.raw_EMG;

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

%% Removing all trials after the raw EMG drops out
raw_EMG_end = xds.raw_EMG_time_frame(end);
drop_out_trials = find(rewarded_start_time >= raw_EMG_end);
if gt(length(drop_out_trials), 1)
    fprintf('The last %0.f trials contain no raw EMG \n', length(drop_out_trials));
    fprintf('The last %0.f trials will be excluded \n', length(drop_out_trials));
    for ii = 1:length(drop_out_trials)
        rewarded_start_time(drop_out_trials(ii)) = NaN;
    end
end

%% Getting the timestamps based on the behavior timings above

trial_timing = struct([]);
for ii = 1:height(rewarded_start_time)
    trial_timing{ii,1} = xds.raw_EMG_time_frame((xds.raw_EMG_time_frame > rewarded_start_time(ii)) & ...
        (xds.raw_EMG_time_frame < rewarded_end_time(ii)));
end

% Subtracting differences between trial timings to make them consecutive
relative_timing = struct([]);
for ii = 1:height(rewarded_start_time)
    if ii == 1
        relative_timing{ii,1} = trial_timing{ii,1} - trial_timing{ii,1}(1,1);
        continue
    end
    relative_timing{ii,1} = trial_timing{ii,1} - trial_timing{ii,1}(1,1) + ...
        relative_timing{ii-1,1}(end,1) + (xds.bin_width  * 2);
end

%% Pull the Raw EMG corresponding to the extracted time frames

% Pull all the raw EMG
raw_EMG = struct([]);
for ii = 1:height(rewarded_start_time)
    raw_EMG{ii,1} = xds_raw_EMG((xds.raw_EMG_time_frame > rewarded_start_time(ii)) & ...
        (xds.raw_EMG_time_frame < rewarded_end_time(ii)),:);
end

% Pull only the selected raw EMG
selected_raw_EMG = struct([]);
for ii = 1:height(rewarded_start_time)
    for jj = 1:length(M)
        selected_raw_EMG{ii,1}(:,jj) = raw_EMG{ii,1}(:,M(jj));
    end
end

%% Concatenating the raw EMG & relative timing

cat_relative_timing = [];
for ii = 1:height(relative_timing)
    cat_relative_timing = cat(1,cat_relative_timing, relative_timing{ii,1});
end

cat_selected_EMG = [];
for ii = 1:height(selected_raw_EMG)
    cat_selected_EMG = cat(1,cat_selected_EMG, selected_raw_EMG{ii,:});
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

consec_raw_EMG = figure;
consec_raw_EMG.Position = [300 300 750 450];

hold on
subplot(211)
plot(cat_relative_timing(1:first_trial_idx,1), cat_selected_EMG(1:first_trial_idx,:), ...
    'linewidth', 1)

% Titling the bottom plot
if ~strcmp(trial_num, 'All')
    Fig_Title = sprintf('First %i Succesful Trials: EMG', trial_num);
else
    Fig_Title = 'All Trials: EMG';
end
title(Fig_Title, 'FontSize', title_font_size)

% Axis Labels
ylabel('Raw EMG',  'FontSize', label_font_size)
xlabel('Time (Sec.)',  'FontSize', label_font_size)

%% Line indicating go cue and rewards (Top Plot)

% Setting the y-axis limits
y_max = max(max(cat_selected_EMG(1:first_trial_idx,:))) + 2;
y_min = min(min(cat_selected_EMG(1:first_trial_idx,:))) - 1;
ylim([y_min, y_max])
ylims = ylim;
% Setting the x-axis limits
xlim([0, relative_end_time(line_number)]);

if isequal(trial_lines, 1) 
    for ii = 1:line_number
        % Dotted green line indicating beginning of measured window
        line([relative_gocue_time(ii) - 0.4, relative_gocue_time(ii) - 0.4], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([relative_gocue_time(ii), relative_gocue_time(ii)], ... 
            [ylims(1), ylims(2)], 'linewidth', 2, 'color', [0 0.5 0]);
        % Dotted red line indicating beginning of measured window
        line([relative_end_time(ii) - xds.meta.TgtHold, relative_end_time(ii) - xds.meta.TgtHold], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color', 'r','linestyle','--');
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], ... 
            [ylims(1), ylims(2)], 'linewidth', 2, 'color', 'r');
    end
end

%% Plot the last number of trials (Bottom Plot)

% Find the length of the last number of trials trials
if ~strcmp(trial_num, 'All')
    last_five_start = relative_timing{length(rewarded_start_time) + 1 - trial_num,1}(1,1);
else
    last_five_start = relative_timing{1,1}(1,1);
end
last_five_end = relative_timing{length(rewarded_start_time),1}(end,1);
last_five_start_idx = find(cat_relative_timing == last_five_start);
last_five_end_idx = find(cat_relative_timing == last_five_end);

subplot(212)
hold on
plot(cat_relative_timing(last_five_start_idx:last_five_end_idx,1), ... 
    cat_selected_EMG(last_five_start_idx:last_five_end_idx,:), ...
    'linewidth', 1)

% Titling the bottom plot
if ~strcmp(trial_num, 'All')
    Fig_Title = sprintf('Last %i Succesful Trials: EMG', trial_num);
else
    Fig_Title = 'All Trials: EMG';
end
title(Fig_Title, 'FontSize', title_font_size)

% Axes Labels
xlabel('Time (Sec.)',  'FontSize', label_font_size)
ylabel('Raw EMG',  'FontSize', label_font_size)

%% Line indicating go cue and rewards (Bottom Plot)

% Setting the y-axis limits
y_max = max(max(cat_selected_EMG(last_five_start_idx:last_five_end_idx,:))) + 2;
y_min = min(min(cat_selected_EMG(last_five_start_idx:last_five_end_idx,:))) - 1;
ylim([y_min, y_max])
ylims = ylim;
% Setting the x-axis limits
xlim([relative_start_time(length(relative_end_time) - line_number + 1), relative_end_time(end)]);

if isequal(trial_lines, 1)
    for ii = length(rewarded_start_time) + 1 - line_number:length(rewarded_start_time)
        % Dotted green line indicating beginning of measured window
        line([relative_gocue_time(ii) - 0.4, relative_gocue_time(ii) - 0.4], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([relative_gocue_time(ii), relative_gocue_time(ii)], ... 
            [ylims(1), ylims(2)], 'linewidth', 2, 'color', [0 0.5 0]);
        % Dotted red line indicating beginning of measured window
        line([relative_end_time(ii) - xds.meta.TgtHold, relative_end_time(ii) - xds.meta.TgtHold], ... 
            [ylims(1), ylims(2)], 'linewidth',2,'color', 'r','linestyle','--');
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], ... 
            [ylims(1), ylims(2)], 'linewidth', 2, 'color', 'r');
    end
end

%% Legend containing EMG names
if isequal(length(M), 12)
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
        'Location', 'NorthWest');
end

if isequal(length(M), 6)
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(4))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(5))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(6))),'EMG_',' ')), ... 
        'NumColumns', 6, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthWest');
end

if isequal(length(M), 4)
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(3))),'EMG_',' ')), ...
        sprintf('%s', strrep(string(xds.EMG_names(M(4))),'EMG_',' ')), ... 
        'NumColumns', 4, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthWest');
end

if isequal(length(M), 2)
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        sprintf('%s', strrep(string(xds.EMG_names(M(2))),'EMG_',' ')), ... 
        'NumColumns', 4, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthWest');
end

if isequal(length(M), 1)
    legend(sprintf('%s', strrep(string(xds.EMG_names(M(1))),'EMG_',' ')), ... 
        'NumColumns', 4, 'FontSize', legend_font_size, 'FontName', font_name, ...
        'Location', 'NorthWest');
end

legend boxoff

%% Save the file if selected
Save_Figs(Fig_Title, Save_File)








    
    
    
    
    
    