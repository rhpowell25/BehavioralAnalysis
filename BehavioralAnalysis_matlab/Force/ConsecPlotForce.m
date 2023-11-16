function ConsecPlotForce(xds, Force_Norm_Factor, trial_num, Save_File)

%% Ending the function if there is no force

if strcmp(xds.meta.task, 'WS')
    disp('Event cannot be force related for this task');
    return
end

if xds.has_force == 0
    disp('No force in this file')
    return
end

%% Some variable extraction & definitions

% Font specifications
label_font_size = 17;
title_font_size = 14;
font_name = 'Arial';

% Do you want to plot the trial lines (Yes = 1; No = 0)
trial_lines = 0;

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
for ii = 1:height(rewarded_gocue_time)
    trial_timing{ii,1} = xds.time_frame((xds.time_frame > rewarded_start_time(ii)) & ... 
        (xds.time_frame < rewarded_end_time(ii)));
end

% Subtracting differences between trial timings to make them consecutive
relative_timing = struct([]);
for ii = 1:height(rewarded_gocue_time)
    if ii == 1
        relative_timing{ii,1} = trial_timing{ii,1} - rewarded_start_time(ii);
        continue
    end
    relative_timing{ii,1} = trial_timing{ii,1} - rewarded_start_time(ii) + ...
        relative_end_time(ii-1) + xds.bin_width;
end

%% Pull the force corresponding to the extracted time frames

% Pull all the Force
Force = struct([]);
for ii = 1:height(rewarded_gocue_time)
    Force{ii,1} = xds.force((xds.time_frame > rewarded_start_time(ii)) & ... 
        (xds.time_frame < rewarded_end_time(ii)),:);
end

%% Normalizing the force
if ~ischar(Force_Norm_Factor)
    for ii = 1:length(rewarded_gocue_time)
        Force{ii,1} = Force{ii,1} / Force_Norm_Factor*100;
    end
else
    for ii = 1:length(rewarded_gocue_time)
        Force{ii,1} = Force{ii,1} / 1000*5; % Millivolt conversion * gain
    end
end

%% Sum the two force transducers
[Sigma_Force] = Sum_Force(xds.meta.task, Force);

%% Concatenating the force & relative timing

cat_relative_timing = [];
for ii = 1:height(relative_timing)
    cat_relative_timing = cat(1,cat_relative_timing, relative_timing{ii,1});
end

cat_force = [];
for ii = 1:height(Sigma_Force)
    cat_force = cat(1,cat_force, Sigma_Force{ii,:});
end

%% Define the number of trials that will be plotted
if ~strcmp(trial_num, 'All')
    line_number = trial_num;
else
    line_number = length(rewarded_gocue_time);
end

%% Plot the first trials (Top Plot)

% Find the length of the first number of trials
if ~strcmp(trial_num, 'All')
    first_trial_length = relative_timing{trial_num,1}(end,1);
else
    first_trial_length = relative_timing{end,1}(end,1);
end
first_trial_idx = find(cat_relative_timing == first_trial_length);

consec_force = figure;
consec_force.Position = [300 300 750 450];

hold on
subplot(211)
plot(cat_relative_timing(1:first_trial_idx,1), cat_force(1:first_trial_idx,:), ...
    'linewidth', 1.5, 'Color', 'k')

% Titling the bottom plot
if ~strcmp(trial_num, 'All')
    Fig_Title = sprintf('First %i Succesful Trials: Force', trial_num);
else
    Fig_Title = 'All Succesful Trials: Force';
end
title(Fig_Title, 'FontSize', title_font_size)

% Axis Labels
ylabel('Force', 'FontSize', label_font_size)
xlabel('Time (Sec.)', 'FontSize', label_font_size)

%% Line indicating go cue and rewards (Top Plot)

% Setting the y-axis limits
y_max = max(max(cat_force(1:first_trial_idx,:))) + 2;
y_min = min(min(cat_force(1:first_trial_idx,:))) - 1;
ylim([y_min, y_max])
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
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', 'r');
    end
end

%% Only label every other tick
figure_axes = gca;
% Set The Font
set(figure_axes,'fontname', font_name);
x_labels = string(figure_axes.XAxis.TickLabels);
y_labels = string(figure_axes.YAxis.TickLabels);
x_labels(2:2:end) = NaN;
y_labels(2:2:end) = NaN;
figure_axes.XAxis.TickLabels = x_labels;
figure_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')

%% Plot the last number of trials (Bottom Plot)

% Find the length of the last number of trials trials
last_five_start = relative_timing{length(rewarded_gocue_time) + 1 - line_number,1}(1,1);
last_five_end = relative_timing{length(rewarded_gocue_time),1}(end,1);
last_five_start_idx = find(cat_relative_timing == last_five_start);
last_five_end_idx = find(cat_relative_timing == last_five_end);

subplot(212)
hold on
plot(cat_relative_timing(last_five_start_idx:last_five_end_idx,1), ... 
    cat_force(last_five_start_idx:last_five_end_idx,:), ...
    'linewidth', 1.5, 'Color' , 'k')

% Titling the bottom plot
if ~strcmp(trial_num, 'All')
    Fig_Title = sprintf('Last %i Succesful Trials: Force', trial_num);
else
    Fig_Title = 'All Succesful Trials: Force';
end
title(Fig_Title, 'FontSize', title_font_size)

% Axes Labels
ylabel('Force', 'FontSize', label_font_size)
xlabel('Time (Sec.)', 'FontSize', label_font_size)

%% Line indicating go cue and rewards (Bottom Plot)

% Setting the y-axis limits
y_max = max(max(cat_force(last_five_start_idx:last_five_end_idx,:))) + 2;
y_min = min(min(cat_force(last_five_start_idx:last_five_end_idx,:))) - 1;
ylim([y_min, y_max])
ylims = ylim;
% Setting the x-axis limits
xlim([relative_start_time(length(relative_end_time) - line_number + 1), relative_end_time(end)]);

if isequal(trial_lines, 1)
    for ii = length(rewarded_gocue_time) + 1 - line_number:length(rewarded_gocue_time)
        % Dotted green line indicating beginning of measured window
        line([relative_gocue_time(ii) - 0.4, relative_gocue_time(ii) - 0.4], [ylims(1), ylims(2)], ...
            'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([relative_gocue_time(ii), relative_gocue_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', [0 0.5 0]);
        % Solid red line indicating the aligned time
        line([relative_end_time(ii), relative_end_time(ii)], [ylims(1), ylims(2)], ...
            'linewidth', 2, 'color', 'r');
    end
end

%% Only label every other tick
figure_axes = gca;
% Set The Font
set(figure_axes,'fontname', font_name);
x_labels = string(figure_axes.XAxis.TickLabels);
y_labels = string(figure_axes.YAxis.TickLabels);
x_labels(2:2:end) = NaN;
y_labels(2:2:end) = NaN;
figure_axes.XAxis.TickLabels = x_labels;
figure_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off')

%% Save the file if selected
Save_Figs(Fig_Title, Save_File)



    
    
    
    
    
    