function Per_Trial_CursorSignal(xds, signal_choice, event, unit_name, cursor_YLims, Save_Figs)

%% File Description:

% This function plots the trial-by-trial cursor position, velocity, or 
% acceleration per target direction / distance of all the succesful trials 
% in an xds file.
% If you set Save_Figs to 0, the figure will not be saved to your desktop.
%
% -- Inputs --
% xds: the xds file
% signal_choice: 'Pos', 'Vel', or 'Acc'
% event: 'trial_gocue', 'window_trial_gocue', 'trial_end', 
% 'window_trial_end', 'force_onset', 'window_force_onset', 'force_max', 
% 'window_force_max', 'window_force_deriv', 'force_deriv', 'cursor_onset', 
% 'window_cursor_onset', 'cursor_veloc', 'window_cursor_veloc', 
% 'cursor_acc', 'window_cursor_acc', 'EMG_max', 'window_EMG_max', 
% 'task_onset', or 'window_task_onset'
% unit_name: 'elec1_1', #, or NaN
% cursor_YLims: [ymax, ymin]
% Save_Figs: 'pdf', 'png', 'fig', or 0

%% End the function if there is no Y-Limit

if isnan(cursor_YLims)
    disp("There is no Y-Limit")
    return
end

%% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

%% Basic Settings, some variable extractions, & definitions

% Pull the binning paramaters
[Bin_Params] = Binning_Parameters;

% Time before & after the event
before_event = Bin_Params.before_event;
after_event = Bin_Params.after_event;

% Window to calculate max firing rate
half_window_length = Bin_Params.half_window_length; % Time (sec.)

if ~contains(event, 'window')
    max_fr_time = 0;
end

if contains(event, 'gocue') || contains(event, 'force_onset')
    % Define the window for the baseline phase
    time_before_gocue = 0.4;
elseif contains(event, 'end')
    % Define the window for the movement phase
    time_before_end = xds.meta.TgtHold;
end

% Extract the cursor signal of chhoice
if strcmp(signal_choice, 'Pos')
    curs_sig = xds.curs_p;
elseif strcmp(signal_choice, 'Vel')
    curs_sig = xds.curs_v;
elseif strcmp(signal_choice, 'Acc')
    curs_sig = xds.curs_a;
end

% Font specifications
label_font_size = 12;
title_font_size = 13;
axes_line_size = 1;
plot_line_size = 3;
font_name = 'Arial';
figure_width = 600;
figure_height = 600;

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dirs = length(target_dirs);

%% Begin the loop through all directions
for jj = 1:num_dirs
    
    %% Times for rewarded trials
    if strcmp(event, 'trial_gocue')
        [rewarded_gocue_time] = GoCueAlignmentTimes(xds, NaN, NaN);
        [rewarded_end_time] = TrialEndAlignmentTimes(xds, NaN, NaN);
        [Alignment_Times] = EventAlignmentTimes(xds, NaN, NaN, event);
    else
        [rewarded_gocue_time] = GoCueAlignmentTimes(xds, target_dirs(jj), target_centers(jj));
        [rewarded_end_time] = TrialEndAlignmentTimes(xds, target_dirs(jj), target_centers(jj));
        [Alignment_Times] = EventAlignmentTimes(xds, target_dirs(jj), target_centers(jj), event);
    end

    if contains(event, 'window')
        % Run the preferred direction window function
        [~, max_fr_time] = ...
        EventWindow(xds, unit_name, target_dirs(jj), target_centers(jj), event);
    end
    
    %% Times between events
    % Find time between the go-cue and reward
    gocue_to_event = Alignment_Times - rewarded_gocue_time;
    event_to_end = rewarded_end_time - Alignment_Times;

    %% Extracting cursor signal & time during successful trials

    % Cursor signal & time measured during each successful trial 
    rewarded_curs_sig = struct([]);
    % Time points during each succesful trial
    timings = struct([]); 
    for ii = 1:length(Alignment_Times)
        temp_start = 0;
        alignment_idx = find(xds.time_frame == Alignment_Times(ii));
        alignment_start_idx = alignment_idx - (before_event / xds.bin_width);
        if alignment_start_idx < 0
            temp_start = alignment_start_idx;
            alignment_start_idx = 1;
        end
        alignment_end_idx = alignment_idx + (after_event / xds.bin_width);
        rewarded_curs_sig{ii,1} = curs_sig(alignment_start_idx : alignment_end_idx, :);
        timings{ii, 1} = xds.time_frame(alignment_start_idx : alignment_end_idx);
        if temp_start ~= 0
            rewarded_curs_sig{ii,1} = cat(1, NaN(abs(temp_start) + 1, width(rewarded_curs_sig{ii})), rewarded_curs_sig{ii,1});
            timings{ii,1} = cat(1, NaN(abs(temp_start) + 1, 1), timings{ii,1});
        end
    end

    %% Find the vector sum of the cursor signal

    z_cur_sig = struct([]);
    for ii = 1:length(rewarded_gocue_time)
        z_cur_sig{ii,1} = zeros(length(rewarded_curs_sig{ii,1}(:,1)),1);
        for dd = 1:length(z_cur_sig{ii,1})
            z_cur_sig{ii,1}(dd,1) = sqrt(rewarded_curs_sig{ii,1}(dd,1).^2 + rewarded_curs_sig{ii,1}(dd,2).^2);
        end
    end

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(rewarded_curs_sig{1,1}));

    %% Plot the decomposed individual cursor signals on the top

    Cursor_figure = figure;
    Cursor_figure.Position = [300 150 figure_width figure_height];
    tiledlayout(2,1);
    nexttile
    hold on
    for ii = 1:length(rewarded_gocue_time)
        plot(absolute_timing, rewarded_curs_sig{ii,1}(:,1), 'b', 'LineWidth',.2, 'LineStyle','--')
        plot(absolute_timing, rewarded_curs_sig{ii,1}(:,2), 'g', 'LineWidth',.2, 'LineStyle','--')
    end
    
    % Setting the axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([cursor_YLims(2), cursor_YLims(1)]);
    ylims = ylim;
    
    if contains(event, 'gocue')
        % Solid dark green line indicating the aligned time
        line([0, 0], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'Color', [0 0.5 0]);
        % Dotted dark green line indicating beginning of measured window
        line([-time_before_gocue, -time_before_gocue], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'Color', [0 0.5 0], 'LineStyle','--');
    elseif contains(event, 'end')
        % Solid red line indicating the aligned time
        line([0, 0], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'color', 'r');
        % Dotted red line indicating beginning of measured window
        line([-time_before_end, -time_before_end], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'color','r','linestyle','--');
    end

    if contains(event, 'window')
        % Dotted purple line indicating beginning of measured window
        line([max_fr_time - half_window_length, max_fr_time - half_window_length], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
        % Dotted purple line indicating end of measured window
        line([max_fr_time + half_window_length, max_fr_time + half_window_length], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
    elseif ~contains(event, 'trial_gocue') && ~contains(event, 'trial_end')
        % Dotted red line indicating beginning of measured window
        line([-0.1, -0.1], [ylims(1), ylims(2)], ...
            'Linewidth', plot_line_size, 'Color', 'r', 'Linestyle','--');
        % Dotted red line indicating end of measured window
        line([0.1, 0.1], [ylims(1), ylims(2)], ...
            'Linewidth', plot_line_size, 'Color', 'r', 'Linestyle','--');
    end
    
    % Define the labels
    if strcmp(signal_choice, 'Pos')
        signal_label = 'Wrist Position';
    elseif strcmp(signal_choice, 'Vel')
        signal_label = 'Wrist Velocity';
    elseif strcmp(signal_choice, 'Acc')
        signal_label = 'Wrist Acceleration';
    end

    % Labeling the axis
    ylabel(signal_label, 'FontSize', label_font_size);
    
    % Titling the top plot
    title(sprintf('%s: %i°, TgtCenter at %0.1f', ... 
        signal_label, target_dirs(jj), target_centers(jj)), 'FontSize', title_font_size)

    % Remove x-axis ticks
    figure_axes = gca;
    figure_axes.LineWidth = axes_line_size;
    x_labels = string(figure_axes.XAxis.TickLabels);
    x_labels(1:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    % Set ticks to outside
    set(gca,'TickDir','out');
    % Remove the top and right tick marks
    set(gca,'box','off')
    % Set The Font
    set(figure_axes,'FontName', font_name);

    %% Plot the individual cursor signals on the bottom

    nexttile
    hold on
    for ii = 1:length(rewarded_gocue_time)
        plot(absolute_timing, z_cur_sig{ii,1}, 'LineWidth',.2)
    end
    
    for ii = 1:length(rewarded_gocue_time)
        cursor_gocue_idx = timings{ii,1} == rewarded_gocue_time(ii);
        cursor_end_idx = timings{ii,1} == rewarded_end_time(ii);
        % Plot the go-cues as dark green dots
        plot(-gocue_to_event(ii), z_cur_sig{ii,1}(cursor_gocue_idx), ...
            'Marker', '.', 'Color', [0 0.5 0], 'Markersize', 15);
        % Plot the trial ends as red dots
        plot(event_to_end(ii), z_cur_sig{ii,1}(cursor_end_idx), ...
            'Marker', '.', 'Color', 'r', 'Markersize', 15);
    end
    
    % Setting the axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([-1, cursor_YLims(1)]);
    ylims = ylim;
    
    if contains(event, 'gocue')
        % Solid dark green line indicating the aligned time
        line([0, 0], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'Color', [0 0.5 0]);
        % Dotted dark green line indicating beginning of measured window
        line([-time_before_gocue, -time_before_gocue], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'Color', [0 0.5 0], 'LineStyle','--');
    elseif contains(event, 'end')
        % Solid red line indicating the aligned time
        line([0, 0], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'color', 'r');
        % Dotted red line indicating beginning of measured window
        line([-time_before_end, -time_before_end], [ylims(1), ylims(2)], ...
            'LineWidth', plot_line_size, 'color','r','linestyle','--');
    end

    if contains(event, 'window')
        % Dotted purple line indicating beginning of measured window
        line([max_fr_time - half_window_length, max_fr_time - half_window_length], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
        % Dotted purple line indicating end of measured window
        line([max_fr_time + half_window_length, max_fr_time + half_window_length], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
    elseif ~contains(event, 'trial_gocue') && ~contains(event, 'trial_end')
        % Dotted red line indicating beginning of measured window
        line([-0.1, -0.1], [ylims(1), ylims(2)], ...
            'Linewidth', plot_line_size, 'Color', 'r', 'Linestyle','--');
        % Dotted red line indicating end of measured window
        line([0.1, 0.1], [ylims(1), ylims(2)], ...
            'Linewidth', plot_line_size, 'Color', 'r', 'Linestyle','--');
    end
    
    % Labeling the axis
    ylabel(signal_label, 'FontSize', label_font_size);
    xlabel('Time (sec.)', 'FontSize', label_font_size);

    % Only label every other tick
    figure_axes = gca;
    figure_axes.LineWidth = axes_line_size;
    x_labels = string(figure_axes.XAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    % Set ticks to outside
    set(gca,'TickDir','out');
    % Remove the top and right tick marks
    set(gca,'box','off')
    % Set The Font
    set(figure_axes,'FontName', font_name);

    % End the event after one loop if showing baseline firing rate
    if strcmp(event, 'trial_gocue')
        return
    end
    
end % End of target loop

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


