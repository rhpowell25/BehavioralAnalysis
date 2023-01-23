function Per_Trial_CursorPos(xds, event, unit_name, cursor_pos_YLims, Save_Figs)

%% End the function if there is no Y-Limit

if isnan(cursor_pos_YLims)
    disp("There is no Y-Limit")
    return
end

%% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

%% Basic Settings, some variable extractions, & definitions

% Event lengths
before_event = 3;
after_event = 3;

% Window to calculate max firing rate
window_size = 0.1;

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

% Font specifications
label_font_size = 12;
title_font_size = 13;
plot_line_size = 3;

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
        [~, max_fr_time, ~] = ...
        EventWindow(xds, unit_name, target_dirs(jj), target_centers(jj), event);
    end
    
    %% Times between events
    % Find time between the go-cue and reward
    gocue_to_event = Alignment_Times - rewarded_gocue_time;
    event_to_end = rewarded_end_time - Alignment_Times;

    %% Extracting cursor position & time during successful trials

    % Cursor position & time measured during each successful trial 
    cursor_p = struct([]); % Cursor position during each successful trial
    timings = struct([]); % Time points during each succesful trial 
    for ii = 1:length(Alignment_Times)
        temp_start = 0;
        alignment_idx = find(xds.time_frame == Alignment_Times(ii));
        alignment_start_idx = alignment_idx - (before_event / xds.bin_width);
        if alignment_start_idx < 0
            temp_start = alignment_start_idx;
            alignment_start_idx = 1;
        end
        alignment_end_idx = alignment_idx + (after_event / xds.bin_width);
        cursor_p{ii,1} = xds.curs_p(alignment_start_idx : alignment_end_idx, :);
        timings{ii, 1} = xds.time_frame(alignment_start_idx : alignment_end_idx);
        if temp_start ~= 0
            cursor_p{ii,1} = cat(1, NaN(abs(temp_start) + 1, width(cursor_p{ii})), cursor_p{ii,1});
            timings{ii,1} = cat(1, NaN(abs(temp_start) + 1, 1), timings{ii,1});
        end
    end

    %% Find the vector sum of the cursor position

    z_cursor_p = struct([]);
    for ii = 1:length(rewarded_gocue_time)
        z_cursor_p{ii,1} = zeros(length(cursor_p{ii,1}(:,1)),1);
        for dd = 1:length(z_cursor_p{ii,1})
            z_cursor_p{ii,1}(dd,1) = sqrt(cursor_p{ii,1}(dd,1).^2 + cursor_p{ii,1}(dd,2).^2);
        end
    end

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(cursor_p{1,1}));

    %% Plot the decomposed individual cursor positions on the top

    figure
    subplot(211);
    hold on
    for ii = 1:length(rewarded_gocue_time)

        plot(absolute_timing, cursor_p{ii,1}(:,1), 'b', 'LineWidth',.2, 'LineStyle','--')
        plot(absolute_timing, cursor_p{ii,1}(:,2), 'g', 'LineWidth',.2, 'LineStyle','--')
    end
    
    % Setting the x-axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([cursor_pos_YLims(2), cursor_pos_YLims(1)]);
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
        line([max_fr_time - window_size, max_fr_time - window_size], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
        % Dotted purple line indicating end of measured window
        line([max_fr_time + window_size, max_fr_time + window_size], ... 
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
    ylabel('Wrist Position', 'FontSize', label_font_size);
    xlabel('Time (sec.)', 'FontSize', label_font_size);
    
    % Titling the top plot
    title(sprintf('Decomposed wrist position: %i°, TgtCenter at %0.1f', ... 
        target_dirs(jj), target_centers(jj)), 'FontSize', title_font_size)

    %% Plot the individual cursor positions on the bottom

    subplot(212);
    hold on
    for ii = 1:length(rewarded_gocue_time)
        plot(absolute_timing, z_cursor_p{ii,1}, 'k', 'LineWidth',.2, 'LineStyle','--')
    end
    
    for ii = 1:length(rewarded_gocue_time)
        cursor_gocue_idx = timings{ii,1} == rewarded_gocue_time(ii);
        cursor_end_idx = timings{ii,1} == rewarded_end_time(ii);
        % Plot the go-cues as dark green dots
        plot(-gocue_to_event(ii), z_cursor_p{ii,1}(cursor_gocue_idx), ...
            'Marker', '.', 'Color', [0 0.5 0], 'Markersize', 15);
        % Plot the trial ends as red dots
        plot(event_to_end(ii), z_cursor_p{ii,1}(cursor_end_idx), ...
            'Marker', '.', 'Color', 'r', 'Markersize', 15);

    end
    
    % Setting the x-axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([-1, cursor_pos_YLims(1)]);
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
        line([max_fr_time - window_size, max_fr_time - window_size], ... 
            [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
        % Dotted purple line indicating end of measured window
        line([max_fr_time + window_size, max_fr_time + window_size], ... 
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
    ylabel('Wrist Position', 'FontSize', label_font_size);
    xlabel('Time (sec.)', 'FontSize', label_font_size);
    
    % Titling the top plot
    title(sprintf('Wrist position: %i°, TgtCenter at %0.1f', ... 
        target_dirs(jj), target_centers(jj)), 'FontSize', title_font_size)

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


