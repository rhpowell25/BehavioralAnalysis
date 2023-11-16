function Per_Trial_Force(xds, event, unit_name, Force_Norm_Factor, force_YLims, Save_File)

%% End the function if there is no Y-Limit

if isnan(force_YLims)
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
    try
        time_before_end = xds.meta.TgtHold;
    catch
        time_before_end = NaN;
    end
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
        [~, max_fr_time] = ...
        EventWindow(xds, unit_name, target_dirs(jj), target_centers(jj), event);
    end
    
    %% Times between events
    % Find time between the go-cue and reward
    gocue_to_event = Alignment_Times - rewarded_gocue_time;
    event_to_end = rewarded_end_time - Alignment_Times;

    %% Extracting force & time during successful trials
    % Force measured during each successful trial
    Force = struct([]); % Force during each successful trial
    timings = struct([]); % Time points during each succesful trial 
    for ii = 1:length(Alignment_Times)
        Force_idx = find(xds.time_frame == Alignment_Times(ii));
        try
            Force{ii,1} = xds.force((Force_idx - (before_event / xds.bin_width)) : (Force_idx + (after_event / xds.bin_width)), :);
            timings{ii, 1} = xds.time_frame(Force_idx - (before_event / xds.bin_width) : ...
                Force_idx + (after_event / xds.bin_width));
        catch
            Force{ii, 1} = xds.force((Force_idx - (before_event / xds.bin_width)) : end, :);
            Force{ii, 1} = cat(1, Force{ii, 1}, NaN(length(Force{ii-1, 1}) - length(Force{ii, 1}), 2));
            timings{ii, 1} = xds.time_frame(Force_idx - (before_event / xds.bin_width) : end);
            timings{ii, 1} = cat(1, timings{ii, 1}, NaN(length(timings{ii-1, 1}) - length(timings{ii, 1}), 1));
        end
    end

    %% Normalizing the force
    if ~ischar(Force_Norm_Factor)
        for ii = 1:length(Alignment_Times)
            Force{ii,1} = Force{ii,1} / Force_Norm_Factor*100;
        end
    else
        for ii = 1:length(Alignment_Times)
            Force{ii,1} = Force{ii,1} / 1000*5; % Millivolt conversion * gain
        end
    end

    %% Sum the two force transducers
    [Sigma_Force] = Sum_Force(xds.meta.task, Force);

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(Sigma_Force{1,1}));

    %% Plot the decomposed individual force positions on the top

    figure
    subplot(211);
    hold on
    for ii = 1:length(rewarded_gocue_time)
        plot(absolute_timing, Force{ii,1}(:,1), 'b', 'LineWidth',.2, 'LineStyle','--')
        plot(absolute_timing, Force{ii,1}(:,2), 'g', 'LineWidth',.2, 'LineStyle','--')
    end
    
    % Setting the x-axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([force_YLims(2), force_YLims(1)]);
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
    ylabel('Force', 'FontSize', label_font_size);
    xlabel('Time (sec.)', 'FontSize', label_font_size);
    
    % Titling the top plot
    Fig_Title = sprintf('Decomposed Force: %i°, TgtCenter at %0.1f', ... 
        target_dirs(jj), target_centers(jj));
    title(Fig_Title, 'FontSize', title_font_size)

    %% Plot the individual force positions on the bottom

    subplot(212);
    hold on
    for ii = 1:length(rewarded_gocue_time)
        plot(absolute_timing, Sigma_Force{ii,1}, 'LineWidth',.2)
    end
    
    for ii = 1:length(rewarded_gocue_time)
        force_gocue_idx = timings{ii,1} == rewarded_gocue_time(ii);
        force_end_idx = timings{ii,1} == rewarded_end_time(ii);
        % Plot the go-cues as dark green dots
        if ~isempty(Sigma_Force{ii,1}(force_gocue_idx))
            plot(-gocue_to_event(ii), Sigma_Force{ii,1}(force_gocue_idx), ...
                'Marker', '.', 'Color', [0 0.5 0], 'Markersize', 15);
        end
        % Plot the trial ends as red dots
        if ~isempty(Sigma_Force{ii,1}(force_end_idx))
            plot(event_to_end(ii), Sigma_Force{ii,1}(force_end_idx), ...
                'Marker', '.', 'Color', 'r', 'Markersize', 15);
        end
    end
    
    % Setting the x-axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([-1, force_YLims(1)]);
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
    ylabel('Force', 'FontSize', label_font_size);
    xlabel('Time (sec.)', 'FontSize', label_font_size);
    
    % Titling the bottom plot
    title(sprintf('Force: %i°, TgtCenter at %0.1f', ... 
        target_dirs(jj), target_centers(jj)), 'FontSize', title_font_size)

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

    %% End the event after one loop if showing baseline firing rate
    if strcmp(event, 'trial_gocue')
        return
    end
    
end % End of target loop


