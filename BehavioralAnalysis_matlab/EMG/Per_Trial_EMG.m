function Per_Trial_EMG(xds, event, EMG_Zero_Factor, EMG_Norm_Factor, muscle_group, Save_File)

%% Display the function being used
disp('Per Trial EMG Function:');

%% Find the EMG index
[M] = EMG_Index(xds, muscle_group);

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

% Binning information
bin_size = xds.bin_width; % Time (sec.)

if ~contains(event, 'window')
    max_amp_time = 0;
end

if contains(event, 'gocue') || contains(event, 'force_onset')
    % Define the window for the baseline phase
    time_before_gocue = 0.4;
elseif contains(event, 'end')
    % Define the window for the movement phase
    time_before_end = xds.meta.TgtHold;
end

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

%% X-axis
EMG_time = (-before_event:bin_size:after_event);
EMG_time = EMG_time(1:end-1) + bin_size/2;

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
        % Get the average EMG & peak amplitude time
        [~, max_amp_time] = EventWindowEMG(xds, muscle_group, target_dirs(jj), target_centers(jj), event);
    end

    %% Times between events
    % Find time between the go-cue and reward
    gocue_to_event = Alignment_Times - rewarded_gocue_time;
    event_to_end = rewarded_end_time - Alignment_Times;

    %% Extracting EMG & time during successful trials

    aligned_EMG = struct([]); % EMG during each successful trial
    timings = struct([]); % Time points during each succesful trial
    for ii = 1:length(rewarded_gocue_time)
        EMG_idx = find(xds.time_frame == Alignment_Times(ii));
        try
            aligned_EMG{ii, 1} = xds.EMG((EMG_idx - (before_event / xds.bin_width) + 1) : (EMG_idx + (after_event / xds.bin_width)), M);
            timings{ii, 1} = xds.time_frame(EMG_idx - (before_event / xds.bin_width) + 1: ...
                EMG_idx + (after_event / xds.bin_width));
        catch
            aligned_EMG{ii, 1} = xds.EMG((EMG_idx - (before_event / xds.bin_width) + 1) : end, M);
            aligned_EMG{ii, 1} = cat(1, aligned_EMG{ii, 1}, NaN(length(aligned_EMG{ii-1, 1}) - length(aligned_EMG{ii, 1}), length(M)));
            timings{ii, 1} = xds.time_frame(EMG_idx - (before_event / xds.bin_width) + 1 : end);
            timings{ii, 1} = cat(1, timings{ii, 1}, NaN(length(timings{ii-1, 1}) - length(timings{ii, 1}), 1));
        end
    end

    %% Putting all succesful trials in one array
    all_trials_EMG = struct([]);
    for ii = 1:length(M)
        all_trials_EMG{ii,1} = zeros(length(aligned_EMG{1,1}), length(rewarded_gocue_time));
        for mm = 1:length(rewarded_gocue_time)
            all_trials_EMG{ii,1}(:,mm) = aligned_EMG{mm, 1}(:, ii);
        end
    end

    %% Zeroing the EMG
    for ii = 1:length(M)
        all_trials_EMG{ii,1} = all_trials_EMG{ii,1} - EMG_Zero_Factor(ii);
    end

    %% Normalizing the average EMG's
    for ii = 1:length(M)
        all_trials_EMG{ii,1} = (all_trials_EMG{ii,1} / EMG_Norm_Factor(ii))*100;
    end

    %% Plot the individual EMG traces on the top

    for ii = 1:length(M)

        % Find the max & min EMG for the y-axis limits
        y_max = round(max(all_trials_EMG{ii,1}, [], 'All')) + 10;
        y_min = round(min(all_trials_EMG{ii,1}, [], 'All')) - 10;

        EMG_figure = figure;
        EMG_figure.Position = [300 300 Plot_Params.fig_size Plot_Params.fig_size / 2];
        hold on

        % Titling the plot
        EMG_name = strrep(string(xds.EMG_names(M(ii))),'EMG_','');
        Fig_Title = sprintf('EMG, %i°, TgtCenter at %0.1f: %s', ...
            target_dirs(jj), target_centers(jj), EMG_name);
        if contains(xds.meta.rawFileName, 'Pre')
            Fig_Title = strcat(Fig_Title, ' (Morning)');
        end
        if contains(xds.meta.rawFileName, 'Post')
            Fig_Title = strcat(Fig_Title, ' (Afternoon)');
        end
        title(Fig_Title, 'FontSize', Plot_Params.title_font_size)

        % Labels
        ylabel('EMG', 'FontSize', Plot_Params.label_font_size);
        xlabel('Time (sec.)', 'FontSize', Plot_Params.label_font_size);

        % Setting the x-axis limits
        if contains(event, 'gocue')
            xlim([-before_event + 2, after_event]);
        elseif contains(event, 'end')
            xlim([-before_event, after_event - 2]);
        else
            xlim([-before_event + 1, after_event - 1]);
        end
        ylim([y_min, y_max]);
        ylims = ylim;

        for pp = 1:width(all_trials_EMG{ii})

            plot(EMG_time, all_trials_EMG{ii}(:,pp))

            EMG_gocue_idx = timings{pp,1} == rewarded_gocue_time(pp);
            EMG_end_idx = timings{pp,1} == rewarded_end_time(pp);
            % Plot the go-cues as dark green dots
            if ~isempty(all_trials_EMG{ii,1}(EMG_gocue_idx,pp))
                plot(-gocue_to_event(pp), all_trials_EMG{ii,1}(EMG_gocue_idx,pp), ...
                    'Marker', '.', 'Color', [0 0.5 0], 'Markersize', 15);
            end
            % Plot the trial ends as red dots
            if ~isempty(all_trials_EMG{ii,1}(EMG_end_idx))
                plot(event_to_end(pp), all_trials_EMG{ii,1}(EMG_end_idx,pp), ...
                    'Marker', '.', 'Color', 'r', 'Markersize', 15);
            end

        end % End of the individual trial loop

        if contains(event, 'gocue')
            % Solid dark green line indicating the aligned time
            line([0, 0], [ylims(1), ylims(2)], ...
                'LineWidth', Plot_Params.mean_line_width, 'Color', [0 0.5 0]);
            % Dotted dark green line indicating beginning of measured window
            line([-time_before_gocue, -time_before_gocue], [ylims(1), ylims(2)], ...
                'LineWidth', Plot_Params.mean_line_width, 'Color', [0 0.5 0], 'LineStyle','--');
        elseif contains(event, 'end')
            % Solid red line indicating the aligned time
            line([0, 0], [ylims(1), ylims(2)], ...
                'LineWidth', Plot_Params.mean_line_width, 'color', 'r');
            % Dotted red line indicating beginning of measured window
            line([-time_before_end, -time_before_end], [ylims(1), ylims(2)], ...
                'LineWidth', Plot_Params.mean_line_width, 'color','r','linestyle','--');
        end
    
        if contains(event, 'window')
            % Dotted purple line indicating beginning of measured window
            line([max_amp_time(ii) - half_window_length, max_amp_time(ii) - half_window_length], ... 
                [ylims(1), ylims(2)], 'linewidth', Plot_Params.mean_line_width, ...
                'color',[.5 0 .5],'linestyle','--');
            % Dotted purple line indicating end of measured window
            line([max_amp_time(ii) + half_window_length, max_amp_time(ii) + half_window_length], ... 
                [ylims(1), ylims(2)], 'linewidth', Plot_Params.mean_line_width, ...
                'color',[.5 0 .5],'linestyle','--');
        elseif ~contains(event, 'trial_gocue') && ~contains(event, 'trial_end')
            % Dotted red line indicating beginning of measured window
            line([-0.1, -0.1], [ylims(1), ylims(2)], ...
                'Linewidth', Plot_Params.mean_line_width, 'Color', 'r', 'Linestyle','--');
            % Dotted red line indicating end of measured window
            line([0.1, 0.1], [ylims(1), ylims(2)], ...
                'Linewidth', Plot_Params.mean_line_width, 'Color', 'r', 'Linestyle','--');
        end

        cc = cc + 1;

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End of the muscle loop

end % End of target loop
