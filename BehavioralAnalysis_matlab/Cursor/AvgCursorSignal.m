function [avg_curs_sig] = ...
    AvgCursorSignal(xds, signal_choice, event, unit_name, cursor_YLims, Plot_Figs, Save_File)

%% File Description:

% This function finds the average cursor position, velocity, or 
% acceleration per target direction / distance of all the succesful trials 
% in an xds file. It can also plot that signal.
% If you set Plot_Figs to 0, the figure will not be plotted.
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
% Plot_Figs: 1 or 0
% Save_Figs: 'pdf', 'png', 'fig', or 0

%% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

%% Basic settings, some variable extractions, & definitions

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

% Font & figure specifications
label_font_size = 15;
title_font_size = 15;
plot_line_size = 3;
figure_width = 750;
figure_height = 250;

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dirs = length(target_dirs);

% Define the output variable
avg_curs_sig = struct([]);

%% Begin the loop through all directions
for jj = 1:num_dirs
    
    %% Times for rewarded trials
    if strcmp(event, 'trial_gocue')
        [Alignment_Times] = EventAlignmentTimes(xds, NaN, NaN, event);
    else
        [Alignment_Times] = EventAlignmentTimes(xds, target_dirs(jj), target_centers(jj), event);
    end

    if contains(event, 'window')
        % Run the firing rate window function
        [~, max_fr_time] = ...
        EventWindow(xds, unit_name, target_dirs(jj), target_centers(jj), event);
    end

    %% Extracting the cursor signal during successful trials

    % Cursor signal measured during each successful trial 
    rewarded_curs_sig = struct([]); % Cursor signal during each successful trial 
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
        if temp_start ~= 0
            rewarded_curs_sig{ii,1} = cat(1, NaN(abs(temp_start) + 1, width(rewarded_curs_sig{ii})),rewarded_curs_sig{ii,1});
        end
    end

    %% Find the vector sum of the cursor signal

    z_curs_sig = struct([]);
    for ii = 1:length(Alignment_Times)
        z_curs_sig{ii,1} = zeros(length(rewarded_curs_sig{ii,1}(:,1)),1);
        for dd = 1:length(z_curs_sig{ii,1})
            z_curs_sig{ii,1}(dd,1) = sqrt(rewarded_curs_sig{ii,1}(dd,1).^2 + rewarded_curs_sig{ii,1}(dd,2).^2);
        end
    end

    %% Put all the recomposed cursor signals in a single matrix
    all_z_curs_sig = zeros(length(z_curs_sig{1,1}), length(z_curs_sig));
    for ii = 1:length(z_curs_sig)
        all_z_curs_sig(:,ii) = z_curs_sig{ii,1}(:,1);
    end

    %% Average the recomposed cursor signals
    avg_curs_sig{jj,1} = zeros(length(all_z_curs_sig),1);
    for ii = 1:length(avg_curs_sig{jj,1})
        avg_curs_sig{jj,1}(ii) = mean(all_z_curs_sig(ii,:));
    end

    %% Find the standard dev of the recomposed cursor signals
    std_z_curs_sig = zeros(length(all_z_curs_sig),1);
    for ii = 1:length(avg_curs_sig{jj,1})
        std_z_curs_sig(ii) = std(all_z_curs_sig(ii,:));
    end

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(rewarded_curs_sig{1,1}));
    
    %% Plot the average & standard deviation cursor signal

    if isequal(Plot_Figs, 1)
        
        Cursor_figure = figure;
        Cursor_figure.Position = [300 300 figure_width figure_height];
        hold on
        
        % Average
        plot(absolute_timing, avg_curs_sig{jj,1}, 'k', 'LineWidth', 2)
        % Standard Dev
        plot(absolute_timing, avg_curs_sig{jj,1} + std_z_curs_sig, ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        plot(absolute_timing, avg_curs_sig{jj,1} - std_z_curs_sig, ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        
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
        xlabel('Time (sec.)', 'FontSize', label_font_size);
        
        % Titling the plot
        Fig_Title = sprintf('Mean %s: %i°, TgtCenter at %0.1f', ... 
            signal_label, target_dirs(jj), target_centers(jj));
        title(Fig_Title, 'FontSize', title_font_size)

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End the plot if-statement

    % End the event after one loop if showing baseline firing rate
    if strcmp(event, 'trial_gocue')
        return
    end
    
end % End of target loop


