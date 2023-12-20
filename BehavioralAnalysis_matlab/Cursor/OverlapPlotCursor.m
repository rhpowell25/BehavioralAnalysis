function OverlapPlotCursor(xds_morn, xds_noon, signal_choice, event, Save_File)

%% File Description:

% This function plots the average cursor position, velocity, or 
% acceleration of two XDS files overtop each other for comparison.
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% signal_choice: 'Pos', 'Vel', or 'Acc'
% event: 'trial_gocue', 'window_trial_gocue', 'trial_end', 
% 'window_trial_end', 'force_onset', 'window_force_onset', 'force_max', 
% 'window_force_max', 'window_force_deriv', 'force_deriv', 'cursor_onset', 
% 'window_cursor_onset', 'cursor_veloc', 'window_cursor_veloc', 
% 'cursor_acc', 'window_cursor_acc', 'EMG_max', 'window_EMG_max', 
% 'task_onset', or 'window_task_onset'
% Save_Figs: 'pdf', 'png', 'fig', or 0

%% Display the function being used
disp('Overlap Cursor Function:');

%% Collect the Y-Limits

cursor_YLims = CursorYLimit(xds_morn, xds_noon, signal_choice);

% End the function if there is no Y-Limit
if isnan(cursor_YLims)
    disp("There is no Y-Limit")
    return
end

%% Basic Settings, some variable extractions, & definitions

% Pull the binning paramaters
[Bin_Params] = Binning_Parameters;

% Time before & after the event
before_event = Bin_Params.before_event;
after_event = Bin_Params.after_event;

if contains(event, 'gocue') || contains(event, 'force_onset')
    % Define the window for the baseline phase
    time_before_gocue = 0.4;
elseif contains(event, 'end')
    % Define the window for the movement phase
    time_before_end = xds_morn.meta.TgtHold;
end

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

%% Load the average cursor signals

[avg_curs_sig_morn] = AvgCursorSignal(xds_morn, signal_choice, event, NaN, cursor_YLims, 0, 0);

[avg_curs_sig_noon] = AvgCursorSignal(xds_noon, signal_choice, event, NaN, cursor_YLims, 0, 0);

%% Extract the target directions & centers
[target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
[target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

%% Check to see if both sessions use a consistent number of targets

% Find matching targets between the two sessions
[Matching_Idxs_Morn, Matching_Idxs_Noon] = ...
    Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon);

% Only use the info of target centers conserved between morn & noon
if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
    disp('Uneven Targets Between Morning & Afternoon');
    target_centers_morn = target_centers_morn(Matching_Idxs_Morn);
    target_dirs_morn = target_dirs_morn(Matching_Idxs_Morn);
    avg_curs_sig_morn = avg_curs_sig_morn(Matching_Idxs_Noon);
    avg_curs_sig_noon = avg_curs_sig_noon(Matching_Idxs_Noon);
end

%% X-axis
cursor_time = linspace(-before_event, after_event, length(avg_curs_sig_morn{1}));

%% Y-axis limits
y_limits = zeros(2,1);
y_max = zeros(length(avg_curs_sig_morn),2);
y_min = zeros(length(avg_curs_sig_morn),2);

for ii = 1:length(avg_curs_sig_morn)
    y_max(ii,1) = max(avg_curs_sig_morn{ii,1});
    y_max(ii,2) = max(avg_curs_sig_noon{ii,1});
    y_min(ii,1) = min(avg_curs_sig_morn{ii,1});
    y_min(ii,2) = min(avg_curs_sig_noon{ii,1});
end

y_limits(1) = min(y_min, [],'all') - 0.125;
y_limits(2) = max(y_max, [],'all') + 0.125;

%% Plot the two overlapped
for ii = 1:length(avg_curs_sig_morn)

    Overlap_Cursor_figure = figure;
    Overlap_Cursor_figure.Position = [300 300 Plot_Params.fig_size Plot_Params.fig_size / 2];

    hold on
    plot(cursor_time, avg_curs_sig_morn{ii,1}, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
    plot(cursor_time, avg_curs_sig_noon{ii,1}, 'LineWidth', 2, 'Color', [.5 0 .5])

    if contains(event, 'gocue')
        % Dotted green line indicating beginning of measured window
        line([-time_before_gocue, -time_before_gocue], [y_limits(1), y_limits(2)], ...
            'linewidth',2,'color',[0 0.5 0],'linestyle','--');
        % Solid green line indicating the aligned time
        line([0, 0], [y_limits(1), y_limits(2)], ...
            'linewidth', 2, 'color', [0 0.5 0]);
    end
    if contains(event, 'trial_end')
        % Solid red line indicating the aligned time
        line([0, 0], [y_limits(1), y_limits(2)], ...
            'linewidth', 2, 'color', 'r');
        % Dotted red line indicating beginning of measured window
        line([-time_before_end, -time_before_end], [y_limits(1), y_limits(2)], ...
            'linewidth',2,'color','r','linestyle','--');
    end
        
    % Setting the axis limits
    if contains(event, 'gocue')
        xlim([-before_event + 2, after_event]);
    elseif contains(event, 'end')
        xlim([-before_event, after_event - 2]);
    else
        xlim([-before_event + 1, after_event - 1]);
    end
    ylim([y_limits(2), y_limits(1)]);

    % Define the labels
    if strcmp(signal_choice, 'Pos')
        signal_label = 'Wrist Position';
    elseif strcmp(signal_choice, 'Vel')
        signal_label = 'Wrist Velocity';
    elseif strcmp(signal_choice, 'Acc')
        signal_label = 'Wrist Acceleration';
    end
    
    % Labeling the axis
    ylabel(signal_label, 'FontSize', Plot_Params.label_font_size);
    xlabel('Time (sec.)', 'FontSize', Plot_Params.label_font_size);
    
    % Titling the plot
    Fig_Title = sprintf('Mean %s: %i°, TgtCenter at %0.1f', ... 
        signal_label, target_dirs_morn(ii), target_centers_morn(ii));
    title(Fig_Title, 'FontSize', Plot_Params.title_font_size)

    % Remove the box of the plot
    box off

    % Define the legend location
    if contains(event, 'gocue')
        legend_location = 'NorthEast';
    else
        legend_location = 'NorthWest';
    end

    % Legend
    legend('Morning', 'Afternoon', 'Location', legend_location, 'FontSize', Plot_Params.legend_size)

    % Remove the legend's outline
    legend boxoff 

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)
 
end


