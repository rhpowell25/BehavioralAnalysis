function [avg_sigma_force] = PlotForce(xds, event, unit_name, Force_Norm_Factor, force_YLims, Plot_Figs, Save_File)

%% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

%% Basic Settings, some variable extractions, & definitions

gain_factor = 5; % 3; %

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

% Font & plotting specifications
[Plot_Params] = Plot_Parameters;

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dirs = length(target_dirs);

% Define the output variable
avg_sigma_force = struct([]);

%% Begin the loop through all directions
for jj = 1:num_dirs
    
    %% Times for rewarded trials
    if strcmp(event, 'trial_gocue')
        [Alignment_Times] = EventAlignmentTimes(xds, NaN, NaN, event);
    else
        [Alignment_Times] = EventAlignmentTimes(xds, target_dirs(jj), target_centers(jj), event);
    end

    if contains(event, 'window')
        % Run the preferred direction window function
        [~, max_fr_time] = ...
        EventWindow(xds, unit_name, target_dirs(jj), target_centers(jj), event);
    end

    %% Extracting force during successful trials
    % Force measured during each successful trial
    Force = struct([]); % Force during each successful trial
    for ii = 1:length(Alignment_Times)
        Force_idx = find(xds.time_frame == Alignment_Times(ii));
        try
            Force{ii,1} = xds.force((Force_idx - (before_event / xds.bin_width)) : (Force_idx + (after_event / xds.bin_width)), :);
        catch
            Force{ii, 1} = xds.force((Force_idx - (before_event / xds.bin_width)) : end, :);
            Force{ii, 1} = cat(1, Force{ii, 1}, NaN(length(Force{ii-1, 1}) - length(Force{ii, 1}), 2));
        end
    end

    %% Normalizing the force
    if ~ischar(Force_Norm_Factor)
        for ii = 1:length(Alignment_Times)
            Force{ii,1} = Force{ii,1} / Force_Norm_Factor*100;
        end
    else
        for ii = 1:length(Alignment_Times)
            Force{ii,1} = Force{ii,1} / 1000*gain_factor; % Millivolt conversion * gain
        end
    end

    %% Sum the two force transducers
    [Sigma_Force] = Sum_Force(xds.meta.task, Force);

    %% Put all the recomposed force info in a single matrix
    all_Sigma_Force = zeros(length(Sigma_Force{1,1}), length(Sigma_Force));
    for ii = 1:length(Sigma_Force)
        all_Sigma_Force(:,ii) = Sigma_Force{ii,1}(:,1);
    end

    %% Average the recomposed force info
    avg_sigma_force{jj,1} = zeros(length(all_Sigma_Force),1);
    for ii = 1:length(avg_sigma_force{jj,1})
        avg_sigma_force{jj,1}(ii) = mean(all_Sigma_Force(ii,:));
    end

    %% Find the standard dev of the recomposed force info

    std_z_force = zeros(length(all_Sigma_Force),1);
    for ii = 1:length(avg_sigma_force{jj,1})
        std_z_force(ii) = std(all_Sigma_Force(ii,:));
    end

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(Sigma_Force{1,1}));
    
    %% Plot the average and standard deviation cursor position

    if isequal(Plot_Figs, 1)
        
        Force_figure = figure;
        Force_figure.Position = [300 300 Plot_Params.fig_size Plot_Params.fig_size / 2];
        hold on
        
        % Average
        plot(absolute_timing, avg_sigma_force{jj,1}, 'k', 'LineWidth', 2)
        % Standard Dev
        plot(absolute_timing, avg_sigma_force{jj,1} + std_z_force, ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        plot(absolute_timing, avg_sigma_force{jj,1} - std_z_force, ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        
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
            line([max_fr_time - half_window_length, max_fr_time - half_window_length], ... 
                [ylims(1), ylims(2)], 'linewidth', Plot_Params.mean_line_width, ...
                'color',[.5 0 .5],'linestyle','--');
            % Dotted purple line indicating end of measured window
            line([max_fr_time + half_window_length, max_fr_time + half_window_length], ... 
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
        
        % Labeling the axis
        ylabel('Force', 'FontSize', Plot_Params.label_font_size);
        xlabel('Time (sec.)', 'FontSize', Plot_Params.label_font_size);
        
        % Titling the top plot
        Fig_Title = sprintf('Mean force: %i°, TgtCenter at %0.1f', ... 
            target_dirs(jj), target_centers(jj));
        title(Fig_Title, 'FontSize', Plot_Params.title_font_size)

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end % End the plot if-statement

    %% End the event after one loop if showing baseline firing rate
    if strcmp(event, 'trial_gocue')
        return
    end
    
end % End of target loop


