function [avg_z_force] = PlotForce(xds, event, unit_name, Force_Norm_Factor, force_YLims, Plot_Figs, Save_Figs)

%% Load the excel file
if ~isnan(unit_name)
    if ~ischar(unit_name)
    
        [xds_output] = Find_Excel(xds);
    
        %% Find the unit of interest
    
        unit = xds_output.unit_names(unit_name);
    
        %% Identify the index of the unit
        N = find(strcmp(xds.unit_names, unit));
    
    else
        N = find(strcmp(xds.unit_names, unit_name));
    end

    %% End the function with NaN output variables if the unit doesnt exist
    if isempty(N)
        fprintf('%s does not exist \n', unit_name);
        return
    end
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
avg_z_force = struct([]);

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
        [~, max_fr_time, ~] = ...
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
            Force{ii,1} = Force{ii,1} / 1000*5; % Millivolt conversion * gain
        end
    end

    %% Sum the two force transducers
    z_Force = struct([]);
    % Loops through force
    for ii = 1:length(Alignment_Times)
        z_Force{ii,1} = Force{ii,1}(:, 2) + Force{ii, 1}(:, 1);
    end

    %% Put all the recomposed cursor positions info in a single matrix
    all_z_Force = zeros(length(z_Force{1,1}), length(z_Force));
    for ii = 1:length(z_Force)
        all_z_Force(:,ii) = z_Force{ii,1}(:,1);
    end

    %% Average the recomposed cursor positions info
    avg_z_force{jj,1} = zeros(length(all_z_Force),1);
    for ii = 1:length(avg_z_force{jj,1})
        avg_z_force{jj,1}(ii) = mean(all_z_Force(ii,:));
    end

    %% Find the standard dev of the recomposed cursor position info

    std_z_force = zeros(length(all_z_Force),1);
    for ii = 1:length(avg_z_force{jj,1})
        std_z_force(ii) = std(all_z_Force(ii,:));
    end

    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(z_Force{1,1}));
    
    %% Plot the average and standard deviation cursor position

    if isequal(Plot_Figs, 1)
        
        Force_figure = figure;
        Force_figure.Position = [300 300 figure_width figure_height];
        hold on
        
        % Average
        plot(absolute_timing, avg_z_force{jj,1}, 'k', 'LineWidth', 2)
        % Standard Dev
        plot(absolute_timing, avg_z_force{jj,1} + std_z_force, ...
            'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
        plot(absolute_timing, avg_z_force{jj,1} - std_z_force, ...
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
        ylabel('Force', 'FontSize', label_font_size);
        xlabel('Time (sec.)', 'FontSize', label_font_size);
        
        % Titling the top plot
        title(sprintf('Mean force: %i°, TgtCenter at %0.1f', ... 
            target_dirs(jj), target_centers(jj)), 'FontSize', title_font_size)

    end % End the plot if-statement

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


