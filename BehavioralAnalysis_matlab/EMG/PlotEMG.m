function [avg_EMG] = PlotEMG(xds, event, unit_name, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, Plot_Figs, Save_Figs)

%% Display the function being used
disp('Plot EMG Function:');

%% Find the EMG index
[M] = EMG_Index(xds, muscle_groups);

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

% Font specifications
label_font_size = 15;
legend_font_size = 12;
title_font_size = 15;
plot_line_size = 3;
figure_width = 750;
figure_height = 250;
font_name = 'Arial';

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dirs = length(target_dirs);

% Define the output variable
avg_EMG = struct([]);

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

    %% Extracting EMG during successful trials

    EMG = struct([]); % EMG during each successful trial
    for ii = 1:length(Alignment_Times)
        EMG_idx = find(xds.time_frame == Alignment_Times(ii));
        try
            EMG{ii, 1} = xds.EMG((EMG_idx - (before_event / xds.bin_width)) : (EMG_idx + (after_event / xds.bin_width)), M);
        catch
            EMG{ii, 1} = xds.EMG((EMG_idx - (before_event / xds.bin_width)) : end, M);
            EMG{ii, 1} = cat(1, EMG{ii, 1}, NaN(length(EMG{ii-1, 1}) - length(EMG{ii, 1}), length(M)));
        end
    end
    
    %% Define the absolute timing
    absolute_timing = linspace(-before_event, after_event, length(EMG{1,1}));

    %% Putting all succesful trials in one array
    all_trials_EMG = struct([]);
    for ii = 1:length(M)
        all_trials_EMG{ii,1} = zeros(length(EMG{1,1}),length(Alignment_Times));
        for mm = 1:length(Alignment_Times)
            all_trials_EMG{ii,1}(:,mm) = EMG{mm, 1}(:, ii);
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
    
    %% Calculating average EMG (Average Across trials)

    cross_trial_std_EMG = struct([]);
    for ii = 1:length(M)
        avg_EMG{jj,1}{ii,1} = zeros(length(EMG{1,1}),1);
        cross_trial_std_EMG{ii,1} = zeros(length(EMG{1,1}),1);
        for mm = 1:length(EMG{1,1})
            avg_EMG{jj,1}{ii,1}(mm) = mean(all_trials_EMG{ii,1}(mm,:));
            cross_trial_std_EMG{ii,1}(mm) = std(all_trials_EMG{ii,1}(mm,:));
        end
    end

    %% Plot the mean EMG
    if isequal(Plot_Figs, 1)
        for ii = 1:length(M)
    
            % Find the max & min EMG for the y-axis limits
            y_max = round(max(all_trials_EMG{ii,1}, [], 'All')) + 10;
            y_min = round(min(all_trials_EMG{ii,1}, [], 'All')) - 10;
    
            EMG_figure = figure;
            EMG_figure.Position = [300 300 figure_width figure_height];
            hold on
    
            % Titling the plot
            EMG_title = strrep(string(xds.EMG_names(M(ii))),'EMG_','');
            title(sprintf('Mean EMG, %i°, TgtCenter at %0.1f: %s', ...
                    target_dirs(jj), target_centers(jj), EMG_title), 'FontSize', title_font_size)
    
            % Labels
            ylabel('EMG', 'FontSize', label_font_size);
            xlabel('Time (sec.)', 'FontSize', label_font_size);
    
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
    
            % Mean EMG
            plot(absolute_timing, avg_EMG{jj,1}{ii,1}, ...
                'LineWidth', 2, 'Color', 'k');
    
            % Standard Deviation
            plot(absolute_timing, avg_EMG{jj,1}{ii,1} + cross_trial_std_EMG{ii,1}, ...
                'LineWidth', 1, 'LineStyle','--', 'Color', 'r');
            plot(absolute_timing, avg_EMG{jj,1}{ii,1} - cross_trial_std_EMG{ii,1}, ...
                'LineWidth', 1, 'LineStyle','--', 'Color', 'r');
    
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
    
            legend(sprintf('%s', strrep(string(xds.EMG_names(M(ii))),'EMG_',' ')), ... 
                    'NumColumns', 1, 'FontSize', legend_font_size, 'FontName', font_name, ...
                    'Location', 'NorthEast');
            legend boxoff
    
        end % End of the muscle loop
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
