function PlotEMG(xds, muscle_group, event, EMG_Zero_Factor, EMG_Norm_Factor, max_EMGYLim, Save_Figs)

%% Display the function being used
disp('Plot EMG Function:');

%% Find the EMG index
[M] = EMG_Index(xds, muscle_group);

%% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

%% Load the average force
avg_EMG = struct([]);
max_amp_time = struct([]);
for jj = 1:length(target_dirs)
    % Get the average EMG & peak amplitude time
    [avg_EMG{jj,1}, max_amp_time{jj,1}] = EventWindowEMG(xds, muscle_group, target_dirs(jj), target_centers(jj), event);
end

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

% Define the figure titles
EMG_title = strings;

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

% Font specifications
label_font_size = 15;
legend_font_size = 12;
title_font_size = 15;
plot_line_size = 3;
figure_width = 750;
figure_height = 250;
font_name = 'Arial';

%% X-axis
EMG_time = (-before_event:bin_size:after_event);
EMG_time = EMG_time(1:end-1) + bin_size/2;

%% Plot the average EMG
cc = 1;

for jj = 1:length(avg_EMG)

    for ii = 1:length(M)

        EMG_name = strrep(xds.EMG_names(M(ii)), 'EMG_', '');

        EMG_figure = figure;
        EMG_figure.Position = [300 300 figure_width figure_height];
        hold on

        %% Zero the average EMG
        avg_EMG{jj,1}(ii,:) = avg_EMG{jj,1}(ii,:) - EMG_Zero_Factor(ii);
        
        %% Normalize the average EMG
        avg_EMG{jj,1}(ii,:) = (avg_EMG{jj,1}(ii,:) / EMG_Norm_Factor(ii))*100;

        %% Plot the average EMG

        plot(EMG_time, avg_EMG{jj,1}(ii,:), 'LineWidth', plot_line_size, 'Color', 'k');

        %% Set the title, labels, axes, & plot lines indicating alignment

        % Titling the plot
        if strcmp(event, 'trial_gocue')
            EMG_title{cc} = strcat(char(EMG_name), {' '}, 'aligned to trial gocue:');
        else
            if contains(event, 'window')
                temp_event = strrep(event, 'window_', '');
            else
                temp_event = event;
            end
            event_title = strcat('aligned to ', {' '}, strrep(temp_event, '_', {' '}), ':');
            EMG_title{cc} = char(strcat(char(EMG_name), {' '}, event_title, {' '}, num2str(target_dirs(jj)), ...
                '°, TgtCenter at', {' '}, num2str(target_centers(jj))));
        end
        if contains(xds.meta.rawFileName, 'Pre')
            EMG_title{cc} = strcat(EMG_title{cc}, ' (Morning)');
        end
        if contains(xds.meta.rawFileName, 'Post')
            EMG_title{cc} = strcat(EMG_title{cc}, ' (Afternoon)');
        end
        title(EMG_title{cc}, 'FontSize', title_font_size)

        % Labels
        ylabel('EMG Amplitude', 'FontSize', label_font_size);
        xlabel('Time (sec.)', 'FontSize', label_font_size);

        % Setting the x-axis limits
        if contains(event, 'gocue') || contains(event, 'onset')
            xlim([-1, after_event]);
        elseif contains(event, 'end')
            xlim([-before_event, 1]);
        else
            xlim([-before_event, after_event]);
        end
        % Setting the y-axis limits
        ylim([0, max_EMGYLim(ii)])
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
            line([max_amp_time{jj}(ii) - half_window_length, max_amp_time{jj}(ii) - half_window_length], ... 
                [ylims(1), ylims(2)], 'linewidth', plot_line_size,'color',[.5 0 .5],'linestyle','--');
            % Dotted purple line indicating end of measured window
            line([max_amp_time{jj}(ii) + half_window_length, max_amp_time{jj}(ii) + half_window_length], ... 
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

        cc = cc + 1;

    end % End of the muscle loop

end % End of target loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = length(findobj('type','figure')):-1:1
        EMG_title{ii} = strrep(EMG_title{ii}, ':', '');
        EMG_title{ii} = strrep(EMG_title{ii}, 'vs.', 'vs');
        EMG_title{ii} = strrep(EMG_title{ii}, 'mg.', 'mg');
        EMG_title{ii} = strrep(EMG_title{ii}, 'kg.', 'kg');
        EMG_title{ii} = strrep(EMG_title{ii}, '.', '_');
        EMG_title{ii} = strrep(EMG_title{ii}, '/', '_');
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(EMG_title{ii})), 'png')
            saveas(gcf, fullfile(save_dir, char(EMG_title{ii})), 'pdf')
            saveas(gcf, fullfile(save_dir, char(EMG_title{ii})), 'fig')
        else
            saveas(gcf, fullfile(save_dir, char(EMG_title{ii})), Save_Figs)
        end
        close gcf
    end
end
