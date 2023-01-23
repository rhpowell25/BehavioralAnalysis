function OverlapPlotEMG(xds_morn, xds_noon, event, zero_EMG, norm_EMG, muscle_groups, Save_Figs)

%% Display the function being used
disp('Overlap EMG Function:');

%% Zero & normalize the EMG
% Zero? (1 = Yes, 0 = No)
zero_method = 'Prctile';
EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG);

% Normalize EMG? (1 = Yes, 0 = No)
norm_prctile = 99;
EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_prctile, norm_EMG);

%% Find the EMG index
[M] = EMG_Index(xds_morn, muscle_groups);

%% Basic Settings, some variable extractions, & definitions

% Define the window for the baseline phase
time_before_gocue = 0.4;
% Define the window for the movement phase
time_before_end = xds_morn.meta.TgtHold;

% Times displayed in the raster
before_event = 3.0;
after_event = 3.0;

% Font & figure specifications
label_font_size = 15;
title_font_size = 15;
legend_font_size = 12;
figure_width = 750;
figure_height = 250;

%% Load the average force

% Average EMG
[average_EMG_morn] = PlotEMG(xds_morn, event, NaN, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0, 0);
[average_EMG_noon] = PlotEMG(xds_noon, event, NaN, EMG_Zero_Factor, EMG_Norm_Factor, muscle_groups, 0, 0);

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
    average_EMG_morn = average_EMG_morn(Matching_Idxs_Noon);
    average_EMG_noon = average_EMG_noon(Matching_Idxs_Noon);
end

%% X-axis
EMG_time = linspace(-before_event, after_event, length(average_EMG_morn{1}{1}));

%% Y-axis limits
%y_limits = zeros(2,1);
%y_max = zeros(length(average_EMG_morn)*length(M),2);
%y_min = zeros(length(average_EMG_morn)*length(M),2);

%cc = 1;
%for jj = 1:length(average_EMG_morn)
%    for pp = 1:length(M)
%        y_max(cc,1) = max(average_EMG_morn{jj,1}{pp,1});
%        y_max(cc,2) = max(average_EMG_noon{jj,1}{pp,1});
%        y_min(cc,1) = min(average_EMG_morn{jj,1}{pp,1});
%        y_min(cc,2) = min(average_EMG_noon{jj,1}{pp,1});
%        cc = cc + 1;
%    end
%end

%y_limits(1) = min(y_min, [],'all') - 0.5;
%y_limits(2) = max(y_max, [],'all') + 0.5;

%% Plot the two overlapped
for jj = 1:length(average_EMG_morn)

    for ii = 1:length(average_EMG_morn{jj,1})

        Overlap_EMG_figure = figure;
        Overlap_EMG_figure.Position = [300 300 figure_width figure_height];
    
        hold on
        plot(EMG_time, average_EMG_morn{jj,1}{ii,1}, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
        plot(EMG_time, average_EMG_noon{jj,1}{ii,1}, 'LineWidth', 2, 'Color', [.5 0 .5])
        
        y_limits = ylim;
        % Titling the plot
        EMG_title = strrep(string(xds_morn.EMG_names(M(ii))),'EMG_','');
        title(sprintf('Mean EMG, %i°, TgtCenter at %0.1f: %s', ...
            target_dirs_morn(jj), target_centers_morn(jj), EMG_title), 'FontSize', title_font_size)
    
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
            
        % Setting the y-axis limits
        %ylim([y_limits(1), y_limits(2)])
        
        ylabel('EMG', 'FontSize', label_font_size);
        xlabel('Time (sec.)', 'FontSize', label_font_size);
        
        % Remove the box of the plot
        box off
    
        % Define the legend location
        if contains(event, 'gocue')
            legend_location = 'NorthEast';
        else
            legend_location = 'NorthWest';
        end
    
        % Setting the x-axis limits
        if contains(event, 'gocue')
            xlim([-before_event + 2, after_event]);
        elseif contains(event, 'end')
            xlim([-before_event, after_event - 2]);
        else
            xlim([-before_event + 1, after_event - 1]);
        end
    
        % Legend
        legend('Morning', 'Afternoon', 'Location', legend_location, 'FontSize', legend_font_size)
    
        % Remove the legend's outline
        legend boxoff

    end
end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for jj = 1:numel(findobj('type','figure'))
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


