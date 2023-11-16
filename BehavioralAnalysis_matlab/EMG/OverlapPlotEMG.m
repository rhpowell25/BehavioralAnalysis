function OverlapPlotEMG(xds_morn, xds_noon, event, zero_EMG, norm_EMG, muscle_group, Save_File)

%% Display the function being used
disp('Overlap EMG Function:');

%% Find the EMG index
[M] = EMG_Index(xds_morn, muscle_group);

%% Zero & normalize the EMG
% Zero? (1 = Yes, 0 = No)
zero_method = 'Prctile';
EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_group, zero_method, zero_EMG);

% Normalize EMG? (1 = Yes, 0 = No)
norm_prctile = 99;
EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_group, norm_prctile, norm_EMG);

%% Basic Settings, some variable extractions, & definitions

% Pull the binning paramaters
[Bin_Params] = Binning_Parameters;

% Time before & after the event
before_event = Bin_Params.before_event;
after_event = Bin_Params.after_event;

% Binning information
bin_size = xds_morn.bin_width; % Time (sec.)

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
legend_font_size = 12;
figure_width = 750;
figure_height = 250;

%% Extract the target directions & centers
[target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
[target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

%% Load the average force

avg_EMG_morn = struct([]);
for jj = 1:length(target_dirs_morn)
    % Get the average EMG & peak amplitude time
    [avg_EMG_morn{jj,1}, ~] = ...
        EventWindowEMG(xds_morn, muscle_group, target_dirs_morn(jj), target_centers_morn(jj), event);
end
avg_EMG_noon = struct([]);
for jj = 1:length(target_dirs_noon)
    % Get the average EMG & peak amplitude time
    [avg_EMG_noon{jj,1}, ~] = ...
        EventWindowEMG(xds_noon, muscle_group, target_dirs_noon(jj), target_centers_noon(jj), event);
end

%% Check to see if both sessions use a consistent number of targets

% Find matching targets between the two sessions
[Matching_Idxs_Morn, Matching_Idxs_Noon] = ...
    Match_Targets(target_dirs_morn, target_dirs_noon, target_centers_morn, target_centers_noon);

% Only use the info of target centers conserved between morn & noon
if ~all(Matching_Idxs_Morn) || ~all(Matching_Idxs_Noon)
    disp('Uneven Targets Between Morning & Afternoon');
    target_centers_morn = target_centers_morn(Matching_Idxs_Morn);
    target_dirs_morn = target_dirs_morn(Matching_Idxs_Morn);
    avg_EMG_morn = avg_EMG_morn(Matching_Idxs_Noon);
    avg_EMG_noon = avg_EMG_noon(Matching_Idxs_Noon);
end

%% X-axis
EMG_time = (-before_event:bin_size:after_event);
EMG_time = EMG_time(1:end-1) + bin_size/2;

%% Plot the two overlapped
for jj = 1:length(avg_EMG_morn)

    for ii = 1:height(avg_EMG_morn{jj,1})

        Overlap_EMG_figure = figure;
        Overlap_EMG_figure.Position = [300 300 figure_width figure_height];
        hold on

        %% Zero the average EMG
        avg_EMG_morn{jj,1}(ii,:) = avg_EMG_morn{jj,1}(ii,:) - EMG_Zero_Factor(ii);
        avg_EMG_noon{jj,1}(ii,:) = avg_EMG_noon{jj,1}(ii,:) - EMG_Zero_Factor(ii);
        
        %% Normalize the average EMG
        avg_EMG_morn{jj,1}(ii,:) = (avg_EMG_morn{jj,1}(ii,:) / EMG_Norm_Factor(ii))*100;
        avg_EMG_noon{jj,1}(ii,:) = (avg_EMG_noon{jj,1}(ii,:) / EMG_Norm_Factor(ii))*100;

        %% Plot the average EMG
        plot(EMG_time, avg_EMG_morn{jj,1}(ii,:), 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
        plot(EMG_time, avg_EMG_noon{jj,1}(ii,:), 'LineWidth', 2, 'Color', [.5 0 .5])
        
        y_limits = ylim;
        % Titling the plot
        Fig_Title = strrep(string(xds_morn.EMG_names(M(ii))),'EMG_','');
        title(sprintf('Mean EMG, %i°, TgtCenter at %0.1f: %s', ...
            target_dirs_morn(jj), target_centers_morn(jj), Fig_Title), 'FontSize', title_font_size)
    
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

        %% Save the file if selected
        Save_Figs(Fig_Title, Save_File)

    end
end



