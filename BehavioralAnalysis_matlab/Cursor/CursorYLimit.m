
function cursor_YLims = CursorYLimit(xds_morn, xds_noon, signal_choice)

%% File Description:

% This function finds the nth percentile of cursor position, velocity, or 
% acceleration between two XDS files to set the Y-axis limits for plotting.
% Depending on the method, the nth percentile will be calculated using
% the entirety of both files or using only the succesful trials.
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% signal_choice: 'Pos', 'Vel', or 'Acc'

%% Display the functions being used
disp('Cursor Y-Limit Function:');

%% What part of the cursor signal do you want to take the percentile of
% All the cursor data ('All_Cursor')
% The cursor data in each succesful trial ('Trial_Cursor')
YLimit_method = 'Trial_Cursor';

%% Basic Settings, some variable extractions, & definitions
% Percentiles for maximums & minimums
max_perc = 99;
min_perc = 1;

axis_expansion = 0;

% Extract the cursor signal of chhoice
if strcmp(signal_choice, 'Pos')
    curs_morn = xds_morn.curs_p;
    curs_noon = xds_noon.curs_p;
elseif strcmp(signal_choice, 'Vel')
    curs_morn = xds_morn.curs_v;
    curs_noon = xds_noon.curs_v;
elseif strcmp(signal_choice, 'Acc')
    curs_morn = xds_morn.curs_a;
    curs_noon = xds_noon.curs_a;
end

% Initialize the output variable
cursor_YLims = zeros(1,2);
    
%% Concatenate the cursor signal

cat_curs = cat(1, curs_morn, curs_noon);

%% Find the vector sum of the cursor signal
z_cat_curs = sqrt(cat_curs(:, 2).^2 + cat_curs(:, 1).^2);

%% Finding the min & max cursor signal throughout both files

if strcmp(YLimit_method, 'All_Cursor')

    % Maximum
    cursor_YLims(1) = prctile(z_cat_curs, max_perc) + axis_expansion;
    % Minimum
    cursor_YLims(2) = prctile(z_cat_curs, min_perc) - axis_expansion;

end

%% Concatenate the morning & afternoon information
if strcmp(YLimit_method, 'Trial_Cursor')

    % Time frame
    noon_time_frame = xds_noon.time_frame + xds_morn.time_frame(end);
    time_frame = cat(1, xds_morn.time_frame, noon_time_frame);
    
    % Bin width
    bin_width = xds_morn.bin_width;
    
    % Trial results
    trial_results = cat(1, xds_morn.trial_result, xds_noon.trial_result);
    % Trial start times
    noon_trial_start_time = xds_noon.trial_start_time + xds_morn.time_frame(end);
    trial_start_times = cat(1, xds_morn.trial_start_time, noon_trial_start_time);
    % Trial go cue times
    noon_trial_gocue_time = xds_noon.trial_gocue_time + xds_morn.time_frame(end);
    trial_gocue_times = cat(1, xds_morn.trial_gocue_time, noon_trial_gocue_time);
    % Trial end times
    noon_trial_end_time = xds_noon.trial_end_time + xds_morn.time_frame(end);
    trial_end_times = cat(1, xds_morn.trial_end_time, noon_trial_end_time);

    %% Index for rewarded trials in the maximum target direction
    
    total_rewarded_idx = find(trial_results == 'R');

    %% Loop to extract only rewarded trials 
    
    % Rewarded start times
    rewarded_start_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_start_time(ii) = trial_start_times(total_rewarded_idx(ii));
    end
    
    % Rewarded go-cues
    rewarded_gocue_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_gocue_time(ii) = trial_gocue_times(total_rewarded_idx(ii));
    end
        
    % Rewarded end times
    rewarded_end_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_end_time(ii) = trial_end_times(total_rewarded_idx(ii));
    end

    %% Pulling the timeframe of the succesful trials
    % Find the rewarded start times in the whole trial time frame
    rewarded_start_idx = zeros(height(rewarded_start_time),1);
    for ii = 1:height(rewarded_start_time)
        rewarded_start_idx(ii) = find(time_frame == rewarded_start_time(ii));
    end
    
    % Find the rewarded end times in the whole trial time frame
    rewarded_end_idx = zeros(height(rewarded_end_time),1);
    for ii = 1:height(rewarded_end_time)
        rewarded_end_idx(ii) = find(time_frame == rewarded_end_time(ii));
    end
    
    timing = struct([]);
    for ii = 1:height(total_rewarded_idx)
        timing{ii,1} = time_frame(rewarded_start_idx(ii) : ... 
            rewarded_end_idx(ii)+(2/bin_width));
    end

    %% Pull the cursor signal corresponding to the extracted time frames
    
    % Cursor signal measured during each successful trial
    rewarded_curs_sig = struct([]);
    for ii = 1:height(total_rewarded_idx)
        if rewarded_end_idx(ii)+(2/bin_width) > height(z_cat_curs)
            rewarded_curs_sig{ii,1} = z_cat_curs(rewarded_start_idx(ii) : end);
            continue
        end
        rewarded_curs_sig{ii,1} = z_cat_curs(rewarded_start_idx(ii) : ... 
            rewarded_end_idx(ii)+(2/bin_width));
    end
    
    %% Finding the min & max cursor signal per trial
    
    max_curs_pertrial = zeros(height(rewarded_curs_sig),1);
    min_curs_pertrial = zeros(height(rewarded_curs_sig),1);
    for ii = 1:height(rewarded_curs_sig)
        min_curs_pertrial(ii,1) = min(rewarded_curs_sig{ii,1}(:,1));
        max_curs_pertrial(ii,1) = max(rewarded_curs_sig{ii,1}(:,1));
    end
    
    %% Finding the min & max cursor signal throughout both files

    % Maximum
    cursor_YLims(1) = prctile(max_curs_pertrial, max_perc) + axis_expansion;
    % Minimum
    cursor_YLims(2) = prctile(min_curs_pertrial, min_perc) - axis_expansion;

end




