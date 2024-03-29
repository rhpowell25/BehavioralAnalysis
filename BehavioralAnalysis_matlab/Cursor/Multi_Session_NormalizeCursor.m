function Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, signal_choice, norm_cursor)

%% File Description:

% This function finds the nth percentile of cursor position, velocity, or 
% acceleration between two XDS files for the purpose of normalization.
% Depending on the method, the nth percentile will be calculated using
% the entirety of both files or using only the succesful trials.
% If you set norm_cursor to 0, the Cursor_Norm_Factor will be 1
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% signal_choice: 'Pos', 'Vel', or 'Acc'
% norm_cursor: 1 or 0

%% End the function if you are not zeroing the cursor signal
if ~isequal(norm_cursor, 1)
    disp('Cursor Signal Will Not Be Normalized')
    Cursor_Norm_Factor = 1;
    return
else
    norm_prctile = 95;
    fprintf('Cursor Signal Will Be Normalized to the %ith percentile \n', norm_prctile);
end

%% What part of the cursor signal do you want to take the percentile of
% All the cursor data ('All_Cursor')
% The cursor data in each succesful trial ('Trial_Cursor')
norm_method = 'Trial_Cursor';

%% Basic Settings, some variable extractions, & definitions

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

%% Concatenate the cursor signal

cat_curs = cat(1, curs_morn, curs_noon);

%% Find the vector sum of the cursor signal
z_cat_curs = sqrt(cat_curs(:, 2).^2 + cat_curs(:, 1).^2);

%% Finding the max cursor signal throughout both files

if stcmp(norm_method, 'All_Cursor')

    % Find the nth percentile of the cursor position
    Cursor_Norm_Factor = prctile(z_cat_curs, norm_prctile);

end

%% Concatenate the morning & afternoon information
if strcmp(norm_method, 'Trial_Cursor')

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
    
    %% Finding the max cursor signal per trial
    
    max_curs_pertrial = zeros(height(rewarded_curs_sig),1);
    for ii = 1:height(rewarded_curs_sig)
        max_curs_pertrial(ii,1) = max(rewarded_curs_sig{ii,1}(:,1));
    end
    
    %% Find the 95th percentile of the max's
    Cursor_Norm_Factor = prctile(max_curs_pertrial, norm_prctile);
end







