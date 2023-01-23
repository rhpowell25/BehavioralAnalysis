function Cursor_Norm_Factor = Multi_Session_NormalizeCursor(xds_morn, xds_noon, norm_cursor)


%% End the function if you are not zeroing the cursor position
if ~isequal(norm_cursor, 1)
    disp('Cursor Will Not Be Normalized')
    Cursor_Norm_Factor = 1;
    return
else
    norm_prctile = 95;
    fprintf('Cursor Will Be Normalized to the %ith percentile \n', norm_prctile);
end

%% Concatenate all the morning & afternoon information

% Time frame
noon_time_frame = xds_noon.time_frame + xds_morn.time_frame(end);
time_frame = cat(1, xds_morn.time_frame, noon_time_frame);

% Bin width
bin_width = xds_morn.bin_width;

% Cursor
cat_Cursor = cat(1, xds_morn.curs_p, xds_noon.curs_p);

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

%% Pull the cursor position corresponding to the extracted time frames

% Pull all the cursor position
cursor_p = struct([]);
for ii = 1:height(total_rewarded_idx)
    if rewarded_end_idx(ii)+(2/bin_width) > height(cat_Cursor)
        cursor_p{ii,1} = cat_Cursor(rewarded_start_idx(ii) : end,:);
        continue
    end
    cursor_p{ii,1} = cat_Cursor(rewarded_start_idx(ii) : ... 
        rewarded_end_idx(ii)+(2/bin_width),:);
end

%% Recompose the cursor position
z_Cursor = struct([]);
% Loops through cursor position
for ii = 1:height(total_rewarded_idx)
    z_Cursor{ii,1} = sqrt(cursor_p{ii,1}(:, 2).^2 + cursor_p{ii, 1}(:, 1).^2);
end

%% Finding the max cursor position per trial

max_Cursor_pertrial = zeros(height(z_Cursor),1);
for ii = 1:height(z_Cursor)
    max_Cursor_pertrial(ii,1) = max(z_Cursor{ii,1}(:,1));
end

%% Find the 95th percentile of the max's
Cursor_Norm_Factor = prctile(max_Cursor_pertrial, norm_prctile);







