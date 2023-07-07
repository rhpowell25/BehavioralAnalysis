function Force_Norm_Factor = Multi_Session_NormalizeForce(xds_morn, xds_noon, norm_force)

%% End the function if you are not normalizing the force
if ~isequal(norm_force, 1)
    if strcmp(norm_force, 'Convert')
        disp('Force Will Be Converted To Cursor Units')
        Force_Norm_Factor = 'Convert'; 
        return
    else
        disp('Force Will Not Be Normalized')
        Force_Norm_Factor = 100;
        return
    end
else
    norm_prctile = 95;
    fprintf('Force Will Be Normalized to the %ith percentile \n', norm_prctile);
end

%% Concatenate all the morning & afternoon information

% Time frame
noon_time_frame = xds_noon.time_frame + xds_morn.time_frame(end);
time_frame = cat(1, xds_morn.time_frame, noon_time_frame);

% Bin width
bin_width = xds_morn.bin_width;

% Force
cat_Force = cat(1, xds_morn.force, xds_noon.force);

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

%% Pull the force corresponding to the extracted time frames

% Pull all the force
Force = struct([]);
for ii = 1:height(total_rewarded_idx)
    Force{ii,1} = cat_Force(rewarded_start_idx(ii) : ... 
        rewarded_end_idx(ii)+(2/bin_width),:);
end

%% Recompose the force
[Sigma_Force] = Sum_Force(xds_morn.meta.task, Force);

%% Finding the max force per trial

max_Force_pertrial = zeros(height(Sigma_Force),1);
for ii = 1:height(Sigma_Force)
    max_Force_pertrial(ii,1) = max(Sigma_Force{ii,1}(:,1));
end

%% Find the 95th percentile of the max's
Force_Norm_Factor = prctile(max_Force_pertrial, norm_prctile);







