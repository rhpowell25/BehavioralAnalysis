function [avg_EMG, max_amp_time] = EventWindowEMG(xds, muscle_group, target_dir, target_center, event)

%% Find the EMG name of interest
[M] = EMG_Index(xds, muscle_group);

%% Catch possible sources of error
% If there is no unit of that name
if isempty(M)
    fprintf('%s does not exist \n', muscle_group);
    avg_EMG = NaN;
    max_amp_time = NaN;
    return
end

%% Basic settings, some variable extractions, & definitions

% Extract the EMG
EMG = xds.EMG(:, M);

% Extract the time-frame
time_frame = xds.time_frame;

% Pull the binning paramaters
[Bin_Params] = Binning_Parameters;

% Time before & after the event
before_event = Bin_Params.before_event;
after_event = Bin_Params.after_event;

% Binning information
bin_size = xds.bin_width; % Time (sec.)

% Window to calculate max firing rate
half_window_size = Bin_Params.half_window_size; % Bins
%half_window_length = Bin_Params.half_window_length;
step_size = Bin_Params.step_size; % Bins

%% Times for rewarded trials

[rewarded_gocue_time] = EventAlignmentTimes(xds, target_dir, target_center, 'trial_gocue');
[rewarded_end_time] = EventAlignmentTimes(xds, target_dir, target_center, 'trial_end');
[Alignment_Times] = EventAlignmentTimes(xds, target_dir, target_center, event);

%% Getting the EMG based on the behavior timings above

array_length = length(-before_event:bin_size:after_event) - 1;

aligned_EMG = struct([]); % EMG during each successful trial
for ii = 1:length(Alignment_Times)
    EMG_idx = find(time_frame == Alignment_Times(ii));
    try
        aligned_EMG{ii, 1} = EMG((EMG_idx - (before_event / bin_size) + 1) : (EMG_idx + (after_event / bin_size)), :);
    catch
        if (EMG_idx - (before_event / bin_size) + 1) > 0
            aligned_EMG{ii, 1} = EMG((EMG_idx - (before_event / bin_size) + 1) : end);
            aligned_EMG{ii, 1} = cat(1, aligned_EMG{ii, 1}, ...
                NaN(array_length - length(aligned_EMG{ii, 1})), length(M));
        else
            aligned_EMG{ii, 1} = EMG(1 : (EMG_idx + (after_event / bin_size)), :);
            aligned_EMG{ii, 1} = cat(1, NaN(array_length - length(aligned_EMG{ii, 1}), length(M)), aligned_EMG{ii, 1});
        end
    end
end

%% Putting all succesful trials in one array
all_trials_EMG = struct([]);
for ii = 1:length(M)
    all_trials_EMG{ii,1} = zeros(length(Alignment_Times), length(aligned_EMG{1,1}));
    for mm = 1:length(Alignment_Times)
        all_trials_EMG{ii,1}(mm,:) = aligned_EMG{mm, 1}(:, ii);
    end
end

% Averaging the EMG
avg_EMG = zeros(length(M), width(all_trials_EMG{1,1}));
for ii = 1:length(M)
    avg_EMG(ii,:) = mean(all_trials_EMG{ii,1});
end

%% Find the trial lengths
gocue_to_event = Alignment_Times - rewarded_gocue_time;
event_to_end = rewarded_end_time - Alignment_Times;

%% Find the 5th percentile of the trial go cue
max_gocue_to_event = prctile(gocue_to_event, 5);

%% Find the 90th percentile of the trial end
max_event_to_end = prctile(event_to_end, 90);

%% Convert the times to fit the bins
max_gocue_idx = round(max_gocue_to_event / bin_size);

% Start 0.2 seconds after trial gocue
if strcmp(event, 'window_trial_gocue')
    max_gocue_idx = -(0.2 / bin_size);
end

max_end_idx = round(max_event_to_end / bin_size);

%% Loop through each muscle
max_amp_time = zeros(length(M),1);
for ii = 1:length(M)
    %% Calculate the floating average
    % This array starts after the 5th percentile go-cue   
    % and ends after the 90th percentile trial end

    try
        array = avg_EMG(ii, length(avg_EMG(ii,:)) / 2 - max_gocue_idx: ...
            length(avg_EMG) / 2 + max_end_idx);
    catch
        array = avg_EMG(ii, length(avg_EMG(ii,:)) / 2 - max_gocue_idx: ...
            end - half_window_size);
    end
    [float_avg, ~, array_idxs] = Sliding_Window(array, half_window_size, step_size);

    %% Find where the highest average was calculated
    max_float_avg = max(float_avg);
    max_amp_idx = find(float_avg == max_float_avg);
    
    max_array_idxs = array_idxs{max_amp_idx(1)};
    
    center_max_amp_idx = max_array_idxs(ceil(end/2)) + length(avg_EMG) / 2 - max_gocue_idx - 1;
    
    %% Print the maximum EMG amplitude in that window
    EMG_name = strrep(xds.EMG_names(M(ii)), 'EMG_', '');
    fprintf("The max %s amplitude is %0.1f mV \n", char(EMG_name), max_float_avg);
    
    %% Display the measured time window
    max_amp_time(ii) = (-before_event) + (center_max_amp_idx*bin_size);
    fprintf("%s peaks at %0.2f seconds \n", char(EMG_name), max_amp_time(ii));

end





















