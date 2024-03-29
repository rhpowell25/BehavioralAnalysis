function [mp_amp, std_mp, pertrial_mp_amp] = EventPeakEMG(xds, muscle_group, event)

%% Find the EMG name of interest
[M] = EMG_Index(xds, muscle_group);

%% Catch possible sources of error
% If there is no unit of that name
if isempty(M)
    fprintf('%s does not exist \n', muscle_group);
    mp_amp = NaN;
    std_mp = NaN;
    pertrial_mp_amp = NaN;
    return
end

%% Basic settings, some variable extractions, & definitions

% Extract the EMG
EMG = xds.EMG(:, M);

% Extract the time-frame
time_frame = xds.time_frame;

% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

% Pull the binning paramaters
[Bin_Params] = Binning_Parameters;

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dir = length(target_dirs);

%% Begin the loop through all directions
for jj = 1:num_dir
    
    %% Times for rewarded trials
    [Alignment_Times] = EventAlignmentTimes(xds, target_dirs(jj), target_centers(jj), event);

    %% Time period for peak EMG amplitude
    if contains(event, 'window')
        [~, max_amp_time] = EventWindowEMG(xds, muscle_group, target_dirs(jj), target_centers(jj), event);
        % Window to calculate max firing rate
        half_window_length = Bin_Params.half_window_length; % Time (sec.)
    else
        max_amp_time = 0; % Time (sec.)
        half_window_length = 0.1; % Time (sec.)
    end

    %% Define the output variables
    if jj == 1
        mp_amp = zeros(num_dir, length(M));
        std_mp = zeros(num_dir, length(M));
        pertrial_mp_amp = struct([]);
    end

    %% Peak EMG amplitude
    for ii = 1:length(Alignment_Times)
        for mm = 1:length((M))
            t_start = Alignment_Times(ii) + max_amp_time(mm) - half_window_length;
            t_end = Alignment_Times(ii) + max_amp_time(mm) + half_window_length;
            EMG_idx = (time_frame >= t_start) & (time_frame <= t_end);
            pertrial_mp_amp{jj,1}(ii,mm) = mean(EMG(EMG_idx,mm));
        end
    end

    %% Defining the output variables
    for mm = 1:length((M))
        % Peak EMG amplitude
        mp_amp(jj,mm) = mean(pertrial_mp_amp{jj,1}(ii,:));
        % Standard Deviation
        std_mp(jj,mm) = std(pertrial_mp_amp{jj,1}(ii,:));
    end
    
end % End of target loop



