function [target_dirs, target_centers, rxn_time, trial_length] = Trial_Length_And_Rxn_Time(xds)

%% File Description:

% This function finds the reaction time (defined as the time after the
% go-cue when the EMG of interest exceeds 2 std of the baseline EMG) and
% the trial length for each succesful trial in each unique target direction
% and target center combination
% The EMG of interest is chosen based on the task / target 
%
% -- Inputs --
% xds: the xds file

%% Basic Settings, some variable extractions, & definitions

% Define period before Go-Cue & end to measure
time_before_gocue = 0.4;

% Find the EMG index
if contains(xds.meta.rawFileName, 'PG')
    muscle_groups = 'Grasp';
    [M] = EMG_Index(xds, muscle_groups);
end

% Extract the target directions & centers
[target_dirs, target_centers] = Identify_Targets(xds);

tgt_Center_header = contains(xds.trial_info_table_header, 'TgtDistance');
tgt_Center_idx = cell2mat(xds.trial_info_table(:, tgt_Center_header));

%% Indexes for rewarded trials in all directions
% Counts the number of directions used
num_dir = length(target_dirs);

%% Begin the loop through all directions
for jj = 1:num_dir

    %% Find the EMG index
    if isequal(xds.meta.task, 'WS')
        if isequal(target_dirs(jj), 0) && strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'Flex';
        end
        if isequal(target_dirs(jj), 90)
            muscle_groups = 'Rad_Dev';
        end
        if isequal(target_dirs(jj), 180) && strcmp(xds.meta.hand, 'Left')
            muscle_groups = 'Exten';
        end
        if isequal(target_dirs(jj), -90)
            muscle_groups = 'Uln_Dev';
        end
        [M] = EMG_Index(xds, muscle_groups);
    end
    
    %% Indexes for rewarded trials
    rewarded_idx = find((xds.trial_result == 'R') & (xds.trial_target_dir == target_dirs(jj)) & ...
        (tgt_Center_idx == target_centers(jj)));

    %% Define the output variables
    if jj == 1
        rxn_time = struct([]);
        trial_length = struct([]);
    end

    %% Loops to extract only rewarded trials
    % End times for succesful trials
    rewarded_end_time = zeros(length(rewarded_idx),1);
    for ii = 1:length(rewarded_idx)
        rewarded_end_time(ii) = xds.trial_end_time(rewarded_idx(ii));
    end

    % Go-cue's for succesful trials
    rewarded_gocue_time = zeros(length(rewarded_idx),1);
    for ii = 1:length(rewarded_idx)
        rewarded_gocue_time(ii) = xds.trial_gocue_time(rewarded_idx(ii));
    end
       
    %% Extracting EMG & time during successful trials

    % EMG & time measured during each successful trial 
    Baseline_EMG = struct([]); % EMG during each successful trial 
    Baseline_timings = struct([]); % Time points during each succesful trial 
    for ii = 1:length(rewarded_idx)
        Baseline_idx = find((xds.time_frame > rewarded_gocue_time(ii) - time_before_gocue) & ...
            (xds.time_frame < rewarded_gocue_time(ii)));
        for mm = 1:length(M)
            Baseline_EMG{ii, 1}(:,mm) = xds.EMG(Baseline_idx, M(mm));
        end
        Baseline_timings{ii, 1} = xds.time_frame(Baseline_idx);
    end

    % EMG & time measured during each successful trial 
    Trial_EMG = struct([]); % EMG during each successful trial 
    Trial_timings = struct([]); % Time points during each succesful trial 
    for ii = 1:length(rewarded_idx)
        Trial_idx = find((xds.time_frame > rewarded_gocue_time(ii)) & ...
            (xds.time_frame < rewarded_end_time(ii)));
        for mm = 1:length(M)
            Trial_EMG{ii, 1}(:,mm) = xds.EMG(Trial_idx, M(mm));
        end
        Trial_timings{ii, 1} = xds.time_frame(Trial_idx);
    end

    %% Find the mean of the EMG
    mean_Baseline_EMG = struct([]);
    for ii = 1:length(Baseline_EMG)
        for mm = 1:length(Baseline_EMG{ii,1})
            mean_Baseline_EMG{ii,1}(mm,1) = mean(Baseline_EMG{ii,1}(mm,:));
        end
    end

    mean_Trial_EMG = struct([]);
    for ii = 1:length(Trial_EMG)
        for mm = 1:length(Trial_EMG{ii,1})
            mean_Trial_EMG{ii,1}(mm,1) = mean(Trial_EMG{ii,1}(mm,:));
        end
    end

    %% Find the standard deviations of the baseline EMG
    std_Baseline_EMG = zeros(length(mean_Baseline_EMG),1);
    for ii = 1:length(mean_Baseline_EMG)
        std_Baseline_EMG(ii) = std(mean_Baseline_EMG{ii,1});
    end

    %% Find the trial EMG that exceeds 2 std of the baseline EMG
    rxn_time_EMG = zeros(length(mean_Trial_EMG),1);
    rxn_time_idx = zeros(length(mean_Trial_EMG),1);
    for ii = 1:length(mean_Trial_EMG)
        rxn_time_EMG_idx = find(mean_Trial_EMG{ii,1} >= ...
            (mean(mean_Baseline_EMG{ii,1}) + 2*std_Baseline_EMG(ii,1)));
        if isempty(rxn_time_EMG_idx) % If nothing exceeds 2 std try 1 std
            rxn_time_EMG_idx = find(mean_Trial_EMG{ii,1} >= ...
                (mean(mean_Baseline_EMG{ii,1}) + std_Baseline_EMG(ii,1)));
            if isempty(rxn_time_EMG_idx) % If nothing exceeds 1 std just find the max of that trial
                rxn_time_EMG_idx = find(max(mean_Baseline_EMG{ii,1}));
            end
        end
        rxn_time_idx(ii) = rxn_time_EMG_idx(1);
        rxn_time_EMG(ii) = Trial_timings{ii,1}(rxn_time_idx(ii));
    end

    %% Test plot to make sure everthing works 
    %figure
    %hold on
    %for ii = 1:length(Baseline_timings)
    %    plot(Baseline_timings{ii,1}, mean_Baseline_EMG{ii,1}, 'Color', 'g')
    %    plot(Trial_timings{ii,1}, mean_Trial_EMG{ii,1}, 'Color', 'r')
    %    scatter(rxn_time_EMG(ii), mean_Trial_EMG{ii,1}(rxn_time_idx(ii),1))
    %end

    %% Defining the output variables

    % Reaction time
    rxn_time{jj,1} = rxn_time_EMG - rewarded_gocue_time;
    % Trial Length
    trial_length{jj,1} = rewarded_end_time - rewarded_gocue_time;

end % End of target loop



