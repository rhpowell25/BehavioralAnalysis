function EMG_Norm_Factor = Multi_Session_NormalizeEMG(xds_morn, xds_noon, muscle_groups, norm_prctile, norm_EMG)

%% File Description:

% This function finds the nth percentile EMG between two concatenated XDS 
% files for the purpose of normalization.
% Depending on norm_method, the nth percentile will be calculated using the
% entirety of both files or using only the succesful trials.
% If you set norm_EMG to 0, the EMG_Norm_Factor will be 1
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% muscle_groups: 'Flex', 'Ext', Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All'
% norm_prctile: 0 to 100
% norm_EMG: 1 or 0

%% Find the EMG index

[M] = EMG_Index(xds_morn, muscle_groups);

%% End the function if you are not normalizing the EMG
if ~isequal(norm_EMG, 1)
    disp('EMG Will Not Be Normalized')
    EMG_Norm_Factor = ones(1, length(M));
    return
else
    fprintf('EMG Will Be Normalized to the %ith percentile \n', norm_prctile);
end

%% What part of the EMG do you want to take the percentile of
% All the EMG data ('All_EMG')
% The EMG in each succesful trial ('Trial_EMG')
norm_method = 'Trial_EMG';

%% Concatenate the EMG

cat_EMG = cat(1, xds_morn.EMG(:,M), xds_noon.EMG(:,M));

%% Finding the max EMG throughout both files

if strcmp(norm_method, 'All_EMG')

    EMG_Norm_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        % Find the percentile of the EMG
        EMG_Norm_Factor(1,ii) = prctile(cat_EMG(ii,:), norm_prctile);
    end

end

%% Concatenate all the morning & afternoon trial related information
if strcmp(norm_method, 'Trial_EMG')

    % Time frame
    noon_time_frame = xds_noon.time_frame + xds_morn.time_frame(end);
    time_frame = cat(1, xds_morn.time_frame, noon_time_frame);
    
    % Bin width
    bin_width = xds_morn.bin_width;

    % Trial results
    trial_results = cat(1, xds_morn.trial_result, xds_noon.trial_result);
    % Trial start times
    noon_trial_start_time = xds_noon.trial_start_time + xds_morn.time_frame(end);
    trial_start_time = cat(1, xds_morn.trial_start_time, noon_trial_start_time);
    % Trial go cue times
    noon_trial_gocue_time = xds_noon.trial_gocue_time + xds_morn.time_frame(end);
    trial_gocue_time = cat(1, xds_morn.trial_gocue_time, noon_trial_gocue_time);
    % Trial end times
    noon_trial_end_time = xds_noon.trial_end_time + xds_morn.time_frame(end);
    trial_end_time = cat(1, xds_morn.trial_end_time, noon_trial_end_time);

    %% Index for rewarded trials
    
    total_rewarded_idx = find(trial_results == 'R');
    
    %% Loop to extract only rewarded trials 
    
    % Rewarded start times
    rewarded_start_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_start_time(ii) = trial_start_time(total_rewarded_idx(ii));
    end
    
    % Rewarded go-cues
    rewarded_gocue_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_gocue_time(ii) = trial_gocue_time(total_rewarded_idx(ii));
    end
        
    % Rewarded end times
    rewarded_end_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_end_time(ii) = trial_end_time(total_rewarded_idx(ii));
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

    %% Pull the EMG corresponding to the extracted time frames

    % Pull all the EMG
    EMG = struct([]);
    for ii = 1:height(M)
        for jj = 1:height(total_rewarded_idx)
            if jj == 1
                EMG{ii} = cat_EMG(rewarded_start_idx(jj) : ... 
                    rewarded_end_idx(jj)+(2/bin_width),ii);
            else
                if rewarded_end_idx(jj)+(2/bin_width) < length(cat_EMG)
                    EMG{ii} = cat(1, EMG{ii}(:,1), cat_EMG(rewarded_start_idx(jj) : ... 
                        rewarded_end_idx(jj)+(2/bin_width),ii));
                else
                    EMG{ii} = cat(1, EMG{ii}(:,1), cat_EMG(rewarded_start_idx(jj) : ... 
                        end,ii));
                end
            end  
        end
    end

    % Find the percentile of all the EMG
    EMG_Norm_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        EMG_Norm_Factor(ii) = prctile(EMG{ii}, norm_prctile);
    end
end







