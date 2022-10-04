function EMG_Norm_Factor = Single_Session_NormalizeEMG(xds, muscle_groups, norm_perctile, norm_EMG)

%% File Description:

% This function finds the nth percentile EMG in an XDS 
% file for the purpose of normalization.
% Depending on norm_method, the nth percentile will be calculated using the
% entirety of the file or using only the succesful trials
% If you set norm_EMG to 0, the EMG_Norm_Factor will be 1
%
% -- Inputs --
% xds: the xds file
% muscle_groups: 'Flex', 'Ext', Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All'
% norm_prctile: 0 to 100
% norm_EMG: 1 or 0

%% Find the EMG index

[M] = EMG_Index(xds, muscle_groups);

%% End the function if you are not zeroing the EMG
if ~isequal(norm_EMG, 1)
    disp('EMG Will Not Be Normalized')
    EMG_Norm_Factor = ones(1, length(M));
    return
else
    fprintf('EMG Will Be Normalized to the %ith percentile \n', norm_perctile);
end

%% What part of the EMG do you want to take the percentile of
% All the EMG data ('All_EMG')
% The EMG in each succesful trial ('Trial_EMG')
norm_method = 'All_EMG';

%% Finding the max EMG throughout the file

if strcmp(norm_method, 'All_EMG')

    % Pull the corresponding EMG
    EMG = xds.EMG(:,M);

    EMG_Norm_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        % Find the percentile of the EMG
        EMG_Norm_Factor(1,ii) = prctile(EMG(ii,:), norm_prctile);
    end

end

%% If normalizing to the succesful trials

if strcmp(norm_method, 'Trial_EMG')

    % Bin width
    bin_width = xds.bin_width;
    
    %% Index for rewarded trials in the maximum target direction
    
    total_rewarded_idx = find(xds.trial_result == 'R');
    
    %% Loop to extract only rewarded trials 
    
    % Rewarded start times
    rewarded_start_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_start_time(ii) = xds.trial_start_time(total_rewarded_idx(ii));
    end
    
    % Rewarded go-cues
    rewarded_gocue_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_gocue_time(ii) = xds.trial_gocue_time(total_rewarded_idx(ii));
    end
        
    % Rewarded end times
    rewarded_end_time = zeros(length(total_rewarded_idx),1);
    for ii = 1:length(total_rewarded_idx)
        rewarded_end_time(ii) = xds.trial_end_time(total_rewarded_idx(ii));
    end
    
    %% Pulling the timeframe of the succesful trials
    % Find the rewarded start times in the whole trial time frame
    rewarded_start_idx = zeros(height(rewarded_start_time),1);
    for ii = 1:height(rewarded_start_time)
        rewarded_start_idx(ii) = find(xds.time_frame == rewarded_start_time(ii));
    end
    
    % Find the rewarded end times in the whole trial time frame
    rewarded_end_idx = zeros(height(rewarded_end_time),1);
    for ii = 1:height(rewarded_end_time)
        rewarded_end_idx(ii) = find(xds.time_frame == rewarded_end_time(ii));
    end
    
    %% Pull the EMG corresponding to the extracted time frames
    
    % Pull all the EMG
    EMG = struct([]);
    for ii = 1:height(M)
        for jj = 1:height(total_rewarded_idx)
            if jj == 1
                EMG{ii} = xds.EMG(rewarded_start_idx(jj) : ... 
                    rewarded_end_idx(jj)+(2/bin_width),M(ii));
            else
                EMG{ii} = cat(1, EMG{ii}(:,1), xds.EMG(rewarded_start_idx(jj) : ... 
                    rewarded_end_idx(jj)+(2/bin_width),M(ii)));
            end  
        end
    end

    % Find the percentile of all the EMG
    EMG_Norm_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        EMG_Norm_Factor(ii) = prctile(EMG{ii}, norm_perctile);
    end

end


