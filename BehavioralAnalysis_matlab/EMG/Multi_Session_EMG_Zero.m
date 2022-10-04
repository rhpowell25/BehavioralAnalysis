function EMG_Zero_Factor = Multi_Session_EMG_Zero(xds_morn, xds_noon, muscle_groups, zero_method, zero_EMG)

%% File Description:

% This function finds the minimum EMG between two concatenated XDS 
% files for the purpose of zeroing.
% Depending on the method, the minima will be calculated by a 200 ms.
% moving average ('Window') or the 5th percentile of the EMG. ('Prctile')
% If you set zero_EMG to 0, the EMG_Zero_Factor will also be 0
%
% -- Inputs --
% xds_morn: the first xds file
% xds_noon: the second xds file
% muscle_groups: 'Flex', 'Ext', Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All'
% zero_method: 'Window', or 'Prctile'
% zero_EMG: 1 or 0

%% Find the EMG index

[M] = EMG_Index(xds_morn, muscle_groups);

%% End the function if you are not zeroing the EMG
if ~isequal(zero_EMG, 1)
    disp('EMG Will Not Be Zeroed')
    EMG_Zero_Factor = zeros(1, length(M));
    return
else
    disp('EMG Will Be Zeroed')
end

%% Extract the EMG
% Morning
EMG_morn = zeros(length(xds_morn.time_frame), length(M));
for ii = 1:length(M)
    EMG_morn(:,ii) = xds_morn.EMG(:,M(ii));
end

% Afternoon
EMG_noon = zeros(length(xds_noon.time_frame), length(M));
for ii = 1:length(M)
    EMG_noon(:,ii) = xds_noon.EMG(:,M(ii));
end

%% Concatenate the morning & afternoon EMG

EMG = zeros(length(xds_morn.time_frame) + length(xds_noon.time_frame), length(M));
for ii = 1:length(M)
    EMG(:,ii) = cat(1, EMG_morn(:, ii), EMG_noon(:, ii));
end

%% Run the moving average
if strcmp(zero_method, 'Window')

    % Moving average window size (200 ms)
    window_size = 200;

    % Window step size (in indices)
    step_size = 1;

    EMG_Zero_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        [sliding_avg, ~, ~] = Sliding_Window(EMG(:, ii), window_size, step_size);
        EMG_Zero_Factor(1,ii) = min(sliding_avg);
    end

end

%% Find the 5th percentile
if strcmp(zero_method, 'Prctile')
    
    EMG_Zero_Factor = zeros(1,length(M));
    for ii = 1:length(M)
        EMG_Zero_Factor(ii) = prctile(EMG(:,ii), 5);
    end
end





