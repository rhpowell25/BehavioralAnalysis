function EMG_Zero_Factor = Single_Session_EMG_Zero(xds, muscle_groups, zero_method, zero_EMG)

%% File Description:

% This function finds the minimum EMG in an XDS file for the purpose of zeroing.
% Depending on the method, the minima will be calculated by a 200 ms.
% moving average ('Window') or the 5th percentile of the EMG. ('Prctile')
% If you set zero_EMG to 0, the EMG_Zero_Factor will also be 0
%
% -- Inputs --
% xds: the xds file
% muscle_groups: 'Flex', 'Ext', Uln_Dev', 'Rad_Dev', 'Both', 'Grasp', 'Custom', or 'All'
% zero_method: 'Window', or 'Prctile'
% zero_EMG: 1 or 0

%% Find the EMG index

[M] = EMG_Index(xds, muscle_groups);

%% End the function if you are not zeroing the EMG
if ~isequal(zero_EMG, 1)
    disp('EMG Will Not Be Zeroed')
    EMG_Zero_Factor = zeros(1, length(M));
    return
else
    disp('EMG Will Be Zeroed')
end

%% Extract the EMG
EMG = zeros(length(xds.time_frame), length(M));
for ii = 1:length(M)
    EMG(:,ii) = xds.EMG(:,M(ii));
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





