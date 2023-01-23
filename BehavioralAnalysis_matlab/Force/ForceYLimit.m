
function force_YLims = ForceYLimit(xds_morn, xds_noon, norm_force)

%% Setting the Y-Limits to 100 if force is to be normalized

if isequal(norm_force, 1)
    force_YLims = [120, -10];
    return
end

%% Ending the function if there is no force

if strcmp(xds_morn.meta.task, 'WS')
    disp('Event cannot be force related for this task');
    force_YLims = NaN;
    return
end

if xds_morn.has_force == 0
    disp('No force in this file')
    force_YLims = NaN;
    return
end

%% Display the functions being used
disp('Force Y-Limit Function:');

%% Basic Settings, some variable extractions, & definitions
% Percentiles for maximums & minimums
max_perc = 100;
min_perc = 5;

axis_expansion = 10;
    
%% Extracting force

force_morn = xds_morn.force;
force_noon = xds_noon.force;
    
%% Find the sum of force

z_force_morn = force_morn(:, 2) + force_morn(:, 1);
z_force_noon = force_noon(:, 2) + force_noon(:, 1);
 
%% Concatenate the arrays
z_force = cat(1, z_force_morn, z_force_noon);

%% Find the minimums & maximums
force_YLims = zeros(1,2);

% Maximum
force_YLims(1) = prctile(z_force, max_perc) + axis_expansion;
% Minimum
force_YLims(2) = prctile(z_force, min_perc) - axis_expansion;

if strcmp(norm_force, 'Convert')
    % Maximum
    force_YLims(1) = prctile(z_force, max_perc) / 1000*5; % Millivolt conversion * gain
    % Minimum
    force_YLims(2) = prctile(z_force, min_perc) / 1000*5; % Millivolt conversion * gain
end











