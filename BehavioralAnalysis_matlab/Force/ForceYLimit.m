
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

gain_factor = 5; % 3; % 

% Percentiles for maximums & minimums
max_perc = 100;
min_perc = 5;

axis_expansion = 0;
    
%% Extracting force

Force_morn = xds_morn.force;
Force_noon = xds_noon.force;
    
%% Recompose the force
[Sigma_Force_morn] = Sum_Force(xds_morn.meta.task, {Force_morn});
Sigma_Force_morn = Sigma_Force_morn{1,1};
[Sigma_Force_noon] = Sum_Force(xds_noon.meta.task, {Force_noon});
Sigma_Force_noon = Sigma_Force_noon{1,1};
 
%% Concatenate the arrays
Sigma_Force = cat(1, Sigma_Force_morn, Sigma_Force_noon);

%% Find the minimums & maximums
force_YLims = zeros(1,2);

% Maximum
force_YLims(1) = prctile(Sigma_Force, max_perc);
% Minimum
force_YLims(2) = prctile(Sigma_Force, min_perc);

if strcmp(norm_force, 'Convert')
    % Maximum
    force_YLims(1) = force_YLims(1) / 1000*gain_factor; % Millivolt conversion * gain
    % Minimum
    force_YLims(2) = force_YLims(2) / 1000*gain_factor; % Millivolt conversion * gain
end

% Expand the axes
force_YLims(1) = force_YLims(1) + axis_expansion;
force_YLims(2) = force_YLims(2) - axis_expansion;










