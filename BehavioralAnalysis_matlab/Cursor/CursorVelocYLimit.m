
function cursor_veloc_YLim = CursorVelocYLimit(xds_morn, xds_noon)

%% Ending the function if there is no cursor position

if strcmp(xds_morn.meta.task, 'PG')
    disp('No Cursor Information In Powergrasp');
    cursor_veloc_YLim = NaN;
    return
end

%% Display the functions being used
disp('Cursor Velocity Y-Limit Function:');

%% Basic Settings, some variable extractions, & definitions
% Percentiles for maximums & minimums
max_perc = 100;
min_perc = 5;

axis_expansion = 0;
    
%% Extracting cursor position

cursor_veloc_morn = xds_morn.curs_v;
cursor_veloc_noon = xds_noon.curs_v;
    
%% Find the vector sum of cursor velocity

z_cursor_veloc_morn = sqrt(cursor_veloc_morn(:, 2).^2 + cursor_veloc_morn(:, 1).^2);
z_cursor_veloc_noon = sqrt(cursor_veloc_noon(:, 2).^2 + cursor_veloc_noon(:, 1).^2);
 
%% Concatenate the arrays
z_cursor_v = cat(1, z_cursor_veloc_morn, z_cursor_veloc_noon);

%% Find the minimums & maximums
cursor_veloc_YLim = zeros(1,2);

% Maximum
cursor_veloc_YLim(1) = prctile(z_cursor_v, max_perc) + axis_expansion;
% Minimum
cursor_veloc_YLim(2) = prctile(z_cursor_v, min_perc) - axis_expansion;


