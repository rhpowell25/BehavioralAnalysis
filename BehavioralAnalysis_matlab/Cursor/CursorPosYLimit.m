
function cursor_YLims = CursorPosYLimit(xds_morn, xds_noon)

%% Ending the function if there is no cursor position

if strcmp(xds_morn.meta.task, 'PG')
    disp('No Cursor Information In Powergrasp');
    cursor_YLims = NaN;
    return
end

%% Display the functions being used
disp('Cursor Y-Limit Function:');

%% Basic Settings, some variable extractions, & definitions
% Percentiles for maximums & minimums
max_perc = 100;
min_perc = 5;

axis_expansion = 0.25;
    
%% Extracting cursor position

cursor_morn = xds_morn.curs_p;
cursor_noon = xds_noon.curs_p;
    
%% Find the vector sum of cursor position

z_cursor_morn = sqrt(cursor_morn(:, 2).^2 + cursor_morn(:, 1).^2);
z_cursor_noon = sqrt(cursor_noon(:, 2).^2 + cursor_noon(:, 1).^2);
 
%% Concatenate the arrays
z_cursor_p = cat(1, z_cursor_morn, z_cursor_noon);

%% Find the minimums & maximums
cursor_YLims = zeros(1,2);

% Maximum
cursor_YLims(1) = prctile(z_cursor_p, max_perc) + axis_expansion;
% Minimum
cursor_YLims(2) = prctile(z_cursor_p, min_perc) - axis_expansion;


