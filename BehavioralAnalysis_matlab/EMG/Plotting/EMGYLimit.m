function max_EMGYLim = EMGYLimit(xds_morn, xds_noon, event, muscle_group, EMG_Zero_Factor, EMG_Norm_Factor)

%% Display the functions being used
disp('EMG Y-Limit Function:');

%% Extract the target directions & centers
[target_dirs_morn, target_centers_morn] = Identify_Targets(xds_morn);
[target_dirs_noon, target_centers_noon] = Identify_Targets(xds_noon);

%% Begin the loop through all directions
avg_EMG_morn = struct([]);
for jj = 1:length(target_dirs_morn)
    [avg_EMG_morn{jj,1}, ~ ]= ...
        EventWindowEMG(xds_morn, muscle_group, target_dirs_morn(jj), target_centers_morn(jj), event);
end

avg_EMG_noon = struct([]);
for jj = 1:length(target_dirs_noon)
    [avg_EMG_noon{jj,1}, ~] = ...
        EventWindowEMG(xds_noon, muscle_group, target_dirs_noon(jj), target_centers_noon(jj), event);
end

%% Zero the average EMG
for jj = 1:length(avg_EMG_morn)
    for ii = 1:height(avg_EMG_morn{1,1})
        avg_EMG_morn{jj,1}(ii,:) = avg_EMG_morn{jj,1}(ii,:) - EMG_Zero_Factor(ii);
    end
end
for jj = 1:length(avg_EMG_noon)
    for ii = 1:height(avg_EMG_noon{1,1})
        avg_EMG_noon{jj,1}(ii,:) = avg_EMG_noon{jj,1}(ii,:) - EMG_Zero_Factor(ii);
    end
end
        
%% Normalize the average EMG
for jj = 1:length(avg_EMG_morn)
    for ii = 1:height(avg_EMG_morn{1,1})
        avg_EMG_morn{jj,1}(ii,:) = (avg_EMG_morn{jj,1}(ii,:) / EMG_Norm_Factor(ii))*100;
    end
end
for jj = 1:length(avg_EMG_noon)
    for ii = 1:height(avg_EMG_noon{1,1})
        avg_EMG_noon{jj,1}(ii,:) = (avg_EMG_noon{jj,1}(ii,:) / EMG_Norm_Factor(ii))*100;
    end
end

%% Finding the maximum EMG amplitude
max_morn_amp = zeros(length(avg_EMG_morn), height(avg_EMG_morn{1,1}));
max_noon_amp = zeros(length(avg_EMG_noon), height(avg_EMG_noon{1,1}));

for ii = 1:length(avg_EMG_morn)
    for mm = 1:height(avg_EMG_morn{1,1})
        max_morn_amp(ii,mm) = max(avg_EMG_morn{ii}(mm,:));
    end
end

for ii = 1:length(avg_EMG_noon)
    for mm = 1:height(avg_EMG_morn{1,1})
        max_noon_amp(ii,mm) = max(avg_EMG_noon{ii}(mm,:));
    end
end

%% Concatenate the morning and afternoon maximums
max_amp = cat(1, max_morn_amp, max_noon_amp);

%% Find the maximum amplitude of each EMG channel
max_EMGYLim = zeros(height(avg_EMG_morn{1,1}),1);
for ii = 1:length(max_EMGYLim)
    max_EMGYLim(ii) = max(max_amp(:,ii));
end
% Round up to the nearest fifth digit
for ii = 1:length(max_EMGYLim)
    max_EMGYLim(ii) = (round(max_EMGYLim(ii) / 5)) * 5 + 10;
end


