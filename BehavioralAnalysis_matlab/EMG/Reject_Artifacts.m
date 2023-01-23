function rejected_EMG = Reject_Artifacts(xds_EMG, std_limit)

%% Reject high amplitude EMG artifacts
rejected_EMG = xds_EMG;

% What standard dev do you want to replace the rejected EMG with?
%reject_replacement = 1;

for jj = 1:width(xds_EMG)
    std_EMG = std(xds_EMG(:,jj));
    high_amps = abs(xds_EMG(:,jj)) > std_limit*std_EMG;
    mean_EMG = mean(xds_EMG(:,jj));
    rejected_EMG(high_amps,jj) = (mean_EMG - std_EMG) + (mean_EMG + std_EMG) .* rand(1,1);
end

