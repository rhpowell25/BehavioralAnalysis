function EMG_Summary(xds, muscle_groups, time_length, Save_File)

%% Basic Settings, some variable extractions, & definitions

% Font specifications
title_font_size = 15;
label_font_size = 10;
figure_width = 700;
figure_height = 750;

%% Find the EMG index

[M] = EMG_Index(xds, muscle_groups);

%% Pull the timeframe, names, & raw EMG of the selected muscles

raw_time_idx = find(round(xds.raw_EMG_time_frame, 4) == time_length);
raw_EMG_timeframe = xds.raw_EMG_time_frame;
EMG_names = xds.EMG_names(M);
raw_EMG = xds.raw_EMG(:,M);
for ii = 1:width(raw_EMG)
    raw_EMG(:,ii) = xds.raw_EMG(:,M(ii));
end

% Define the sampling frequency (Hz)
samp_freq = 1 / (raw_EMG_timeframe(end) / length(raw_EMG_timeframe));

%% Calculate the Welch Power Spectral Density of the raw EMG
% How many indexes do you want each window to be?
window_length = 100000;

[raw_EMG_welch_PSD, welch_freq] = pwelch(raw_EMG, window_length, [], [], samp_freq);
%[raw_EMG_welch_PSD, welch_freq] = pwelch(raw_EMG, [], [], [], samp_freq);

%% Run a Notch filter to remove 60 Hz noise

% Design the notch filer
notch_filter = designfilt('bandstopiir', 'FilterOrder', 4, ...
           'HalfPowerFrequency1', 59,'HalfPowerFrequency2', 61, ...
           'DesignMethod','butter','SampleRate', samp_freq);

% Define the filtered EMG
notched_EMG = zeros(length(raw_EMG), width(raw_EMG));

for ii = 1:width(raw_EMG)
    % Apply the Notch Filter
    notched_EMG(:,ii) = filtfilt(notch_filter, raw_EMG(:, ii));
end

%% Calculate the Welch Power Spectral Density of the notched EMG
%[notched_EMG_welch_PSD, notched_welch_freq] = pwelch(notched_EMG, 100000, [], [], samp_freq);
%[raw_EMG_welch_PSD, welch_freq] = pwelch(raw_EMG, [], [], [], samp_freq);

%% Reject the high amplitude artifacts

% Reject EMG that surpass more than 8 standard deviations
std_limit = 8;

rejected_EMG = notched_EMG;

mean_raw_EMG = zeros(width(notched_EMG),1);
std_raw_EMG = zeros(width(notched_EMG),1);
for ii = 1:width(notched_EMG)
    std_raw_EMG(ii,1) = std(notched_EMG(:,ii));
    mean_raw_EMG(ii,1) = mean(notched_EMG(:,ii));
    high_raw_amps = abs(notched_EMG(:,ii)) > std_limit*std_raw_EMG(ii,1);
    num_reject = length(find(high_raw_amps == 1));
    fprintf("%0.1f Samples Rejected \n", round(num_reject));
    rejected_EMG(high_raw_amps,ii) = (mean_raw_EMG(ii,1) - std_raw_EMG(ii,1)) + ... 
        (mean_raw_EMG(ii,1) + std_raw_EMG(ii,1)) .* rand(1,1);
end

%% High pass filter, rectify, and low pass filter the EMG

% Construct filter off 1/2 the sampling frequency (to prevent aliasing)
nyquist_num = 2;

% High pass 4th order Butterworth band pass filter (50 Hz)
[b_high, a_high] = butter(4, nyquist_num*50/samp_freq, 'high');
highpassed_EMG = filtfilt(b_high, a_high, rejected_EMG);

% Full wave rectification
rect_EMG = abs(highpassed_EMG);

% Low pass 4th order Butterworth band pass filter (10 Hz)
[b_low, a_low] = butter(4, nyquist_num*10/samp_freq, 'low');
lowpassed_EMG = filtfilt(b_low, a_low, rect_EMG);

%% Plot the raw EMG & welch filter

for ii = 1:width(raw_EMG)

    EMG_sum_fig = figure;
    EMG_sum_fig.Position = [200 50 figure_width figure_height];
    Fig_Title = strcat('EMG Summary:', strrep(string(EMG_names(ii)),'EMG_',' '));
    sgtitle(Fig_Title, 'FontSize', (title_font_size + 5), 'Interpreter', 'None');

    % Raw EMG
    subplot(4,1,1)
    hold on
    title('Raw EMG:', 'FontSize', title_font_size)
    xlabel('Time (sec)', 'FontSize', label_font_size)
    ylabel('Amplitude', 'FontSize', label_font_size)
    plot(raw_EMG_timeframe(1:raw_time_idx), raw_EMG(1:raw_time_idx, ii), 'Color', 'k')
    x_max = max(raw_EMG_timeframe(1:raw_time_idx));
    x_min = min(raw_EMG_timeframe(1:raw_time_idx));
    xlim([x_min x_max]);

    % Log Scale Welch Power Spectral Density
    subplot(4,1,2)
    hold on
    title('Log Scale Welch PSD', 'FontSize', title_font_size)
    xlabel('Frequency (Hz)', 'FontSize', label_font_size)
    ylabel('Power / Frequency (dB/Hz)', 'FontSize', label_font_size)
    plot(welch_freq, 10*log10(raw_EMG_welch_PSD(:,ii)), 'Color', 'k')
    figure_axes = gca;
    figure_axes.XScale = 'log';
    x_max = max(welch_freq);
    x_min = min(welch_freq);
    xlim([x_min x_max]);

    % Notched Raw EMG
    subplot(4,1,3)
    hold on
    title('Notched & High-Pass Filtered:', 'FontSize', title_font_size)
    xlabel('Time (sec)', 'FontSize', label_font_size)
    ylabel('Amplitude', 'FontSize', label_font_size)
    plot(raw_EMG_timeframe(1:raw_time_idx), highpassed_EMG(1:raw_time_idx, ii), 'Color', 'k')
    x_max = max(raw_EMG_timeframe(1:raw_time_idx));
    x_min = min(raw_EMG_timeframe(1:raw_time_idx));
    xlim([x_min x_max]);

    % Rectified & Low-Pass Filtered EMG
    subplot(4,1,4)
    hold on
    title('Rectified & Low-Pass Filtered', 'FontSize', title_font_size)
    xlabel('Time (sec)', 'FontSize', label_font_size)
    ylabel('Amplitude', 'FontSize', label_font_size)
    plot(raw_EMG_timeframe(1:raw_time_idx), lowpassed_EMG(1:raw_time_idx,ii), 'Color', 'k')
    x_max = max(raw_EMG_timeframe(1:raw_time_idx));
    x_min = min(raw_EMG_timeframe(1:raw_time_idx));
    xlim([x_min x_max]);

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end





