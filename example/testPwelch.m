%% Test Pwelch CLASS CPP
clc; clear;

% load data
data = readmatrix('rawdata.csv');
% parameters pwelch with Hamming
wlength = 0.5; % we need it in seconds, so 256/512
pshift = 0.25; % overlapping 128 / sampleRate 512 ---> 0.25
wshift = 0.0625; % new signlas after these seconds
mlength = 1; % moving average biggest windows
SampleRate = 512;
wconv = 'backward';

% apply psd
[matlab_psd, f] = (proc_spectrogram(data, wlength, wshift, pshift, SampleRate));
nfreqs = size(matlab_psd, 2);
nwindow = size(matlab_psd, 1);
nchannels = size(matlab_psd, 3);
matlab_psd = log(matlab_psd);

% load psd rosneuro
psd_rosneuro = readmatrix('psd_window.csv'); % [(windows*freqs) x channels]
% give the same structure of matlab_psd
cont_window = 1;
rosneuro_psd = nan(size(matlab_psd));
for start_window = 1:nfreqs:size(psd_rosneuro,1)
    end_window = start_window + nfreqs - 1;
    rosneuro_psd(cont_window,:,:) = psd_rosneuro(start_window:end_window,:);
    cont_window = cont_window + 1;
end


%% Plot data of a selected channel CLASS CPP
nsegments = size(matlab_psd, 1);

% select channels and frequencies
channelsSelected = [1,2,3];
freqsSelected = [1,2,3, 127, 128, 129];

% iterate over channels
for idx_ch = channelsSelected
    figure;
    t = 1:1:nsegments;
    idx_left = 1;
    idx_right = 2;
    for idx_freq = 1:length(freqsSelected)
        c_freq = freqsSelected(idx_freq);
        subplot(length(freqsSelected), 2, idx_left);
        hold on;
        plot(t, matlab_psd(:,c_freq, idx_ch), '--r', 'LineWidth', 1);
        plot(t, rosneuro_psd(:,c_freq, idx_ch), '-.b', 'LineWidth', 1);
        legend('matlab', 'rosneuro');
        title(['Pwelch | channel=' num2str(idx_ch) ' | freq=' num2str(freqsSelected(idx_freq))]);
        hold off;

        subplot(length(freqsSelected), 2, idx_right);
        grid on;
        plot(t, abs(matlab_psd(:,c_freq, idx_ch) - rosneuro_psd(:,c_freq, idx_ch)))
        title('Difference')

        idx_left = idx_left+2;
        idx_right = idx_right + 2;
    end

    sgtitle('Evaluation Pwelch algorithm using class cpp')
end

% display the biggest difference in the values
disp(['Biggest difference using the class cpp: ' num2str(max(abs(rosneuro_psd - matlab_psd), [], 'all'))])


%% Test pwelch with the ROS NODE
load('psd_nodePwelch.csv'); % to have this file you need to do:
                            % rostopic echo -p /neurodata_psd > ~/psd_nodePwelch.csv
psd_nodePwelch = psd_nodePwelch(:,9:end); % [windows x (channels*freqs)]

% change the structured as before
rosneuroNode_psd = nan(size(matlab_psd));
for idx_window = 1:nwindow
    for idx_freq = 1:nfreqs
        start_ = (idx_freq-1)*nchannels + 1;
        end_ = start_ + nchannels -1;
        rosneuroNode_psd(idx_window, idx_freq, :) = psd_nodePwelch(idx_window, start_:end_);
    end
end

% display the biggest difference in the values
disp(['Biggest difference using the ros node: ' num2str(max(abs(rosneuroNode_psd - matlab_psd), [], 'all'))])

%% Plot data of a selected channel ROS NODE
nsegments = size(matlab_psd, 1);

% select channels and frequencies
channelsSelected = [1,2,3];
freqsSelected = [1,2,3, 127, 128, 129];

% iterate over channels
for idx_ch = channelsSelected
    figure;
    t = 1:1:nsegments;
    idx_left = 1;
    idx_right = 2;
    for idx_freq = 1:length(freqsSelected)
        c_freq = freqsSelected(idx_freq);
        subplot(length(freqsSelected), 2, idx_left);
        hold on;
        plot(t, matlab_psd(:,c_freq, idx_ch), '--r', 'LineWidth', 1);
        plot(t, rosneuroNode_psd(:,c_freq, idx_ch), '-.b', 'LineWidth', 1);
        legend('matlab', 'rosneuro');
        title(['Pwelch | channel=' num2str(idx_ch) ' | freq=' num2str(freqsSelected(idx_freq))]);
        hold off;

        subplot(length(freqsSelected), 2, idx_right);
        grid on;
        plot(t, abs(matlab_psd(:,c_freq, idx_ch) - rosneuroNode_psd(:,c_freq, idx_ch)))
        title('Difference')

        idx_left = idx_left+2;
        idx_right = idx_right + 2;
    end

    sgtitle('Evaluation Pwelch algorithm using ros node pwelch')
end

%% Check if taken same features
chans = [3,3,5,5]; % they are the channels, in this case also the index 
freqs = [8, 10,12,14]; % they are the freqs, so you need to compute the index
all_freqs = 0:2:256;
allFeatures = [];

for idx_window = 1:size(matlab_psd,1)
    features = [];
    for idx = 1:length(chans)
        c_channel = chans(idx);
        c_freq = freqs(idx);
        c_idx_freq = find(all_freqs == c_freq);
        features = cat(1, features, matlab_psd(idx_window,c_idx_freq,c_channel));
    end
    allFeatures = cat(1, allFeatures, features);
end

% load rosneuro features
load('/home/paolo/rosneuro_ws/src/rosneuro_pwelch/test/features_rosneuro.csv');

