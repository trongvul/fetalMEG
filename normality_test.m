path = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\Consciousness_24_ICA\C001-500A-Q-40SSW_Consciousness_20160728_02_FAIRY\3_fMEG\initial_cluster.mat";
dataStruct = load(path);
data = dataStruct.initial_cluster;

fs = 610.3516;
% Filter Design (1-15 Hz)
% We use a 4th order Butterworth filter
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

% Apply Filter to each channel
% 'filtfilt' requires data to be along the first dimension, 
% so we transpose to [Samples x Channels], filter, then transpose back.
data_filtered = filtfilt(b, a, data')';

assert(size(data_filtered,1)==numChannels,'Wrong number of channels')
[normData, ~] = data2gaussian(data_filtered');
normData = normData';
%SIGMA = cov(data_filtered');
assert(size(normData,1)==numChannels,'Wrong number of channels')
SIGMA = cov(normData);
omegaVal = soinfo_from_covmat(SIGMA, length(normData));

%Select one channel to test (e.g., Channel 1)
channelData = normData(1, :);

figure;
qqplot(channelData);
title('Q-Q Plot of Filtered EEG (Channel 1)');
grid on;


% Test each of the 10 channels
% h = zeros(1, 10);
% p = zeros(1, 10);
% 
% for ch = 1:10
%     [h(ch), p(ch)] = jbtest(data_filtered(ch, :));
% end
% 
% fprintf('Number of channels failing normality (p < 0.05): %d/10\n', sum(h));