% O-Information analysis of fMEG
clearvars
close all
rng(329014901)

% try
%     parpool
% catch
%     fprintf('Seems like parallel pool is already open\n')
% end

if exist('OCTAVE_VERSION', 'builtin') && ~exist('fcdf', 'file')
  pkg('load', 'statistics');
end


% %Load data
load("\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\Consciousness_24_ICA\C001-500A-Q-40SSW_Consciousness_20160728_02_FAIRY\3_fMEG\initial_cluster.mat")
data = initial_cluster;
[numChannels,numSamples] = size(data);

%fs = 610.3516; % Sampling rate (Hz)
%t = (0:numSamples-1) / fs; % Create time vector in seconds

assert(size(data,1)==numChannels,'Wrong number of channels')
SIGMA = cov(zscore(data'));
JH = -log(det(SIGMA))/2; % EQ. 3 in Varley et al. 2024
Omega = soinfo_from_covmat(SIGMA,length(data));


figure;
imagesc(SIGMA); 
colorbar;
title('Channel Covariance (Z-scored)');
xlabel('Channel ID'); 
ylabel('Channel ID');
axis square;

windowSize = 2 * 610; % 2-second windows
stepSize = 1 * 610;   % 1-second overlap
nWindows = floor((numSamples - windowSize) / stepSize);

omega_time = zeros(nWindows, 1);

for i = 1:nWindows
    idx = (1:windowSize) + (i-1)*stepSize;
    window_data = data(:, idx);
    window_sigma = cov(zscore(window_data'));
    
    % Use your soinfo function here
    omega_time(i) = soinfo_from_covmat(window_sigma, windowSize);
end

figure;
plot(omega_time);
hold on; yline(0, 'r--'); % Line separating Synergy from Redundancy
title('O-Information Dynamics (\Omega)');
xlabel('Window Index'); ylabel('\Omega (bits)');

% % Extract the variable from the loaded .mat file
% 
% Basic plot of all 10 channels
% plot(data(1,:));
% title('10-Channel MEG Data');
% xlabel('Samples');
% ylabel('fT');