clearvars;
close all;
rng(1234);

if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Analysis of spontaneous (resting state) complexity
baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
%baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;
m = 3;
t = 32; %s
tau = round((t/1000)/(1/fs)); % convert from ms to sample

trialDuration = 1.5; % seconds
trialSamples  = round(trialDuration * fs); 

[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

compiledData = {}; 
fprintf('Starting spontaneous Entropy-Average analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir  && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "SPONT", 'IgnoreCase', true) 
            
            % Load Data
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 
            
            % Slice into trials
            numChannels = size(initial_cluster, 1);
            totalSamples = size(initial_cluster, 2);
            startIndices = 1 : trialSamples : (totalSamples - trialSamples);
            numTrials = length(startIndices);
            
            % Filter entire continuous data first (more efficient)
            data_filt = filtfilt(b, a, initial_cluster')';

            % Array to store entropy of each channel
            channelPerms = nan(numChannels, 1);

            % Array to store entropy of each trial within a channel
            trialPerms = nan(numTrials, 1);

            for ch = 1:numChannels
                parfor t = 1:numTrials
                    s_idx = startIndices(t);
                    segment = data_filt(ch, s_idx : s_idx + trialSamples - 1);
                    trialPerms(t) = PermEn(segment, m, tau);
                end
                % Average PermEn across trials
                channelPerms(ch) = mean(trialPerms, 'omitnan');
            end
            
            % Average PermEn over channels
            permEn = mean(channelPerms, 'omitnan');
            
            % Store Results
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = permEn;
            
            fprintf('Processed %s: SSW = %d, Mean PermEn = %.4f\n', ...
                subjName, ssw, permEn);        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'PermEn'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('PermEn_SSW_IUGR_SPONT.mat', 'results');
fprintf('Analysis saved to PermEn_SSW_IUGR_SPONT.mat\n');