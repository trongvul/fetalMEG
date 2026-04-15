clearvars;
close all;
rng(1234);

%if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;

trialDuration = 0.7; % seconds
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

            trialMatrix = nan(numChannels,trialSamples, numTrials);
            
            for t = 1:numTrials
                s_idx = startIndices(t);
                segment = data_filt(:, s_idx : s_idx + trialSamples - 1);
                trialMatrix(:,:,t) = segment;
            end

            % Average LZ over trials
            evoked_response = mean(trialMatrix, 3, 'omitnan'); 


            channelLZs = nan(numChannels, 1);

            for ch = 1:numChannels
                signal = evoked_response(ch,:);
                LZ = LZ76(signal >= median(signal));
                % Normalized LZ
                channelLZs(ch) = LZ * log2(length(signal)) / length(signal);
            end
            
            avgLZ = mean(channelLZs,"omitnan");
            
            % Store Results
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = avgLZ;
            
            fprintf('Processed %s: SSW = %s, Mean LZ = %.4f\n', ...
                subjName, ssw, avgLZ);        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'LZ'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('LZ_SSW_IUGR_SPONT_AvgTrials_07.mat', 'results');
fprintf('Analysis saved to LZ_SSW_IUGR_SPONT_AvgTrials_07.mat\n');