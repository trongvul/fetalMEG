clearvars
close all

if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Analysis of LZ of visual auditory stimulus response 
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
preSamples  = round(0.1 * fs); 
postSamples = round(1.4 * fs); 

[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "LIGHT", 'IgnoreCase', true) 

            % Load Data and Triggers
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 

            numChannels = size(initial_cluster, 1);
            totalSamples = size(initial_cluster, 2);
            
            isL1 = strcmp({alltriggers.type}, 'L1');
            l1_samples = [alltriggers(isL1).sample];
                        
            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';

             % Array to store entropy of each channel
            channelLZs = nan(numChannels, 1);
            
            parfor ch = 1:numChannels
            % Init matrix to store LZ in each trial (epoch)
            trialLZs = nan(totalSamples, 1);
                for t = 1:length(l1_samples)
                    s_idx = l1_samples(t);
                    if (s_idx - preSamples) > 0 && (s_idx + postSamples) <= size(initial_cluster, 2)
                        segment = data_filt(ch, s_idx - preSamples : s_idx + postSamples);
                        baseline_val = mean(segment(:, 1:preSamples), 2);
                        segment = segment - baseline_val;
                        LZ = LZ76(segment >= median(segment));
                        trialLZs(t) = LZ * log2(length(segment)) / length(segment);
                    end
                end
                 % Average LZ over trials
                channelLZs(ch) = mean(trialLZs, 'omitnan');
            end
            
            % Average LZ over channels
            avgLZ = mean(channelLZs, 'omitnan');

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = avgLZ;

            fprintf('Processed %s: SSW = %d, LZ = %.2f\n', subjName, ssw, avgLZ);        
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
save('AvgLZ_SSW_IUGR_LIGHT_15_baseline.mat', 'results');
fprintf('Analysis saved to AvgLZ_SSW_IUGR_LIGHT_15_baseline.mat\n');