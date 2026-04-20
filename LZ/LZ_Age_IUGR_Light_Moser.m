clearvars
close all

%if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Analysis of LZ of visual auditory stimulus response 
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
postSamples = round(0.7 * fs); 
epochLen = postSamples + 1;

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
            numTrials = length(l1_samples);
                        
            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';

            rms_signal = rms(data_filt, 1);

            trialMatrix = nan(numTrials, epochLen);
            
            for t = 1:numTrials
                s_idx = l1_samples(t);
                if (s_idx + postSamples) <= size(initial_cluster, 2)
                    segment = rms_signal(s_idx : s_idx + postSamples);
                    trialMatrix(t,:) = segment;
                end
            end
             
            % Average LZ over trials
            evoked_response = mean(trialMatrix, 1, 'omitnan'); 
            
            LZ = LZ76(evoked_response >= median(evoked_response));
            % Normalized LZ
            normLZ = LZ * log2(length(evoked_response)) / length(evoked_response);
            
            % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = normLZ;

            fprintf('Processed %s: SSW = %s, LZ = %.2f\n', subjName, ssw, normLZ);        
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
save('LZ_SSW_IUGR_LIGHT_Moser_07.mat', 'results');
fprintf('Analysis saved to LZ_SSW_IUGR_LIGHT_Moser_07.mat\n');