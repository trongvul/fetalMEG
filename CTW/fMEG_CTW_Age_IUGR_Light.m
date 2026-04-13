clearvars;
close all;
rng(1234);

% Open parallel pool if not already open
if isempty(gcp('nocreate')), parpool; end

% Analysis of Averaged Stimulus Response with Quality Control Plot
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
preSamples  = round(0.1 * fs); 
postSamples = round(1.4 * fs); 
timeVector  = (-preSamples:postSamples) / fs; % Time axis in seconds

vmo = 20; % Context depth 
numSymbols = 2; % Binary complexity

% Design filter
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

            isL1 = strcmp({alltriggers.type}, 'L1');
            l1_samples = [alltriggers(isL1).sample];
            numTrials = length(l1_samples);
            
            % Segment and Baseline Correct
            epochLen = preSamples + postSamples + 1;
            trialMatrix = nan(10, epochLen, numTrials);
            
            % Filter the raw data
            data_filt = filtfilt(b, a, initial_cluster')';

            % Init matrix for CTW in each trial
            trialEntropies = nan(numTrials, 1);
            
            parfor t = 1:numTrials
                s_idx = l1_samples(t);
                if (s_idx - preSamples) > 0 && (s_idx + postSamples) <= size(initial_cluster, 2)
                    segment = data_filt(:, s_idx - preSamples : s_idx + postSamples);
                    baseline_val = mean(segment(:, 1:preSamples), 2);
                    segment = segment - baseline_val;
                    trialEntropies(t) = CTWEntropyRate(segment, numSymbols, vmo);
                end
            end
            
            % Compute avarage CTW across all trials
            ctwVal = mean(trialEntropies, 'omitnan');

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = ctwVal;

            fprintf('Processed %s: SSW = %d, CTW = %.2f\n', subjName, ssw, ctwVal);        

        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'CTW'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('CTW_SSW_IUGR_Light.mat', 'results');
fprintf('Analysis saved to CTW_SSW_IUGR_Light.mat\n');