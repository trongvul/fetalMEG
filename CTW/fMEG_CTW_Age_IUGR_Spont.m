clearvars
close all

% Open parallel pool if not already open
if isempty(gcp('nocreate')), parpool; end

% Analysis of Spontaneous (Resting State) Complexity (Average of Trial Entropies)
baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;

trialDuration = 1.5; % seconds
trialSamples  = round(trialDuration * fs); 
vmo = 10; 
numSymbols = 2; 

% Design filter
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

compiledData = {}; 
fprintf('Starting spontaneous Entropy-Average analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SPONT", 'IgnoreCase', true) 
            
            % Load Data
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 
            
            % Slice into trials
            totalSamples = size(initial_cluster, 2);
            startIndices = 1 : trialSamples : (totalSamples - trialSamples);
            numTrials = length(startIndices);
            
            % Filter entire continuous data first (more efficient)
            data_filt = filtfilt(b, a, initial_cluster')';
            
            % % Array to store entropy of each trial
            trialEntropies = nan(numTrials, 1);
            
            % Loop through trials and compute CTW
            parfor t = 1:numTrials
                s_idx = startIndices(t);
                segment = data_filt(:, s_idx : s_idx + trialSamples - 1);
                trialEntropies(t) = CTWEntropyRate(segment, numSymbols, vmo);
            end

            % Compute avarage CTW across all trials
            ctwVal = mean(trialEntropies, 'omitnan');
            
            % Store Results
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = ctwVal;
            
            fprintf('Processed %s: SSW = %d, Mean CTW = %.4f\n', ...
                subjName, ssw, ctwVal);        
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
save('CTW_SSW_IUGR_SPONT.mat', 'results');
fprintf('Analysis saved to CTW_SSW_IUGR_SPONT.mat\n');