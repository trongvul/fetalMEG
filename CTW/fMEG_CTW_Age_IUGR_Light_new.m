clearvars;
close all;
rng(1234);

% Open parallel pool if not already open
%if isempty(gcp('nocreate')), parpool; end

baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
%baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
preSamples  = round(0.1 * fs); 
postSamples = round(1.4 * fs); 
epochLen    = preSamples + postSamples + 1;

vmo = 20; % Context depth 
numSymbols = 2; % Binary complexity

% Design filter
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

compiledData = {}; 
fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "LIGHT", 'IgnoreCase', true) 
            
            % Load Data
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 
            
            isL1 = strcmp({alltriggers.type}, 'L1');
            l1_samples = [alltriggers(isL1).sample];
            numTrials = length(l1_samples);
            
            % Filter the raw data once (entire continuous signal)
            data_filt = filtfilt(b, a, initial_cluster')';
            
            % Pre-allocate the 3D matrix [Channels x Time x Trials]
            trialMatrix = nan(size(data_filt, 1), epochLen, numTrials);
            
            validTrialCount = 0;
            for t = 1:numTrials
                s_idx = l1_samples(t);
                % Check boundaries
                if (s_idx - preSamples) > 0 && (s_idx + postSamples) <= size(data_filt, 2)
                    segment = data_filt(:, s_idx - preSamples : s_idx + postSamples);
                    
                    % Baseline Correction
                    baseline_val = mean(segment(:, 1:preSamples), 2);
                    segment = segment - baseline_val;
                    
                    % Store in 3D Matrix
                    trialMatrix(:, :, t) = segment;
                    validTrialCount = validTrialCount + 1;
                end
            end
            
            % Remove any trials that were NaN due to boundary issues
            trialMatrix = trialMatrix(:, :, 1:validTrialCount);
            
            % Pass 3D Tensor to CTWEntropyRate
            if validTrialCount > 0
                ctwVal = CTWEntropyRate(trialMatrix, numSymbols, vmo);
                
                % Store results
                compiledData{end+1, 1} = subjName; %#ok<AGROW>
                compiledData{end, 2} = ssw;
                compiledData{end, 3} = ctwVal;
                fprintf('Processed %s: SSW = %d, CTW = %.2f (from %d trials)\n', ...
                    subjName, ssw, ctwVal, validTrialCount);
            end
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Final Table formatting
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'CTW'});
    results = rmmissing(results);
    save('CTW_SSW_IUGR_Light_new.mat', 'results');
    fprintf('Analysis saved. Total valid subjects: %d\n', height(results));
end