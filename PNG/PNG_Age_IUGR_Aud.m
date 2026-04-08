clearvars
close all

if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Analysis of PNG of auditory stimulus response 

%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
preSamples  = round(0.1 * fs); 
postSamples = round(1.4 * fs);

tau = 4; % tau tag value (in ms)
sampint = tau/(1000/fs);
nsymb = 64;

[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "AUD", 'IgnoreCase', true) 

            % Load Data and Triggers
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 

            numChannels = size(initial_cluster, 1);
            totalSamples = size(initial_cluster, 2);
            
            isStim = strcmp({alltriggers.type}, 'DEV');
            stimSamples = [alltriggers(isStim).sample];
            numTrials = length(stimSamples);
                        
            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';
            
            % Init matrix to store entropy in each trial (epoch)
            trialPNGs = nan(numTrials, 1);
            for t = 1:numTrials
                s_idx = stimSamples(t);
                if (s_idx - preSamples) > 0 && (s_idx + postSamples) <= size(initial_cluster, 2)
                    segment = data_filt(:, s_idx - preSamples : s_idx + postSamples);
                    baseline_val = mean(segment(:, 1:preSamples), 2);
                    segment = segment - baseline_val;
                    % Compute PNG for a segment in each channel
                    PNG = PNGentropy(segment, sampint, nsymb, 'name');
                    % Average PNG across channels
                    trialPNGs(t) = mean(PNG, 'omitnan'); 
                end
            end
            
            % Avarage PNG across trials
            avgPNG = mean(trialPNGs, 'omitnan');

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = avgPNG;

            fprintf('Processed %s: SSW = %d, PNG = %.2f\n', subjName, ssw, avgPNG);        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'PNG'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('PNG_SSW_IUGR_AUD_DEV.mat', 'results');
fprintf('Analysis saved to PNG_SSW_IUGR_AUD_DEV.mat\n');