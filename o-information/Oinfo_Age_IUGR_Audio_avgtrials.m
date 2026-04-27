clearvars
close all

% Analysis of o-information of auditory stimulus response 
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
postSamples = round(1.5 * fs);
epochLen = postSamples + 1;

[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "AUDIO", 'IgnoreCase', true) 

            % Load Data and Triggers
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 
            
            isStim = strcmp({alltriggers.type}, 'STD');
            stimSamples = [alltriggers(isStim).sample];
            numTrials = length(stimSamples);

            %trialOmegas = nan(length(l1_samples), 1);
            
            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';

            trialMatrix = nan(10, epochLen, numTrials);

            for t = 1:numTrials
                s_idx = stimSamples(t);
                if  (s_idx + postSamples) <= size(initial_cluster, 2)
                    segment = data_filt(:, s_idx : s_idx + postSamples);
                    trialMatrix(:,:,t) = segment;
                    %[~, SIGMA] = data2gaussian(segment');
                    %[oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(segment, 2));
                    %trialOmegas(t) = oinfo / sinfo;
                end
            end
            
            % Average LZ over trials
            evoked_response = mean(trialMatrix, 3, 'omitnan'); 

            [~, SIGMA] = data2gaussian(evoked_response');
            [oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(segment, 2));
            omegaVal = oinfo / sinfo;

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = omegaVal;

            fprintf('Processed %s: SSW = %s, Omega = %.2f\n', subjName, ssw, omegaVal);        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'Omega'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('AvgOmega_SSW_IUGR_AUD_STD_norm_avgtrials.mat', 'results');
fprintf('Analysis saved to AvgOmega_SSW_IUGR_AUD_STD_norm_avgtrials.mat\n');