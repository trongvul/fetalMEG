clearvars
close all

% Analysis of o-information of auditory stimulus response 
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 

[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    %condtion = nan;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) 
            if contains(subjName, "AUDIO", 'IgnoreCase', true)
                condition = "AUDIO";
            elseif contains(subjName, "LIGHT", 'IgnoreCase', true)
                condition = "LIGHT";
            elseif contains(subjName, "SPONT", 'IgnoreCase', true)
                condition = "SPONT";
            end

            % Load Data and Triggers
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 

            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';

            [coeff, score, latent, tsquared, explained] = pca(data_filt');

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = condition;
            compiledData{end, 3} = ssw;
            compiledData{end, 4} = explained(1);
            compiledData{end, 5} = explained(2);
            compiledData{end, 6} = explained(3);

            fprintf('Processed %s: Condition = %s, SSW =%s, Explained1 = %.2f, Explained2 = %.2f, Explained3 = %.2f\n', subjName, condition, ssw, explained(1), explained(2), explained(3));        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'Condition', 'SSW', 'Explained1', 'Explained2', 'Explained3'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('PCA_results.mat', 'results');
fprintf('Analysis saved to PCA_results.mat\n');