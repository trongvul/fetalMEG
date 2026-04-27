clearvars
close all

%if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Analysis of spontaneous (resting state) complexity (average of trial entropies)
%baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;

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
          
            % Filter entire continuous data first (more efficient)
            data_filt = filtfilt(b, a, initial_cluster')';
           
            [~, SIGMA] = data2gaussian(data_filt');
            [oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(data_filt, 2));
            omegaVal = oinfo / sinfo;
          
            % Store Results
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = omegaVal;
            
            fprintf('Processed %s: SSW = %s, Mean Omega = %.4f\n', ...
                subjName, ssw, omegaVal);        
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
save('AvgOmega_SSW_IUGR_SPONT_notrials_norm.mat', 'results');
fprintf('Analysis saved to AvgOmega_SSW_IUGR_SPONT_notrials_norm.mat\n');