clearvars
close all
%if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open
% Analysis of spontaneous (resting state) complexity (average of trial entropies)
% Uses the cluster-validated 10-channel signal (<subj>_selected10.mat), which
% only exists for ACCEPTED subjects -- so this loop analyses accepted subjects only.

baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;
trialDuration = 1.5; % seconds
trialSamples  = round(trialDuration * fs);
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
compiledData = {};
nSkipped = 0;   % accepted-file missing (rejected / not processed)
fprintf('Starting spontaneous Entropy-Average analysis (accepted subjects only)...\n');
for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir  && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "SPONT", 'IgnoreCase', true)

            % Cluster-validated 10-channel file (accepted subjects only)
            selFile = fullfile(baseDir, subjName, '3_fMEG', [subjName '_selected10.mat']);
            if ~isfile(selFile)
                % No selected-10 file -> subject was rejected or not processed -> skip
                nSkipped = nSkipped + 1;
                fprintf('Skipped %s: no selected10 file (not accepted).\n', subjName);
                continue
            end

            % Load the selected 10-channel signal and the SSW label
            S   = load(selFile, 'sel');
            sel = S.sel;
            sig = sel.signal;                 % 10 x nSamples, RAW selected channels
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat'));   % -> ssw

            % Filter entire continuous data first (more efficient)
            data_filt = filtfilt(b, a, sig')';

            [~, SIGMA] = data2gaussian(data_filt');
            [oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(data_filt, 2));
            
            % Store Results
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = oinfo;
            fprintf('Processed %s: SSW = %s, Omega = %d\n', ...
                subjName, ssw, oinfo);
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end
% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'Omega'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid (accepted) subjects. Skipped %d non-accepted.\n', ...
        height(results), nSkipped);
else
    fprintf('\nNo subjects matched the criteria.\n');
end
% Save the result
save('AvgOmega_SSW_IUGR_SPONT_new.mat', 'results');
fprintf('Analysis saved to AvgOmega_SSW_IUGR_SPONT_new.mat\n');