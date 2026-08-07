clearvars
close all
% Analysis of o-information of auditory stimulus response
% Uses the cluster-validated 10-channel signal (<subj>_selected10.mat), which
% only exists for ACCEPTED subjects -- so this loop analyses accepted subjects only.

baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516;
postSamples = round(1.5 * fs);
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
% Use a cell array to accumulate results
compiledData = {};
nSkipped = 0;   % accepted-file missing (rejected / not processed)
fprintf('Starting analysis (accepted subjects only)...\n');
for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "AUDIO", 'IgnoreCase', true)

            % Cluster-validated 10-channel file (accepted subjects only)
            selFile = fullfile(baseDir, subjName, '3_fMEG', [subjName '_selected10.mat']);
            if ~isfile(selFile)
                % No selected-10 file -> subject was rejected or not processed -> skip
                nSkipped = nSkipped + 1;
                fprintf('Skipped %s: no selected10 file (not accepted).\n', subjName);
                continue
            end

            % Load the selected 10-channel signal, triggers, and SSW label
            S   = load(selFile, 'sel');
            sel = S.sel;
            sig = sel.signal;                 % 10 x nSamples, RAW selected channels
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat'));
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat'));

            isL1 = strcmp({alltriggers.type}, 'DEV');
            l1_samples = [alltriggers(isL1).sample];
            trialOmegas = nan(length(l1_samples), 1);

            % Filter the data
            data_filt = filtfilt(b, a, sig')';
            nSamplesSig = size(sig, 2);

            for t = 1:length(l1_samples)
                s_idx = l1_samples(t);
                if (s_idx + postSamples) <= nSamplesSig
                    segment = data_filt(:, s_idx : s_idx + postSamples);
                    [~, SIGMA] = data2gaussian(segment');
                    [oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(segment, 2));
                    trialOmegas(t) = oinfo ;
                end
            end

            avgOmega = mean(trialOmegas, 'omitnan');
            % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = avgOmega;
            fprintf('Processed %s: SSW = %s, Omega = %.2f\n', subjName, ssw, avgOmega);
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
save('AvgOmega_SSW_IUGR_AUD_DEV_new.mat', 'results');
fprintf('Analysis saved to AvgOmega_SSW_IUGR_AUD_DEV_new.mat\n');