% Analysis of Averaged Stimulus Response with Quality Control Plot
baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
preSamples  = round(0.1 * fs); 
postSamples = round(0.6 * fs); 
timeVector  = (-preSamples:postSamples) / fs; % Time axis in seconds

% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');


for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "LIGHT", 'IgnoreCase', true) 

            % 1. Load Data and Triggers
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 

            
            isL1 = strcmp({alltriggers.type}, 'L1');
            l1_samples = [alltriggers(isL1).sample];
            
            % 2. Segment and Baseline Correct
            epochLen = preSamples + postSamples + 1;
            trialMatrix = nan(10, epochLen, length(l1_samples));
            
            for t = 1:length(l1_samples)
                s_idx = l1_samples(t);
                if (s_idx - preSamples) > 0 && (s_idx + postSamples) <= size(initial_cluster, 2)
                    segment = initial_cluster(:, s_idx - preSamples : s_idx + postSamples);
                    baseline_val = mean(segment(:, 1:preSamples), 2);
                    trialMatrix(:,:,t) = segment - baseline_val;
                end
            end
            
            % 3. Average across trials (Evoked Response)
            evoked_response = mean(trialMatrix, 3, 'omitnan'); 
            
            % % 4. Quality Control Plot (Visual Check)
            % figure('Name', ['Subject: ' subjName], 'Color', 'w');
            % plot(timeVector, evoked_response, 'LineWidth', 1.5);
            % hold on;
            % % Draw a patch for the 500ms light flash duration
            % patch([0 0.5 0.5 0], [min(evoked_response(:)) min(evoked_response(:)) ...
            %       max(evoked_response(:)) max(evoked_response(:))], ...
            %       [0.9 0.9 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % grid on;
            % xline(0, '--k', 'Flash Start');
            % xlabel('Time (s)');
            % ylabel('Amplitude (fT)');
            % title(['Averaged Visual Response (n=' num2str(length(l1_samples)) ' trials)']);
            % legend({'Channels 1-10'}, 'Location', 'northeastoutside');
            
            % 5. Filter and Compute O-Info
            [b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
            evoked_filt = filtfilt(b, a, evoked_response')';

            [~, SIGMA] = data2gaussian(evoked_filt');
            omegaVal = soinfo_from_covmat(SIGMA, size(evoked_filt, 2));

             % Store in cell array: {SubjectName, SSW, Omega}
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = omegaVal;

            fprintf('Processed %s: SSW = %d, Omega = %.2f\n', subjName, ssw, omegaVal);        

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
save('Fetal_Omega_SSW_IUGR_LIGHT_AVG.mat', 'results');
fprintf('Analysis saved to Fetal_Omega_SSW_IUGR_LIGHT_AVG.mat\n');