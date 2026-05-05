clearvars
close all

baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 
postSamples = round(1.5 * fs); 
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

compiledData = {}; 
fprintf('Starting PCA-based analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "LIGHT", 'IgnoreCase', true) 
            load(fullfile(baseDir, subjName, '3_fMEG', 'initial_cluster.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'alltriggers.mat')); 
            load(fullfile(baseDir, subjName, '0_info', 'ssw.mat')); 
            
            isL1 = strcmp({alltriggers.type}, 'L1');
            l1_samples = [alltriggers(isL1).sample];
            trialOmegas = nan(length(l1_samples), 1);
            
            % Filter the data
            data_filt = filtfilt(b, a, initial_cluster')';
            
            % Compute PCA on the whole filtered recording
            % 'coeff' are the weights to transform sensors to components
            [coeff, score, latent] = pca(data_filt'); 
            
            % OPTIONAL: Keep only components explaining 95% of variance
            % cumulativeVar = cumsum(latent) / sum(latent);
            % numComponents = find(cumulativeVar >= 0.95, 1);
            % coeff = coeff(:, 1:numComponents);

            for t = 1:length(l1_samples)
                s_idx = l1_samples(t);
                if (s_idx + postSamples) <= size(data_filt, 2)
                    % Extract segment in sensor space
                    segment_sensor = data_filt(:, s_idx : s_idx + postSamples);
                    
                    % 3. Transform segment to PCA space
                    % (Channels x Time) -> (Components x Time)
                    segment_pca = (segment_sensor' * coeff)';
                    
                    % Compute O-info in Component Space
                    [~, SIGMA] = data2gaussian(segment_pca');
                    [oinfo, sinfo] = soinfo_from_covmat(SIGMA, size(segment_pca, 2));
                    trialOmegas(t) = oinfo;
                end
            end
            
            avgOmega = mean(trialOmegas, 'omitnan');
            compiledData{end+1, 1} = subjName; %#ok<AGROW>
            compiledData{end, 2} = ssw;
            compiledData{end, 3} = avgOmega;
            fprintf('Processed %s (PCA): SSW = %s, Omega = %.2f\n', subjName, ssw, avgOmega);        
        end
    catch ME
        fprintf('Failed subject %s: %s\n', subjName, ME.message);
    end
end

if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'Omega'});
    results = rmmissing(results);
end

save('AvgOmega_PCA_Space_LIGHT.mat', 'results');