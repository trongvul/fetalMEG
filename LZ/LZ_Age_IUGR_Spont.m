% O-Information analysis of fMEG with respect to conceptional age
clearvars
close all
rng(329014901)

%if isempty(gcp('nocreate')), parpool; end % Open parallel pool if not already open

% Setup paths
baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
numChannels = 10;
fs = 610.3516;

% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;

    if startsWith(subjName, '.'), continue; end
    try
        if allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true) && contains(subjName, "Spont", 'IgnoreCase', true)
            % Define Paths
            subjPath = fullfile(allSubjectsDir(i).folder, subjName);
            clusterPath = fullfile(subjPath, '3_fMEG', 'initial_cluster.mat');
            sswPath = fullfile(subjPath, '0_info', 'ssw.mat');

            if isfile(clusterPath) && isfile(sswPath)
                % Load MEG data
                dataStruct = load(clusterPath);
                data = dataStruct.initial_cluster;

                fs = 610.3516;
                % Filter Design (1-15 Hz)
                % We use a 4th order Butterworth filter
                [b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

                % Apply Filter to each channel
                % 'filtfilt' requires data to be along the first dimension, 
                % so we transpose to [Samples x Channels], filter, then transpose back.
                data_filtered = filtfilt(b, a, data')';
             
                % Calculate Lempel-Ziv entropy for each channel
                C = zeros(numChannels, 1);
                parfor ch = 1:numChannels
                    ch_data = data_filtered(ch,:);
                    C(ch) = LZ76(ch_data >= median(ch_data));
                end
                
                avg_LZ = mean(C);
                
                % Normalize LZ
                LZ  = avg_LZ * log2(length(data_filtered)) / length(data_filtered);
                
                % Load conceptional age (SSW)
                sswStruct = load(sswPath);
                currentSSW = sswStruct.ssw;

                % Store in cell array: {SubjectName, SSW, Omega}
                compiledData{end+1, 1} = subjName; %#ok<AGROW>
                compiledData{end, 2} = currentSSW;
                compiledData{end, 3} = LZ;

                fprintf('Processed %s | SSW: %d | Lempel-Ziv entropy: %.2f\n', ...
                            subjName, currentSSW, LZ);
            end
        end
    catch ME
        fprintf('Error at index %d (%s): %s\n', i, subjName, ME.message);    
    end
end

% Convert to Table at the very end
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'LZ'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Found %d valid subjects.\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
end

% Save the result
save('LZ_SSW_IUGR_SPONT.mat', 'results');
fprintf('Analysis saved to LZ_SSW_IUGR_SPONT.mat\n');

