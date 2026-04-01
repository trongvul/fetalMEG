% O-Information analysis of Neonatal fMEG with respect to Post-Birth Age
clearvars
close all
rng(329014901)

% Load Age Metadata CSV
ageTable = readtable('NeonatalTable.csv'); 

% Setup paths
baseDir = "\\10.174.1.22\data\FETAL_USER\FAIRY_THEM_ALL_2024_withICA\Consciousness_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
numChannels = 10;
fs = 610.3516;

% Use a cell array to accumulate results
compiledData = {}; 

fprintf('Starting Neonatal analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    if startsWith(subjName, '.'), continue; end
    
    try
        if allSubjectsDir(i).isdir && contains(subjName, "NEO", 'IgnoreCase', true) && contains(subjName, "SPONT", 'IgnoreCase', true)
            
            % 1. Extract ID from folder name
            idTokens = regexp(subjName, 'C0*(\d+)', 'tokens');
            
            if ~isempty(idTokens)
                currentID = str2double(idTokens{1}{1});
                
                % 2. Find matching rows in CSV for this ID
                % Check if column name is ID (case sensitive)
                allMatchesIdx = find(ageTable.ID == currentID);
                
                if ~isempty(allMatchesIdx)
                    % 3. Check if the MEG data file exists
                    subjPath = fullfile(allSubjectsDir(i).folder, subjName);
                    clusterPath = fullfile(subjPath, '3_fMEG', 'initial_cluster.mat');
                    
                    if isfile(clusterPath)
                        % 4. Load and Process Data
                        dataStruct = load(clusterPath);
                        data = dataStruct.initial_cluster;
                        
                        % Filter Design (1-15 Hz)
                        [b, a] = butter(4, [1 15] / (fs/2), 'bandpass');
                        data_filtered = filtfilt(b, a, data')';
                        
                        % Calculate Omega
                        [~, SIGMA] = data2gaussian(data_filtered');
                        omegaVal = soinfo_from_covmat(SIGMA, length(data_filtered));

                        % 5. Get Age from the FIRST match found
                        rowIdx = allMatchesIdx(1); % Corrected logic
                        currentAge = ageTable.age(rowIdx);
                      
                        % 6. Store result as a single row to prevent index errors
                        newRow = {subjName, currentAge, omegaVal};
                        compiledData = [compiledData; newRow]; 
                        
                        fprintf('Processed %s | ID: %d | Age: %d | Omega: %.2f\n', ...
                            subjName, currentID, currentAge, omegaVal);
                    end
                end
            end
        end
    catch ME
        fprintf('Error at %s: %s\n', subjName, ME.message);    
    end
end

% Post-Processing and Statistical Analysis
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectName', 'Age', 'Omega'});
    
    % Perform Pearson Correlation
    [R, P] = corrcoef(results.Age, results.Omega);
    r_val = R(1,2);
    p_val = P(1,2);
    
    % Plotting
    figure('Color', 'w');
    scatter(results.Age, results.Omega);
    hLine = lsline; 
    set(hLine, 'LineWidth', 1, 'Color', 'r');
    
    xlabel('Post-birth Age');
    ylabel('O-Information');
    title(sprintf('Neonatal Development (Control): Age vs O-Info\nr = %.2f, p = %.4f', r_val, p_val));
    grid on;
    
    save('Neonatal_OInfo_Results.mat', 'results');
    fprintf('\nAnalysis complete. Correlation: r = %.2f (p = %.4f)\n', r_val, p_val);
else
    fprintf('\nNo subjects matched the "NEO" criteria with valid CSV data.\n');
end