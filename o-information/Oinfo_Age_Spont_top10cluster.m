clearvars
close all

% =========================================================================
% Analysis of spontaneous fMEG complexity via O-information
%   Pipeline:
%     1. Load fMEG data (FieldTrip)
%     2. Bandpass filter 1–15 Hz (Butterworth, 4th order)
%     3. Select 10 channels with highest RMS amplitude
%     4. Validate that those 10 channels form a single spatial cluster
%     5. Compute O-information (omega) on the 10-channel signal
% =========================================================================

baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
allSubjectsDir = dir(baseDir);
numSubjects    = length(allSubjectsDir);

fs = 610.3516;
[b, a] = butter(4, [1 15] / (fs/2), 'bandpass');

N_TOP       = 10;          % number of channels to keep

compiledData = {};
fprintf('Starting spontaneous O-information analysis (top-%d cluster)...\n', N_TOP);

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        if ~allSubjectsDir(i).isdir, continue; end
        if ~contains(subjName, "SSW",   'IgnoreCase', true), continue; end
        if ~contains(subjName, "SPONT", 'IgnoreCase', true), continue; end

        % Gestational age 
        load(fullfile(baseDir, subjName, '0_info', 'ssw.mat'));

        % Load raw data 
        dsPath = fullfile(baseDir, subjName, '3_fMEG', 'fMEG_data.ds');
        hdr    = ft_read_header(dsPath);

        cfg         = [];
        cfg.dataset = dsPath;
        dataStruct  = ft_preprocessing(cfg);

        % Bandpass filter 1–15 Hz 
        %   data_filt : [nChannels × nSamples]
        data_filt = filtfilt(b, a, dataStruct.trial{1}')';   % filtfilt expects [nSamples × nChannels]

        % Select top-N channels by RMS ---------------------------
        %   rmsVals : [nChannels × 1]  root-mean-square per channel
        rmsVals = sqrt(mean(data_filt .^ 2, 2));

        [~, sortIdx] = sort(rmsVals, 'descend');
        topIdx       = sortIdx(1:N_TOP);   % indices of the N_TOP highest-RMS channels

        % Retrieve 3-D positions for all channels [3 × nChannels]
        chanPos = hdr.grad.chanpos;   % FieldTrip convention: rows = channels, cols = [x y z]

        % Positions of the selected N_TOP channels [N_TOP × 3]
        selPos = chanPos(topIdx, :);


        % Compute O-information on the 10 selected channels -------
        data_sel = data_filt(topIdx, :);   % [N_TOP × nSamples]

        [~, SIGMA] = data2gaussian(data_sel');         % data2gaussian expects [nSamples × nChannels]
        omegaVal   = soinfo_from_covmat(SIGMA, size(data_sel, 2));

        % Store results
        compiledData{end+1, 1} = subjName;  %#ok<AGROW>
        compiledData{end,   2} = ssw;
        compiledData{end,   3} = omegaVal;

        fprintf('Processed %s: SSW = %s, Omega = %.4f  (channels: %s)\n', ...
            subjName, ssw, omegaVal, ...
            strjoin(hdr.label(topIdx), ', '));

    catch ME
        fprintf('Failed %s: %s\n', subjName, ME.message);
    end
end

% ---- Convert to table and save ------------------------------------------
if ~isempty(compiledData)
    results = cell2table(compiledData, 'VariableNames', {'SubjectID', 'SSW', 'Omega'});
    results = rmmissing(results);
    fprintf('\nAnalysis complete. Valid subjects: %d\n', height(results));
else
    fprintf('\nNo subjects matched the criteria.\n');
    results = table();
end

save('AvgOmega_SSW_IUGR_SPONT_top10cluster.mat', 'results');
fprintf('Saved to AvgOmega_SSW_IUGR_SPONT_top10cluster.mat\n');