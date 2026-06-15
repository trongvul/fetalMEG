clearvars;
close all;
rng(329014901);

baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

allSubjectsDir = dir(baseDir);
numSubjects = length(allSubjectsDir);
fs = 610.3516; 

compiledData = {}; 

fprintf('Starting analysis...\n');

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;
    try
        dsPath = fullfile(baseDir, subjName, '3_fMEG', 'fMEG_data.ds');
        hdr = ft_read_header(dsPath);
        [b, a] = butter(4, [1 15] / (hdr.Fs/2), 'bandpass');

        randIdx = randperm(hdr.nChans, numChansToPick);
        cfg = [];
        cfg.dataset = dsPath;
        cfg.channel = hdr.label(randIdx);
        cfg.feedback = 'no';
        dataStruct = ft_preprocessing(cfg);
        
        % Filtering
        data_filt = filtfilt(b, a, dataStruct.trial{1}')';

        % O-Info Calculation
        [~, SIGMA] = data2gaussian(data_filt');
        omegaMatrix(i, s) = soinfo_from_covmat(SIGMA, size(data_filt, 2));


        
    catch ME

    end


end

