%% Batch fMEG Analysis Loop

% Specify the parent directory containing subject folders
parent_dir = ('/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA');
output_base_dir = ('/server/fo2-22/data/FETAL_USER/Dimitrios/Test_fmeg');

if parent_dir == 0
    error('No directory selected');
end

% Initialize component counters
contiguous = 0;
two_islands = 0;
three_islands = 0;
three_plus_islands = 0;

% Get list of all folders in parent directory
subject_folders = dir(parent_dir);
subject_folders = subject_folders([subject_folders.isdir]); % Keep only directories
subject_folders(ismember({subject_folders.name}, {'.', '..'})) = []; % Remove . and ..

fprintf('Found %d subject folders. Processing...\n', length(subject_folders));

%% Loop through each subject folder
for s = 1:length(subject_folders)
    subject_name = subject_folders(s).name;
    subject_path = fullfile(parent_dir, subject_name);
   
    % Check if 3_fMEG folder exists
    fmeg_dir = fullfile(subject_path, '3_fMEG');
    if ~isfolder(fmeg_dir)
        fprintf('Skipping %s: 3_fMEG folder not found\n', subject_name);
        continue;
    end
   
    % Check if .ds file exists
    ds_files = dir(fullfile(fmeg_dir, '*.ds'));
    if isempty(ds_files)
        fprintf('Skipping %s: No .ds file found in 3_fMEG\n', subject_name);
        continue;
    end
   
    fprintf('\nProcessing %s...\n', subject_name);
   
    try
        % Point to the .ds directory, not the .meg4 file itself
        datadir = fullfile(fmeg_dir, ds_files(1).name);
       
        cfg = [];
        cfg.dataset = datadir;
       
        data = ft_preprocessing(cfg);
       
        %% find top (10) RMS in (1 - 10) Hz band (parameters to adjust)
       
        nbutt = 4; % 4th order Butterworth
        lof = 1;
        hif = 10;
        ncluster = 10;
        [B,A] = butter(nbutt,[lof hif]./(data.fsample/2));
       
       
        fMEG = data.trial{1};
        fMEG_filt = filtfilt(B,A,fMEG')';
        rms_fMEG = rms(fMEG_filt,2);
        [~,tmp] = sort(rms_fMEG,'descend');
        top = tmp(1:ncluster);
        where = data.hdr.grad.chanpos;
        winners = where(top,:);
        howfar = squareform(pdist(winners));
       
        % Build adjacency matrix (logical, exclude self-loops)
        threshold = 4; % all channels within 4 mm are neighbors --- this is not arbitrary, taken from local minimum in histogram of pdist(where)
        adj = (howfar < threshold) & ~eye(ncluster);

        % Create graph and find connected components
        G = graph(adj);
        bins = conncomp(G);

        % Check if all points are in one component
        num_components = numel(unique(bins));
        all_connected = num_components == 1;
       
        fprintf('  All connected: %d | Components: %d\n', all_connected, num_components);
       
        % Count components
        if num_components == 1
            contiguous = contiguous + 1;
        elseif num_components == 2
            two_islands = two_islands + 1;
        elseif num_components == 3
            three_islands = three_islands + 1;
        else
            three_plus_islands = three_plus_islands + 1;
        end
       
        % Plotting (not recommended in loop)
        figure('Visible', 'off')
        [~,PCs] = pca(where);
        switch sign(mean(PCs(:,1)))
            case -1
                scatter(PCs(:,2),-PCs(:,1))
                scatter(PCs(top,2),-PCs(top,1),100,'r','filled')
            case 1
                scatter(PCs(:,2),PCs(:,1))
                scatter(PCs(top,2),PCs(top,1),100,'r','filled')
        end
        title('fMEG channels (red = brain cluster candidate)')
        axis off
        hold on
       
        output_folder = fullfile(output_base_dir, ['test_' subject_name]);
        if ~isfolder(output_folder)
            try
                mkdir(output_folder);
            catch ME
                fprintf('  Failed to create folder %s: %s\n', output_folder, ME.message);
                continue;
            end
        end
        saveas(gcf, fullfile(output_folder, 'graph.jpg'));
        close(gcf);
       
    catch ME
        fprintf('  Error processing %s: %s\n', subject_name, ME.message);
    end
end

fprintf('\n=== Batch processing complete ===\n');
fprintf('\nComponent Summary:\n');
fprintf('  Contiguous (1 component):  %d subjects\n', contiguous);
fprintf('  Two islands (2 components):  %d subjects\n', two_islands);
fprintf('  Three islands (3 components):  %d subjects\n', three_islands);
fprintf('  3+ islands (4+ components):  %d subjects\n', three_plus_islands);