%% Batch fMEG cluster-validity check + save 10-channel signal per subject
% For every SSW subject in baseDir:
%   - rank channels by smoothed RMS, take top-Kpool, reduce to ncluster by
%     centroid erosion, gate on contiguity (1 component = accept).
%
% Saved INTO each subject's own folder  baseDir/<subj>/3_fMEG/ :
%   <subj>_selected10.mat   the 10-channel signal + metadata   (accepted only)
%   <subj>_clustercheck.png before/after plot                  (accepted & rejected)
%
% Saved centrally (summaryDir):
%   cluster_validity_results.csv / .mat   verdict for every subject
%
% The saved signal is RAW (unfiltered) selected-channel data -- filtering here
% was only for RMS-based selection. Filter downstream as your next analysis needs.

clear; clc;

%% -------------------- Parameters --------------------
P.nbutt     = 4;    % Butterworth order
P.lof       = 1;    % low cutoff  (Hz)
P.hif       = 10;   % high cutoff (Hz)
P.Kpool     = 15;   % candidate pool size
P.ncluster  = 10;   % final number of channels
P.threshold = 4;    % neighbor distance (cm)
P.ksmooth   = 7;    % kNN smoothing for RMS ranking

baseDir    = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";
summaryDir = fullfile(pwd, 'cluster_validity_summary');   % central summary location
if ~exist(summaryDir, 'dir'); mkdir(summaryDir); end

allSubjectsDir = dir(baseDir);
numSubjects    = length(allSubjectsDir);

%% -------------------- Batch loop --------------------
results = struct('name', {}, 'status', {}, 'ncomp', {}, ...
                 'ncomp_pool', {}, 'nfinal', {}, 'message', {});

for i = 1:numSubjects
    subjName = allSubjectsDir(i).name;

    if ~(allSubjectsDir(i).isdir && contains(subjName, "SSW", 'IgnoreCase', true))
        continue
    end

    subjFMEGdir = fullfile(baseDir, subjName, '3_fMEG');
    datadir     = fullfile(subjFMEGdir, 'fMEG_data.ds');

    entry = struct('name', subjName, 'status', 'error', ...
                   'ncomp', NaN, 'ncomp_pool', NaN, 'nfinal', NaN, 'message', '');

    try
        if ~exist(datadir, 'dir')
            error('fMEG_data.ds not found at %s', datadir);
        end

        % Load raw data
        cfg         = [];
        cfg.dataset = char(datadir);
        cfg.channel = 'MEG';
        dataStruct  = ft_preprocessing(cfg);

        % Assess
        out = assessSubject(dataStruct, P);

        entry.ncomp      = out.ncomp;
        entry.ncomp_pool = out.ncomp_pool;
        entry.nfinal     = out.nfinal;
        entry.status     = ternary(out.accept, 'accept', 'reject');

        % Make sure the target directory exists (it should -- .ds lives there)
        if ~exist(subjFMEGdir, 'dir'); mkdir(subjFMEGdir); end

        % --- Save the before/after plot (accepted AND rejected) ---
        plotPath = fullfile(subjFMEGdir, [subjName '_clustercheck.png']);
        saveClusterPlot(out, subjName, entry.status, P, plotPath);

        % --- Save the 10-channel SIGNAL for accepted subjects ---
        if out.accept
            sel            = struct();
            sel.subject    = subjName;
            sel.label      = out.labels;      % 10x1 cell, channel names
            sel.signal     = out.signal;      % ncluster x nSamples, RAW signal
            sel.time       = out.time;        % 1 x nSamples, time vector (s)
            sel.fsample    = out.fsample;     % sampling rate (Hz)
            sel.coords     = out.coords;      % ncluster x 3, sensor positions
            sel.idx        = out.top;         % indices into the aligned MEG set
            sel.rms        = out.rms;         % raw RMS (ranking)
            sel.rms_smooth = out.rms_smooth;  % smoothed RMS (ranking)
            sel.params     = P;
            save(fullfile(subjFMEGdir, [subjName '_selected10.mat']), 'sel', '-v7.3');
        end

        fprintf('[%3d] %-24s %-7s  (final comps: %d, pool comps: %d)\n', ...
                i, subjName, upper(entry.status), out.ncomp, out.ncomp_pool);

    catch ME
        entry.status  = 'error';
        entry.message = ME.message;
        fprintf(2, '[%3d] %-24s ERROR   (%s)\n', i, subjName, ME.message);
    end

    results(end+1) = entry; %#ok<SAGROW>
end

%% -------------------- Summary --------------------
statuses = {results.status};
nAccept  = sum(strcmp(statuses, 'accept'));
nReject  = sum(strcmp(statuses, 'reject'));
nError   = sum(strcmp(statuses, 'error'));
nTotal   = numel(results);

fprintf('\n================ SUMMARY ================\n');
fprintf('Subjects processed : %d\n', nTotal);
fprintf('  Qualified (accept) : %d\n', nAccept);
fprintf('  Rejected  (scatter): %d\n', nReject);
fprintf('  Errors             : %d\n', nError);
if (nAccept + nReject) > 0
    fprintf('  Qualified rate     : %.1f%% (of %d processed)\n', ...
            100 * nAccept / (nAccept + nReject), nAccept + nReject);
end

rej = {results(strcmp(statuses, 'reject')).name};
if ~isempty(rej); fprintf('\nRejected subjects:\n'); fprintf('  %s\n', rej{:}); end
err = {results(strcmp(statuses, 'error')).name};
if ~isempty(err); fprintf('\nErrored subjects:\n');  fprintf('  %s\n', err{:}); end

resultsTable = struct2table(results);
writetable(resultsTable, fullfile(summaryDir, 'cluster_validity_results.csv'));
save(fullfile(summaryDir, 'cluster_validity_results.mat'), 'results', 'resultsTable', 'P');

fprintf('\nPer-subject files written to each  baseDir/<subj>/3_fMEG/\n');
fprintf('Summary written to:\n  %s\n', summaryDir);
disp(' '); disp(resultsTable);

%% ======================================================================
%  Local functions
%% ======================================================================
function out = assessSubject(dataStruct, P)
% Per-subject cluster-validity assessment; also returns the selected signal.

    % align data channels with grad coordinates
    grad      = dataStruct.hdr.grad;
    [tf, loc] = ismember(dataStruct.label, grad.label);
    keep      = tf;
    chanpos   = grad.chanpos(loc(keep), :);
    fMEG      = dataStruct.trial{1}(keep, :);   % aligned RAW signal (chan x time)
    labels    = dataStruct.label(keep);
    nchan     = size(chanpos, 1);

    if nchan < P.ncluster
        error('Only %d aligned channels (< ncluster = %d).', nchan, P.ncluster);
    end

    % filter + per-channel RMS (filtering used ONLY for selection)
    [B, A]    = butter(P.nbutt, [P.lof P.hif] ./ (dataStruct.fsample/2));
    fMEG_filt = filtfilt(B, A, fMEG')';
    rms_chan  = rms(fMEG_filt, 2);

    % smooth RMS over the sensor graph (ranking only)
    D = squareform(pdist(chanpos));
    k = min(P.ksmooth, nchan);
    [~, order] = sort(D, 2, 'ascend');
    knn_idx    = order(:, 1:k);
    rms_smooth = mean(rms_chan(knn_idx), 2);

    % top-Kpool pool
    [~, ord] = sort(rms_smooth, 'descend');
    pool     = ord(1:min(P.Kpool, nchan));

    % reduce pool -> ncluster by centroid erosion
    sel = pool(:);
    while numel(sel) > P.ncluster
        C        = mean(chanpos(sel,:), 1);
        dcent    = sqrt(sum((chanpos(sel,:) - C).^2, 2));
        [~, far] = max(dcent);
        sel(far) = [];
    end
    top = sel(:);

    % contiguity gate
    Adj        = (D < P.threshold) & ~eye(nchan);
    ncomp      = numel(unique(conncomp(graph(Adj(top,  top)))));
    ncomp_pool = numel(unique(conncomp(graph(Adj(pool, pool)))));

    out.accept     = (ncomp == 1);
    out.ncomp      = ncomp;
    out.ncomp_pool = ncomp_pool;
    out.nfinal     = numel(top);
    out.top        = top;
    out.pool       = pool;
    out.labels     = labels(top);
    out.coords     = chanpos(top, :);
    out.rms        = rms_chan(top);
    out.rms_smooth = rms_smooth(top);
    out.chanpos    = chanpos;              % all aligned coords, for plotting
    out.signal     = fMEG(top, :);         % 10-channel RAW signal for next analysis
    out.time       = dataStruct.time{1};   % time vector (s)
    out.fsample    = dataStruct.fsample;   % sampling rate (Hz)
end

function saveClusterPlot(out, subjName, status, P, filepath)
% Build a before/after plot (invisible figure) and save it to filepath.
    [~, PCs] = pca(out.chanpos);
    if mean(PCs(:,1)) < 0; PCs(:,1) = -PCs(:,1); end
    x = PCs(:,2); y = PCs(:,1);

    pool = out.pool; top = out.top; dropped = setdiff(pool, top);

    fig = figure('Visible', 'off', 'Position', [100 100 1080 520]);

    subplot(1,2,1); hold on;
    scatter(x,       y,       36, [0.75 0.75 0.75]);
    scatter(x(pool), y(pool), 100, 'r', 'filled');
    title({sprintf('Before: top-%d RMS pool', P.Kpool), ...
           sprintf('%d component%s', out.ncomp_pool, repmat('s',1,out.ncomp_pool>1))});
    axis equal off; hold off;

    subplot(1,2,2); hold on;
    scatter(x,          y,          36, [0.75 0.75 0.75]);
    scatter(x(dropped), y(dropped),  90, [0.95 0.6 0.2], 'filled');
    scatter(x(top),     y(top),     110, 'r', 'filled');
    if strcmp(status, 'accept')
        verdict = 'ACCEPT (1 cluster)';  tcol = [0 0.5 0];
    else
        verdict = sprintf('REJECT (%d components)', out.ncomp);  tcol = [0.8 0 0];
    end
    title({sprintf('After: reduced to %d', P.ncluster), verdict}, 'Color', tcol);
    axis equal off; hold off;

    sgtitle(sprintf('%s   (%d \\rightarrow %d)', ...
            strrep(subjName,'_','\_'), P.Kpool, P.ncluster), 'Interpreter', 'tex');

    try
        exportgraphics(fig, filepath, 'Resolution', 150);  % R2020a+
    catch
        saveas(fig, filepath);                              % fallback
    end
    close(fig);
end

function r = ternary(cond, a, b)
    if cond; r = a; else; r = b; end
end