%% fMEG head-position channel selection
% Reads a CTF .ds recording, ranks channels by band-limited RMS smoothed
% over the sensor layout, and selects a contiguous cluster of ncluster
% channels by region growing (contiguous by construction, fixed count).
%
% Pipeline:
%   1. Load .ds with FieldTrip
%   2. Align MEG data channels with grad coordinates
%   3. Band-pass filter (zero-phase) and compute per-channel RMS
%   4. Smooth the RMS map over the sensor graph (ranking only)
%   5. Region-grow a contiguous cluster from the strongest seed
%   6. Validate connectivity as a sanity check
%   7. Plot the array with the selected cluster highlighted
%
% Requires: FieldTrip (on path), Signal Processing Toolbox (butter/filtfilt/rms),
%           Statistics Toolbox (pdist), plus pca for the plot.

clear; clc;

%% -------------------- Parameters (adjust as needed) --------------------
%% -------------------- Parameters (adjust as needed) --------------------
%baseDir = "/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA";

%allSubjectsDir = dir(baseDir);
%numSubjects = length(allSubjectsDir);

% I008-AER-AUDIO-01-34SSW_Rossitza_20100215_01_FAIRY
% I040-AVER-AUDIO-01-28SSW_Rossitza_20101108_01_FAIRY

datadir   = '/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA/I008-AER-AUDIO-01-34SSW_Rossitza_20100215_01_FAIRY/3_fMEG/fMEG_data.ds';  % point to the .ds directory, not .meg4

nbutt     = 4;    % Butterworth order
lof       = 1;    % low cutoff  (Hz)
hif       = 10;   % high cutoff (Hz)
ncluster  = 10;   % number of channels to select

threshold = 4;    % neighbor distance (cm): channels within this are neighbors
                  % taken from the local minimum in the histogram of pdist(where)

ksmooth   = 7;    % kNN smoothing: channel itself + (ksmooth-1) neighbors
                  % keep modest so you don't blur across the quiet band
                  % between two field lobes

%% -------------------- 1. Load --------------------
cfg         = [];
cfg.dataset = datadir;
cfg.channel = 'MEG';          % restrict to MEG so RMS ranking and coordinate
                              % lookup are on the same channel set
data = ft_preprocessing(cfg);

%% -------------------- 2. Align data channels with grad coordinates -----
% RMS is computed on rows of data.trial{1} (order = data.label).
% Coordinates come from grad.chanpos (order = grad.label).
% These must refer to the same channels in the same order.
grad = data.hdr.grad;

[tf, loc] = ismember(data.label, grad.label);
if ~all(tf)
    warning('%d data channels have no matching grad coordinate; dropping them.', ...
            sum(~tf));
end

keep      = tf;                       % channels present in both
chanpos   = grad.chanpos(loc(keep),:);% coordinates aligned to data rows
fMEG      = data.trial{1}(keep,:);    % data aligned to those coordinates
labels    = data.label(keep);
nchan     = size(chanpos,1);

assert(size(fMEG,1) == nchan, 'Channel/coordinate mismatch after alignment.');

%% -------------------- 3. Filter and per-channel RMS --------------------
[B,A]     = butter(nbutt, [lof hif] ./ (data.fsample/2));
fMEG_filt = filtfilt(B, A, fMEG')';    % zero-phase, channels x time preserved
rms_chan  = rms(fMEG_filt, 2);         % one RMS value per channel

%% -------------------- 4. Smooth RMS over the sensor graph --------------
% kNN average: replace each channel's RMS with the mean over itself and its
% nearest neighbors. Used for RANKING ONLY -- downstream analysis uses the
% original unsmoothed channel data.
D = squareform(pdist(chanpos));        % full pairwise sensor distances (cm)

k = min(ksmooth, nchan);
[~, order]  = sort(D, 2, 'ascend');    % nearest first (col 1 = self, dist 0)
knn_idx     = order(:, 1:k);           % nchan x k neighbor indices
rms_smooth  = mean(rms_chan(knn_idx), 2);

%% -------------------- 4b. Naive selection (for comparison) ------------
% The original approach: top ncluster channels by raw RMS, no smoothing and
% no contiguity constraint. This is the "before" set -- it may be split.
[~, ord_raw] = sort(rms_chan, 'descend');
top_raw      = ord_raw(1:min(ncluster, nchan));
% (connectivity of this set is checked in step 6, once Adj is built)

%% -------------------- 5. Region-growing selection ---------------------
% Neighbor graph over ALL channels, so we know the cluster frontier.
Adj = (D < threshold) & ~eye(nchan);

[~, seed] = max(rms_smooth);           % start at the strongest (smoothed) channel
cluster   = seed;

while numel(cluster) < ncluster
    nbrs = find(any(Adj(cluster,:), 1));   % channels bordering the cluster
    nbrs = setdiff(nbrs, cluster);         % ...not already in it
    if isempty(nbrs)
        warning(['Connected component around the seed has only %d channels; ' ...
                 'cannot reach %d.'], numel(cluster), ncluster);
        break
    end
    [~, kbest]     = max(rms_smooth(nbrs));% add strongest frontier channel
    cluster(end+1) = nbrs(kbest); %#ok<SAGROW>
end

top = cluster(:);                      % selected channel indices (into aligned set)

%% -------------------- 6. Validate connectivity (sanity check) ----------
% Region-grown set
bins          = conncomp(graph(Adj(top, top)));
all_connected = numel(unique(bins)) == 1;
if all_connected
    fprintf('Region-grown: %d channels form one connected cluster.\n', numel(top));
else
    warning('Region-grown set split into %d components (unexpected).', numel(unique(bins)));
end

% Naive raw top-N set (the "before")
bins_raw   = conncomp(graph(Adj(top_raw, top_raw)));
ncomp_raw  = numel(unique(bins_raw));
if ncomp_raw > 1
    fprintf('Naive raw top-%d: %d connected components (split!).\n', ncluster, ncomp_raw);
else
    fprintf('Naive raw top-%d: 1 connected component.\n', ncluster);
end

% Report which channels were chosen
fprintf('Selected channels:\n');
for i = 1:numel(top)
    fprintf('  %s  (RMS %.3g, smoothed %.3g)\n', ...
            labels{top(i)}, rms_chan(top(i)), rms_smooth(top(i)));
end

%% -------------------- 7. Plot: before vs after ------------------------
% Flatten the 3D array onto its two principal axes for a 2D view, with a
% consistent orientation across subjects. Both panels share this projection
% so positions are directly comparable.
[~, PCs] = pca(chanpos);
if mean(PCs(:,1)) < 0
    PCs(:,1) = -PCs(:,1);              % flip so orientation is stable
end
x = PCs(:,2);  y = PCs(:,1);          % 2D plotting coordinates

figure('Position', [100 100 1000 450]);

% --- Before: naive top-N by raw RMS (may be split) ---
subplot(1,2,1); hold on;
scatter(x,          y,          36, [0.6 0.6 0.6]);          % all channels
scatter(x(top_raw), y(top_raw), 100, 'r', 'filled');         % naive selection
title(sprintf('Before: raw top-%d RMS (%d component%s)', ...
      ncluster, ncomp_raw, repmat('s', 1, ncomp_raw > 1)));
axis equal off; hold off;

% --- After: region-grown contiguous cluster ---
subplot(1,2,2); hold on;
scatter(x,        y,        36, [0.6 0.6 0.6]);               % all channels
scatter(x(top),   y(top),   100, 'r', 'filled');             % selected cluster
scatter(x(seed),  y(seed),  120, 'b', 'filled');             % seed
title(sprintf('After: region-grown cluster (%d channels)', numel(top)));
axis equal off; hold off;

sgtitle('fMEG channel selection: naive vs contiguous');