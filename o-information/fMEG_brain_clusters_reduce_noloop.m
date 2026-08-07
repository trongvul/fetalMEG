%% fMEG head-position channel selection -- reduction variant
% Alternative to region growing:
%   1. Rank channels by band-limited, smoothed RMS
%   2. Take the top Kpool (=15) as a candidate pool
%   3. Reduce the pool to ncluster (=10) by repeatedly removing the channel
%      farthest from the pool centroid (centroid recomputed each step)
%   4. Use contiguity as a GATE: if the final 10 form one connected cluster
%      -> accept; if they are still scattered -> REJECT the dataset
%
% No region growing here -- selection is purely RMS ranking + geometric
% compaction, and the contiguity test is a quality check, not a repair step.
%
% Requires: FieldTrip, Signal Processing Toolbox, Statistics Toolbox, pca.

clear; clc;

% I008-AER-AUDIO-01-34SSW_Rossitza_20100215_01_FAIRY

%% -------------------- Parameters --------------------
datadir   = '/server/fo2-22/data/FETAL_USER/FAIRY_THEM_ALL_2024_withICA/IUGR_24_ICA/I008-AER-AUDIO-01-34SSW_Rossitza_20100215_01_FAIRY/3_fMEG/fMEG_data.ds';  % point to the .ds directory, not .meg4

nbutt     = 4;    % Butterworth order
lof       = 1;    % low cutoff  (Hz)
hif       = 10;   % high cutoff (Hz)

Kpool     = 15;   % candidate pool size (top-RMS channels)
ncluster  = 10;   % final number of channels after reduction

threshold = 4;    % neighbor distance (cm): from the local min in pdist histogram
ksmooth   = 7;    % kNN smoothing for the RMS ranking (channel + neighbors)

%% -------------------- 1. Load --------------------
cfg         = [];
cfg.dataset = datadir;
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

%% -------------------- 2. Align data channels with grad coordinates -----
grad      = data.hdr.grad;
[tf, loc] = ismember(data.label, grad.label);
if ~all(tf)
    warning('%d data channels have no matching grad coordinate; dropping them.', sum(~tf));
end
keep    = tf;
chanpos = grad.chanpos(loc(keep),:);
fMEG    = data.trial{1}(keep,:);
labels  = data.label(keep);
nchan   = size(chanpos,1);
assert(size(fMEG,1) == nchan, 'Channel/coordinate mismatch after alignment.');

%% -------------------- 3. Filter and per-channel RMS --------------------
[B,A]     = butter(nbutt, [lof hif] ./ (data.fsample/2));
fMEG_filt = filtfilt(B, A, fMEG')';
rms_chan  = rms(fMEG_filt, 2);

%% -------------------- 4. Smooth RMS over the sensor graph (ranking) ----
D = squareform(pdist(chanpos));            % full pairwise sensor distances (cm)
k = min(ksmooth, nchan);
[~, order] = sort(D, 2, 'ascend');         % nearest first (col 1 = self)
knn_idx    = order(:, 1:k);
rms_smooth = mean(rms_chan(knn_idx), 2);   % ranking only; data stays unsmoothed

%% -------------------- 5. Top-Kpool candidate pool ---------------------
[~, ord] = sort(rms_smooth, 'descend');
pool     = ord(1:min(Kpool, nchan));       % indices into the aligned set

%% -------------------- 6. Reduce Kpool -> ncluster by centroid erosion --
% Remove the most eccentric channel one at a time, recomputing the centroid
% after each removal so a departing outlier re-centers the remaining set.
% (mean(X,1) is exactly what kmeans(X,1) returns as its single centroid.)
sel = pool(:);
while numel(sel) > ncluster
    C     = mean(chanpos(sel,:), 1);                 % centroid of current set
    dcent = sqrt(sum((chanpos(sel,:) - C).^2, 2));   % distance to centroid
    [~, far] = max(dcent);                           % most eccentric member
    sel(far) = [];                                   % drop it
end
top = sel(:);

%% -------------------- 7. Contiguity GATE (accept / reject) ------------
Adj    = (D < threshold) & ~eye(nchan);
bins   = conncomp(graph(Adj(top, top)));
ncomp  = numel(unique(bins));
accept = (ncomp == 1);

if accept
    fprintf('ACCEPT: the %d reduced channels form one connected cluster.\n', ncluster);
else
    fprintf(2, 'REJECT: reduced channels are scattered (%d components).\n', ncomp);
    fprintf(2, '        -> flag this dataset for exclusion.\n');
end

% Component sizes, useful for judging how bad a rejection is
csize = accumarray(bins(:), 1);
fprintf('Component sizes: %s\n', mat2str(sort(csize', 'descend')));

% Which channels survived reduction
fprintf('Final channels:\n');
for i = 1:numel(top)
    fprintf('  %s  (RMS %.3g, smoothed %.3g)\n', ...
            labels{top(i)}, rms_chan(top(i)), rms_smooth(top(i)));
end

%% -------------------- 8. Plot: before vs after reduction --------------
[~, PCs] = pca(chanpos);
if mean(PCs(:,1)) < 0
    PCs(:,1) = -PCs(:,1);
end
x = PCs(:,2);  y = PCs(:,1);

dropped = setdiff(pool, top);              % the Kpool-ncluster removed channels

% Contiguity of the pool itself (before reduction), for the left-panel title
ncomp_pool = numel(unique(conncomp(graph(Adj(pool, pool)))));

figure('Position', [100 100 1080 520]);

% --- Before: the top-Kpool candidate pool ---
subplot(1,2,1); hold on;
scatter(x,        y,        36, [0.75 0.75 0.75]);           % all channels
scatter(x(pool),  y(pool), 100, 'r', 'filled');             % candidate pool
title({sprintf('Before: top-%d RMS pool', Kpool), ...
       sprintf('%d component%s', ncomp_pool, repmat('s',1,ncomp_pool>1))});
axis equal off; hold off;

% --- After: reduced set, dropped channels shown, with verdict ---
subplot(1,2,2); hold on;
scatter(x,          y,          36, [0.75 0.75 0.75]);           % all channels
scatter(x(dropped), y(dropped),  90, [0.95 0.6 0.2], 'filled');  % dropped from pool
scatter(x(top),     y(top),     110, 'r', 'filled');             % final cluster
if accept
    verdict = 'ACCEPT (1 cluster)';       tcol = [0 0.5 0];
else
    verdict = sprintf('REJECT (%d components)', ncomp);  tcol = [0.8 0 0];
end
title({sprintf('After: reduced to %d', ncluster), verdict}, 'Color', tcol);
legend({'all channels', sprintf('dropped (%d)', numel(dropped)), ...
        sprintf('final (%d)', numel(top))}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal');
axis equal off; hold off;

sgtitle(sprintf('fMEG channel reduction  %d \\rightarrow %d', Kpool, ncluster));