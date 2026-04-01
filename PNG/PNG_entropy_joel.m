% REMEMBER TO SET THE TAU, ~ LINE 20

close all
clearvars 
rng(5702)
try
    parpool(8)
catch
    fprintf('Parpool already open ..?\n')
end
load('infant_anesth_chan_locs.mat');
topoplot_output_dir = './Infant_PNG';
if ~exist(topoplot_output_dir, 'dir')
    mkdir(topoplot_output_dir);
end

if ~exist('./recyclemepls', 'dir')
    mkdir('./recyclemepls');
end
tic

tau = 4; % tau tag value (in ms)
fs = 250;
sampint = tau/(1000/fs);
nsymb = 64;

%% now for GENERAL baseline vs GENERAL anesthesia  ~~~~~~~~~~~~~
%  INPUT folders
baseline_general_input_dir = 'I:\infant_eeg_data\3_baseline_general_checkedPostproc_reref';
general_input_dir = 'I:\infant_eeg_data\4_anesth_general_checkedPostproc_reref';

in_dirs = {baseline_general_input_dir, general_input_dir};
for i = 1:length(in_dirs)
    if ~exist(in_dirs{i}, 'dir')
        error('does directory exist %s', in_dirs{i});
    end
end

% OUTPUT DIRS
baseline_general_output_dir = './Infant_PNG/BaselineGeneral/';
general_output_dir = './Infant_PNG/GeneralAnesthesia/';
png_matrix_baseline_general = './Infant_PNG/BaselineGeneral/';
png_matrix_general = './Infant_PNG/GeneralAnesthesia/';
out_dirs = {baseline_general_output_dir, general_output_dir, png_matrix_baseline_general, png_matrix_general};
for i = 1:length(out_dirs)
    if ~exist(out_dirs{i}, 'dir')
        mkdir(out_dirs{i});
    end
end

%% == loop through BASELINE general
baseline_general_files = dir(fullfile(baseline_general_input_dir, '*.set'));
for i = 1:length(baseline_general_files)
    EEG = pop_loadset('filename', baseline_general_files(i).name, 'filepath', baseline_general_input_dir);
    assert(EEG.srate == fs,'Wrong sampling rate')
    if EEG.nbchan ~= 19
        warning('how many channels in the recording???: %s', baseline_general_files(i).name);
    end
    if ~exist('./recyclemepls', 'dir')
        mkdir('./recyclemepls');
    end

    png = PNGentropy(EEG.data,sampint,nsymb,baseline_general_files(i).name(1:end-4));
    fprintf('\n Baseline general %1.2f%% complete\n',i/length(baseline_general_files)*100)
    save(sprintf('%s/%s_nsymb=%i_int=%i_%s',baseline_general_output_dir,baseline_general_files(i).name,nsymb,sampint,'baseline_general_PNG.mat'), 'png');
end

toc
%% ANESTHESIA general now
% loop through ANESTHESIA general
general_files = dir(fullfile(general_input_dir, '*.set'));
for i = 1:length(general_files)
    EEG = pop_loadset('filename', general_files(i).name, 'filepath', general_input_dir);

    assert(EEG.srate == fs,'Wrong sampling rate')

    if EEG.nbchan ~= 19
        warning('how many channels in the recording???: %s', baseline_general_files(i).name);
    end
    png = PNGentropy(EEG.data,sampint,nsymb,general_files(i).name(1:end-4));
    fprintf('\n General anesthesia %1.2f%% complete\n',i/length(general_files)*100)
    save(sprintf('%s/%s_nsymb=%i_int=%i_%s',general_output_dir,general_files(i).name,nsymb,sampint,'general_anesthesia__PNG.mat'), 'png');
end

toc

%% Paired t-test per channel across participants. Also Multiple comparisons correct

% get common files betw the 2 groups ~~~~~~~~~~~~~~~~~
baseline_files = dir(fullfile(baseline_general_output_dir, sprintf('*int=%i_baseline_general_PNG.mat',sampint)));
general_files = dir(fullfile(general_output_dir,  sprintf('*int=%i_general_anesthesia__PNG.mat',sampint)));

getID = @(filename) strsplit(filename, '_');
base_ids = cell(1, length(baseline_files));
for i = 1:length(baseline_files)
    parts = getID(baseline_files(i).name);
    base_ids{i} = parts{2};
end
general_ids = cell(1, length(general_files));
for i = 1:length(general_files)
    parts = getID(general_files(i).name);
    general_ids{i} = parts{2};
end

[matched_ids, base_idx, general_idx] = intersect(base_ids, general_ids);
matched_general = general_files(general_idx);
matched_baseline = baseline_files(base_idx);

disp('Matched files:');
for i = 1:length(matched_ids)
    fprintf('ID %s: %s <-> %s\n', matched_ids{i}, matched_baseline(i).name, matched_general(i).name);
end

num_participants = length(matched_ids);
num_channels = 19;

baseline_data = NaN(num_channels, num_participants);
general_data = NaN(num_channels, num_participants);

for i = 1:num_participants
    base = load(fullfile(baseline_general_output_dir, matched_baseline(i).name));
    general = load(fullfile(general_output_dir, matched_general(i).name));
    baseline_data(:,i) = base.png(:,1);  % shouldn't be more than one column
    general_data(:,i) = general.png(:,1); % shouldn't be more than one column
end


for i = 1:min(3, num_participants)  % Check first 3
    fprintf('Participant %d: baseline mean=%.4f, general mean=%.4f\n', ...
        i, mean(baseline_data(:,i)), mean(general_data(:,i)));
end



% arrays for stats
tvals = NaN(num_channels, 1);
pvals = NaN(num_channels, 1);
dof = num_participants - 1;

% t test for all chans
for ch = 1:num_channels
    [h, p, ci, stats] = ttest(baseline_data(ch,:), general_data(ch,:));
    tvals(ch) = stats.tstat;
    pvals(ch) = p;
end

baseline_grand = mean(baseline_data, 1);  % 1 × num_participants
general_grand = mean(general_data, 1);    % 1 × num_participants

[~, p_grand, ci_grand, stats_grand] = ttest(baseline_grand, general_grand);
fprintf('\nGrand average (across all channels): t = %.3f, p = %.4f\n', stats_grand.tstat, p_grand);

grand_avg_general = struct('mean_diff', mean(baseline_grand) - mean(general_grand), ...
                           't', stats_grand.tstat, ...
                           'p', p_grand, ...
                           'df', stats_grand.df, ...
                           'ci', ci_grand);
save('grand_avg_general_PNG.mat', 'grand_avg_general');


qvals = mafdr(pvals,'bhfdr',1);
significant_channels = qvals <= 0.05;

disp('Channel  Raw p-value  FDR p-value  Significant');
n = length(pvals);
for ch = 1:n
    fprintf('%4d      %7.5f     %7.5f      %s\n', ...
        ch, pvals(ch), qvals(ch), char('*' * significant_channels(ch)));
end

%% T-values Topoplot
figure('Position', [100 60 700 700]);
topoplot(tvals, CL, ...
    'maplimits', [-2, 2], ...
    'colormap', parula(64), ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'whitebk', 'on');
title(sprintf('PNG: Baseline vs. General Anesthesia \n t-values, tau = %d (%d ms)', sampint, tau), 'FontSize', 18, 'FontWeight', 'bold');

c = colorbar;
ylabel(c, 't-value', 'FontSize', 24, 'FontWeight', 'bold');
set(gcf,'color','w'); set(gca,'color','w');
set(gca,'linewidth',3); set(gca,'TickDir','out');
set(gca,'FontSize',16); axis equal;
c.LineWidth = 2; c.FontSize = 20;
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]);
filename_base = 'PNG_general_tvalues';
print(gcf, fullfile(topoplot_output_dir, [filename_base sprintf('int=%i.png',sampint)]), '-dpng', '-r300');

%% P-values Topoplot (FDR corrected)
logqvals = -log10(qvals);
figure('Position', [100 60 700 700]);
topoplot(logqvals, CL, ...
    'maplimits', [0, 2], ...
    'colormap', parula(64), ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'whitebk', 'on');
title(sprintf('PNG: Baseline vs. General Anesthesia \n p-values, tau = %d (%d ms)', sampint, tau), 'FontSize', 18, 'FontWeight', 'bold');

c = colorbar;
ylabel(c, '-log10(FDR p)', 'FontSize', 24, 'FontWeight', 'bold');
set(gcf,'color','w'); set(gca,'color','w');
set(gca,'linewidth',3); set(gca,'TickDir','out');
set(gca,'FontSize',16); axis equal;
c.LineWidth = 2; c.FontSize = 18;
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]);
filename_base = 'baseline_general_vs_general_anesthesia_pvalues_FDR';
print(gcf, fullfile(topoplot_output_dir, [filename_base sprintf('int=%i.png',sampint)]), '-dpng', '-r300');

save('results_general_png.mat', 'tvals', 'pvals', 'qvals', 'significant_channels');

%% baseline SPINAL vs anesthesia SPINAL
% for INPUT folders containing the data. check they exist
baseline_spinal_input_dir = 'I:\infant_eeg_data\1_baseline_spinal_checkedPostproc_reref';
spinal_input_dir = 'I:\infant_eeg_data\2_anesth_spinal_checkedPostproc_reref';
in_dirs = {baseline_spinal_input_dir, spinal_input_dir}; %input dirs
for i = 1:length(in_dirs)
    if ~exist(in_dirs{i}, 'dir')
        error('This directory doesn''t exist: %s', in_dirs{i});
    end
end

%- set / create OUTPUT directories for PNG values
baseline_spinal_output_dir = './Infant_PNG/BaselineSpinal/';
spinal_output_dir = './Infant_PNG/SpinalAnesthesia/';
png_matrix_baseline_spinal = './Infant_PNG/BaselineSpinal/';
png_matrix_spinal = './Infant_PNG/SpinalAnesthesia/';
out_dirs = {baseline_spinal_output_dir, spinal_output_dir, png_matrix_baseline_spinal, png_matrix_spinal};
for i = 1:length(out_dirs)
    if ~exist(out_dirs{i}, 'dir')
        mkdir(out_dirs{i});
    end
end

%% == BASELINE spinal first
baseline_spinal_files = dir(fullfile(baseline_spinal_input_dir, '*.set'));
for i = 1:length(baseline_spinal_files)
    try
        EEG = pop_loadset('filename', baseline_spinal_files(i).name, 'filepath', baseline_spinal_input_dir);
    catch
        if contains(baseline_spinal_files(i).name,'._id_') % this isn't a real file
            continue
        end
    end

    assert(EEG.srate == fs,'Wrong sampling rate')

    if EEG.nbchan ~= 19
        warning('how many channels in the recording???: %s', baseline_general_files(i).name);
    end
    png = PNGentropy(EEG.data,sampint,nsymb,baseline_spinal_files(i).name(1:end-4));
    fprintf('\n Baseline spinal %1.2f%% complete\n',i/length(baseline_spinal_files)*100)
    save(sprintf('%s/%s_nsymb=%i_int=%i_%s',baseline_spinal_output_dir,baseline_spinal_files(i).name,nsymb,sampint,'baseline_spinal__PNG.mat'), 'png');
end

%% ANESTHESIA spinal then
% == loop through ANESTHESIA spinal
spinal_files = dir(fullfile(spinal_input_dir, '*.set'));
for i = 1:length(spinal_files)
    EEG = pop_loadset('filename', spinal_files(i).name, 'filepath', spinal_input_dir);
    assert(EEG.srate == fs,'Wrong sampling rate')

    if EEG.nbchan ~= 19
        warning('how many channels in the recording???: %s', baseline_general_files(i).name);
    end
    png = PNGentropy(EEG.data,sampint,nsymb,spinal_files(i).name(1:end-4));
    fprintf('\n Spinal anesthesia %1.2f%% complete\n',i/length(spinal_files)*100)
    save(sprintf('%s/%s_nsymb=%i_int=%i_%s',spinal_output_dir,spinal_files(i).name,nsymb,sampint,'spinal_anesthesia_PNG.mat'), 'png');
end

%% Paired t-test per channel across participants. Also Multiple comparisons correct

% get common files betw the 2 groups ~~~~~~~~~~~~~~~~~
baseline_files = dir(fullfile(baseline_spinal_output_dir, sprintf('*int=%i_baseline_spinal__PNG.mat',sampint)));
spinal_files = dir(fullfile(spinal_output_dir,  sprintf('*int=%i_spinal_anesthesia_PNG.mat',sampint)));

getID = @(filename) strsplit(filename, '_'); %split by
base_ids = cell(1, length(baseline_files));
for i = 1:length(baseline_files)
    parts = getID(baseline_files(i).name);
    base_ids{i} = parts{2};
end
spinal_ids = cell(1, length(spinal_files));
for i = 1:length(spinal_files)
    parts = getID(spinal_files(i).name);
    spinal_ids{i} = parts{2};
end

[matched_ids, base_idx, spinal_idx] = intersect(base_ids, spinal_ids); % IDs per group > common in both groups

% Get matched files only
matched_baseline = baseline_files(base_idx);
matched_spinal = spinal_files(spinal_idx);

disp('Matched files:');
for i = 1:length(matched_ids)
    fprintf('ID %s: %s <-> %s\n', matched_ids{i}, matched_baseline(i).name, matched_spinal(i).name);
end

num_participants = length(matched_ids);
num_channels = EEG.nbchan;

baseline_data = NaN(num_channels, num_participants);
spinal_data = NaN(num_channels, num_participants);

for i = 1:num_participants
    base = load(fullfile(baseline_spinal_output_dir, matched_baseline(i).name));
    spinal = load(fullfile(spinal_output_dir, matched_spinal(i).name));
    baseline_data(:,i) = base.png;  % shouldn't be more than one column
    spinal_data(:,i) = spinal.png;  % shouldn't be more than one column
end

tvals = NaN(num_channels, 1);
pvals = NaN(num_channels, 1);

for ch = 1:num_channels % t test for all chans
    [h, p, ci, stats] = ttest(baseline_data(ch,:), spinal_data(ch,:));
    tvals(ch) = stats.tstat;
    pvals(ch) = p;
end

baseline_grand = mean(baseline_data, 1);
spinal_grand = mean(spinal_data, 1);

[~, p_grand, ci_grand, stats_grand] = ttest(baseline_grand, spinal_grand);
fprintf('\nGrand average (across all channels): t = %.3f, p = %.4f\n', stats_grand.tstat, p_grand);
grand_avg_spinal = struct('mean_diff', mean(baseline_grand) - mean(spinal_grand), ...
                          't', stats_grand.tstat, ...
                          'p', p_grand, ...
                          'df', stats_grand.df, ...
                          'ci', ci_grand);
save('grand_avg_spinal_PNG.mat', 'grand_avg_spinal');

qvals = mafdr(pvals,'bhfdr',1);

significant_channels = qvals <= 0.05;

% Display results
disp('Channel  Raw p-value  FDR p-value  Significant');
n = length(pvals);
for ch = 1:n
    fprintf('%4d      %7.5f     %7.5f      %s\n', ...
        ch, pvals(ch), qvals(ch), char('*' * significant_channels(ch)));
end
%% figures generation - t values plot
figure('Position', [100 60 700 700]);  % Add figure size
topoplot(tvals, CL, ...
    'maplimits', [-2, 2], ...
    'colormap', parula(64), ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'whitebk', 'on');
title(sprintf('PNG: Baseline vs. Spinal Anesthesia \n t-values, tau = %d (%d ms)', sampint, tau), 'FontSize', 18, 'FontWeight', 'bold');

c = colorbar;  % Assign to variable
ylabel(c, 't-value', 'FontSize', 24, 'FontWeight', 'bold');  % Add colorbar label
set(gcf,'color','w'); set(gca,'color','w');  % Add these
set(gca,'linewidth',3); set(gca,'TickDir','out');  % Add these
set(gca,'FontSize',16); axis equal;  % Add these
c.LineWidth = 2; c.FontSize = 20;  % Add colorbar formatting
pos = get(gca, 'Position');  % Add position adjustment
set(gca, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]);  % Add position adjustment
filename_base = 'baseline_spinal_vs_spinal_anesthesia_tvalues';
print(gcf, fullfile(topoplot_output_dir, [filename_base sprintf('int=%i.png',sampint)]), '-dpng', '-r300');

%% === plot sig. channels p < 0.05) ~~~~
logfdrpvals = -log10(qvals);

figure('Position', [100 60 700 700]);  % Add this
topoplot(logfdrpvals, CL, ...
    'maplimits', [0, 2], ...
    'colormap', parula(64), ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'whitebk', 'on');
title(sprintf('PNG: Baseline vs. Spinal Anesthesia \n p-values, tau = %d (%d ms)', sampint, tau), 'FontSize', 18, 'FontWeight', 'bold');

c = colorbar;
ylabel(c, '-log10(FDR p)', 'FontSize', 24, 'FontWeight', 'bold');
set(gcf,'color','w'); set(gca,'color','w');
set(gca,'linewidth',3); set(gca,'TickDir','out');
set(gca,'FontSize',16); axis equal;
c.LineWidth = 2; c.FontSize = 18;
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1), pos(2)-0.05, pos(3), pos(4)]);

% Save in both formats
filename_base = 'baseline_spinal_vs_spinal_anesthesia_pvalues_FDR';
print(gcf, fullfile(topoplot_output_dir, [filename_base sprintf('int=%i.png',sampint)]), '-dpng', '-r300');

save('results_spinal_png.mat', 'tvals', 'pvals', 'qvals', 'significant_channels');

%%

function [png] = PNGentropy(dat0, sampint, nsymb, name)

dat = dat0(:,1:sampint:end,:); % downsample
clear dat0
warning('off', 'MATLAB:DELETE:FileNotFound');

png = nan(size(dat,1),1);
n = size(dat,2)*size(dat,3);
num = log(n)/log(nsymb); % same normalization as LZ76 (entropy rate)
%num = log(n)/log(nsymb + 1); % ask Pedro if we should do it this way?
normfact = num/n;
dpi =  1000;
dpm = dpi * 39.3701; % Convert inches to meters

% Store filenames for cleanup
filenames = cell(size(dat,1), 1);

% remove byte size of header by saving a blank image
try 
    file_info = dir('./recyclemepls/PNG_REF_IMAGE.png');
    refval = file_info.bytes;
catch
    refimg = uint8(zeros(10,10));
    filename = './recyclemepls/PNG_REF_IMAGE.png';
    imwrite(refimg, filename, 'png', 'XResolution', dpm, 'YResolution', dpm, 'ResolutionUnit', 'meter');
    file_info = dir('./recyclemepls/PNG_REF_IMAGE.png');
    refval = file_info.bytes;
end

parfor ich = 1:size(dat,1)
    X = squeeze(dat(ich,:,:));
    tmp = discretize(X,quantile(X(:),0:1/nsymb:1));
    img = uint8((tmp - 1) / (nsymb - 1) * 255);

    % Generate unique filename for each iteration
    % @Dimitrios update the file path for your computer 
    filename = sprintf('./recyclemepls/%s_CH%i_deleteme.png', name, ich);
    filenames{ich} = filename;

    % Write with DPI metadata
    imwrite(img', filename, 'png', 'XResolution', dpm, 'YResolution', dpm, 'ResolutionUnit', 'meter');
    % Get file size
    file_info = dir(filenames{ich});
    png(ich) = (file_info.bytes - refval) * normfact;
    fprintf('.')
end

% Clean up files after parfor loop
recycle('off'); % Turn off recycle for cleanup
for i = 1:length(filenames)
    if exist(filenames{i}, 'file')
        delete(filenames{i});
    end
end
recycle('on'); % Turn recycle back on
fprintf('\n')
close all
end

