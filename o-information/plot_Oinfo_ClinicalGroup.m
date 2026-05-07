% Load all data
filePath = fullfile('..', 'spreadsheet', 'AvgOmega_Model_LIGHT.csv');
stimType = 'LIGHT';
results = readtable(filePath);

% Define the group order
conditions = {'Control', 'IUGR', 'SGA'};

% Ensure the Group column is categorical with the specified order
results.Group = categorical(results.Group, conditions);

% Colors: Control (Purple), IUGR (Yellow), SGA (Blue)
plotColors = [
    0.4940, 0.1840, 0.5560; % Purple (Control)
    0.9290, 0.6940, 0.1250; % Yellow (IUGR)
    0.0000, 0.4470, 0.7410; % Blue   (SGA)
];

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;

% Draw underlying distribution (Boxchart)
b = boxchart(results.Group, results.Omega, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.1, 'MarkerStyle', 'none');

% Overlay individual data points (Swarmchart)
s = swarmchart(results.Group, results.Omega, 40, results.Group, 'filled', 'MarkerFaceAlpha', 0.6);

% Apply the custom color palette
colororder(plotColors);

% Formatting the plot
title(stimType, 'FontSize', 14);
ylabel('O-information');
% Set axes properties for clean look
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'none');

% Display summary statistics in the Command Window
stats = groupsummary(results, 'Group', 'mean', 'Omega');
disp('Summary Statistics:');
disp(stats);

