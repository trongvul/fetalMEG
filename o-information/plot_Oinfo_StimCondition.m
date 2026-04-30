% Load all data
spontRes = load('AvgOmega_SSW_IUGR_SPONT.mat');
lightRes = load("AvgOmega_SSW_IUGR_LIGHT_15.mat");
stdRes   = load("AvgOmega_SSW_IUGR_AUD_STD_15.mat");
%% 
devRes   = load("AvgOmega_SSW_IUGR_AUD_DEV_15.mat");

% Put tables into a cell array for easy looping
tables = {spontRes.results, lightRes.results, stdRes.results, devRes.results};
conditions = {'SPONT', 'LIGHT', 'AUD-STD', 'AUD-DEV'};

% Prepare Master Table for plotting
allData = table();
for i = 1:length(tables)
    tempTab = tables{i};
    % Add a column to identify the condition
    tempTab.Condition = repmat(string(conditions{i}), height(tempTab), 1);
    allData = [allData; tempTab];
end

% Ensure Condition is categorical and in the correct order for the plot
allData.Condition = categorical(allData.Condition, conditions);

% Settings
% Colors: SPONT (Purple), LIGHT (Yellow), AUD (Blue) - Adjusted to match your condition order
% 1. SPONT, 2. LIGHT, 3. AUD-STD, 4. AUD-DEV
plotColors = [
    0.4940, 0.1840, 0.5560; % Purple (Spont)
    0.9290, 0.6940, 0.1250; % Yellow (Light)
    0.0000, 0.4470, 0.7410; % Blue (Std)
    0.3010, 0.7450, 0.9330  % Lighter Blue (Dev)
];

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;

% Draw underlying distribution (Boxchart)
% We use 'allData' directly. Boxchart handles grouping by categorical x-axis.
b = boxchart(allData.Condition, allData.Omega, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.1, 'MarkerStyle', 'none');

% Overlay individual data points (Swarmchart)
% We pass the categorical column as X and Omega as Y
s = swarmchart(allData.Condition, allData.Omega, 40, allData.Condition, 'filled', 'MarkerFaceAlpha', 0.5);

% Apply the custom color palette
colororder(plotColors);

% Formatting
title('O-Information across stimulus conditions', 'FontSize', 14);
ylabel('Omega');
%xlabel('Stimulus Type');
%grid on;
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'none');

% Optional: Statistical Annotation
% Calculate means to display or use for reference
stats = groupsummary(allData, 'Condition', 'mean', 'Omega');
disp(stats);