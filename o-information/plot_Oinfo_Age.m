clearvars
close all

%filePath = fullfile('..', 'spreadsheet', 'AvgOmega_PCA_Space_LIGHT.csv');
%results = readtable(filePath);

load('AvgOmega_SSW_IUGR_SPONT_top10cluster.mat')
stimType = 'SPONT';

if iscell(results.SSW) || isstring(results.SSW)
    xData = str2double(results.SSW);
else
    xData = results.SSW;
end
yData = results.Omega;
% Create the figure
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;

% Scatter plot of the data
s = scatter(xData, yData, 'filled', 'MarkerFaceAlpha', 0.5);

% Add a Trend Line (Linear Fit)
if height(results) > 1
    % Fit: Omega = p1*SSW + p2
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    [y_fit, delta] = polyval(p, x_range, S);
    
    plot(x_range, y_fit, 'r-');
   
    % Calculate Correlation for the title
    [r, pval] = corr(xData, yData, 'Type', 'Spearman');

    if pval < 0.001
        subtitle(sprintf('rho = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    else
        subtitle(sprintf('rho = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    end
end

% Add reference line for Synergy vs Redundancy
%yline(0, '--k', 'LineWidth', 1, 'Label', 'Redundancy (+) / Synergy (-)');

% Labeling
xlabel('Gestational age (days)');
ylabel('O-information');

% Clean up axes
set(gca, 'FontSize', 12, 'TickDir', 'out');
title(stimType, 'FontWeight', 'bold');