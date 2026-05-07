filePath = fullfile('..', 'spreadsheet', 'AvgOmega_BarleyScore.csv');
results = readtable(filePath);

groups = {'Control', 'IUGR', 'SGA'};
stimType = 'Pooled';

% Colors: Control (Purple), IUGR (Yellow), SGA (Blue)
plotColors = {
    [0.4940, 0.1840, 0.5560], % Purple (Control)
    [0.9290, 0.6940, 0.1250], % Yellow (IUGR)
    [0.0000, 0.4470, 0.7410]  % Blue   (SGA)
    };


% Conversion of SSW to numeric
% if iscell(results.MDI) || isstring(results.MDI)
%     allX = str2double(results.MDI);
% else
%     allX = results.MDI;
% end
allX = results.MDI;
allY = results.Omega; 
allIDs = string(results.ID); % Ensure you have an 'id' column

figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

for i = 1:3
    subplot(1, 3, i);
    hold on;

    % Filter data for current sex
    isData = strcmp(results.Group, groups{i});
    xData = allX(isData);
    yData = allY(isData);
    idData = allIDs(isData);

    if isempty(xData) || all(isnan(xData)), continue; end

    % Scatter plot
    scatter(xData, yData, 40, 'filled', 'MarkerFaceColor', plotColors{i}, 'MarkerFaceAlpha', 0.5);

    % Trend Line
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    y_fit = polyval(p, x_range);
    plot(x_range, y_fit, 'r--', 'LineWidth', 1.5);

    % Stats
    [r, pval] = corr(xData, yData, 'Type', 'Spearman', 'Rows', 'complete');

    xlabel('Mental Developement Index (MDI)');
    ylabel('O-Information');
    title(groups{i}, 'FontSize', 14);    

    if pval < 0.001
        subtitle(sprintf('rho = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    else
        subtitle(sprintf('rho = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    end

    set(gca, 'FontSize', 11);
end

sgtitle(stimType, 'FontSize', 16, 'FontWeight', 'bold');
