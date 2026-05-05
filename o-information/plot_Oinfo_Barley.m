clearvars
close all

filePath = fullfile('..', 'spreadsheet', 'AvgOmega_BarleyScore.csv');
results = readtable(filePath);
condition = 'SGA';
results = rmmissing(results);

%if iscell(results.MDI) || isstring(results.MDI)
%    xData = str2double(results.MDI);
%else
%    xData = results.MDI;
%end
xData = results.MDI;
yData = results.Omega;

isData = strcmp(results.Group, condition); %& results.Sex == 1;
xData = xData(isData);
yData = yData(isData);

% Create the figure
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;

% Scatter plot of the data
s = scatter(xData, yData);

% Add a Trend Line (Linear Fit)
if height(results) > 1
    % Fit: Omega = p1*MDI/PDI + p2
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    [y_fit, delta] = polyval(p, x_range, S);
    
    plot(x_range, y_fit, 'r-');
    
    % Calculate Correlation for the title
    [r, pval] = corr(xData, yData);

    if pval < 0.001
        subtitle(sprintf('r = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    else
        subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    end
end

% Labeling
xlabel('Mental Development Index (MDI)');
ylabel('O-information');

% Clean up axes
set(gca, 'FontSize', 12, 'TickDir', 'out');
%title(stimType, 'FontWeight', 'bold');


%results.Group = categorical(results.Group);
%results.Group = reordercats(results.Group, {'Control', 'SGA', 'IUGR'});

lme = fitlme(results, 'Omega ~ 1 + MDI*Group + MDI*Sex + (1|ID)');
disp(lme);