load('AvgOmega_SSW_IUGR_SPONT.mat');

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
s = scatter(xData, yData);

% Add a Trend Line (Linear Fit)
if height(results) > 1
    % Fit: Omega = p1*SSW + p2
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    [y_fit, delta] = polyval(p, x_range, S);
    
    plot(x_range, y_fit, 'r-');
    
    % Add a shaded 95% Confidence Interval
    %fill([x_range, fliplr(x_range)], [y_fit+2*delta, fliplr(y_fit-2*delta)], ...
    %    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Calculate Correlation for the title
    [r, pval] = corr(xData, yData);
    title_str = sprintf('Development of Fetal Brain Information Dynamics\n(r = %.2f, p = %.3f)', r, pval);
else
    title_str = 'Fetal Brain Information Dynamics';
end

% Add reference line for Synergy vs Redundancy
yline(0, '--k', 'LineWidth', 1, 'Label', 'Redundancy (+) / Synergy (-)');

% 5. Labeling
xlabel('Conceptional age (weeks)');
ylabel('O-information');
title(title_str);
%legend('Location', 'northeast');

% 6. Clean up axes
set(gca, 'FontSize', 12, 'TickDir', 'out');