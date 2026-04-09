load('PCA_results.mat');

cond_cols = {'AUDIO', 'LIGHT', 'SPONT'};

% Master conversion of SSW to numeric (do this ONCE)
if iscell(results.SSW) || isstring(results.SSW)
    allX = str2double(results.SSW);
else
    allX = results.SSW;
end
allY = results.Explained1; % Master Y data

figure('Color', 'w', 'Position', [100, 100, 1500, 500]);

for i = 1:3
    subplot(1, 3, i); % Create three separate panels
    hold on;
    
    % Use logical indexing on the MASTER data
    isData = strcmp(results.Condition, cond_cols{i});
    
    % Extract local copies so we don't overwrite the masters
    xData = allX(isData);
    yData = allY(isData);
    
    % Skip if no data for this condition
    if isempty(xData) || all(isnan(xData)), continue; end
    
    % Scatter plot
    scatter(xData, yData, 40, 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Trend Line
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    [y_fit, delta] = polyval(p, x_range, S);
    
    plot(x_range, y_fit, 'r-', 'LineWidth', 2);
    
    % Confidence Interval
    %fill([x_range, fliplr(x_range)], [y_fit+2*delta, fliplr(y_fit-2*delta)], ...
    %    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    % Stats for title
    %[r, pval] = corr(xData, yData, 'Rows', 'complete');
    [r, pval] = corr(xData, yData, 'Type', 'Spearman', 'Rows', 'complete');
  
    xlabel('Gestational age (weeks)');
    ylabel('PC1 Variance (%)');
    %title(sprintf('%s\nr = %.2f, p = %.3f', cond_cols{i}, r, pval));
    title(cond_cols{i}, 'FontSize', 14);    
    if i == 3
        subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    else
        subtitle(sprintf('r = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    end
    %grid on;
    set(gca, 'FontSize', 11);
end
