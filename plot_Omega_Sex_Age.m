results = readtable('spreadsheet/AvgOmega_Model_AUD_DEV.csv');
sex = {1, 2};

% Master conversion of SSW to numeric 
if iscell(results.GA_days) || isstring(results.GA_days)
    allX = str2double(results.GA_days);
else
    allX = results.GA_days;
end
allY = results.Omega; % Master Y data

figure('Color', 'w', 'Position', [100, 100, 1500, 500]);

for i = 1:2
    subplot(1, 2, i); % Create three separate panels
    hold on;
    
    % Use logical indexing on the MASTER data
    %isData = strcmp(results.Sex, sex{i});
    isData = results.Sex == sex{i};
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
    ylabel('Omega');
   
    if i == 1
        title("Male", 'FontSize', 14);
        subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    else
        title("Female", 'FontSize', 14);
        subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');

        %subtitle(sprintf('r = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    end
    %subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');

    %grid on;
    set(gca, 'FontSize', 11);
end
