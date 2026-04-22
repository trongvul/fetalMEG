results = readtable('spreadsheet/CTW_Model_AUD_DEV_Moser.csv');
stimType = 'LIGHT';
sex = {1, 2};
sexNames = {'Male', 'Female'};

% Define custom colors (RGB values between 0 and 1)
colorMale = [0.2, 0.4, 0.8];  
colorFemale = [0.9, 0.3, 0.5]; 

% Master conversion of SSW to numeric
if iscell(results.SSW) || isstring(results.SSW)
    allX = str2double(results.SSW);
else
    allX = results.SSW;
end
allY = results.CTW; 
allIDs = string(results.ID); % Ensure you have an 'id' column

figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

for i = 1:2
    subplot(1, 2, i);
    hold on;
    
    % Filter data for current sex
    isData = results.Sex == sex{i};
    xData = allX(isData);
    yData = allY(isData);
    idData = allIDs(isData);
    
    if isempty(xData) || all(isnan(xData)), continue; end

    % Select the color based on the loop index
    if i == 1
        currentCol = colorMale;
    else
        currentCol = colorFemale;
    end
    
    % Plot Within-Subject Lines
    uniqueIDs = unique(idData);
    for j = 1:length(uniqueIDs)
        currentID = uniqueIDs(j);
        idx = (idData == currentID);
        
        % Only connect if there is more than 1 measurement
        if sum(idx) > 1
            % Sort by xData (Gestational Age) to ensure lines don't zig-zag back
            [x_sorted, sortIdx] = sort(xData(idx));
            y_vals = yData(idx);
            y_sorted = y_vals(sortIdx);
            
            plot(x_sorted, y_sorted, '-', 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1);
        end
    end
   
    % Scatter plot (Overlayed on top of lines)
    scatter(xData, yData, 40, 'filled', 'MarkerFaceColor', currentCol, 'MarkerFaceAlpha', 0.6);
    
    % Trend Line
    [p, S] = polyfit(xData, yData, 1);
    x_range = linspace(min(xData), max(xData), 100);
    y_fit = polyval(p, x_range);
    plot(x_range, y_fit, 'r--', 'LineWidth', 1.5);
    
    % Stats
    [r, pval] = corr(xData, yData, 'Type', 'Spearman', 'Rows', 'complete');
  
    xlabel('Gestational age (weeks)');
    ylabel('CTW');
    title(sexNames{i}, 'FontSize', 14);    
    
    if pval < 0.001
        subtitle(sprintf('r = %.2f, p < 0.001', r), 'FontAngle', 'italic');
    else
        subtitle(sprintf('r = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
    end
    
    set(gca, 'FontSize', 11);
end
sgtitle(stimType, 'FontWeight', 'bold');