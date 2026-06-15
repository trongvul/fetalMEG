clearvars
close all


filePath = fullfile('..', 'spreadsheet', 'AvgOmega_Barley_LIGHT_full.csv');
results = readtable(filePath);

% Filter for the last measurement per subject ---
% Sort by ID and then by SSW (Gestational Age)
results = sortrows(results, {'ID', 'SSW'});

% Find the last occurrence of each unique ID
[~, lastIdx] = unique(results.ID, 'last');

% Create a filtered table containing only the final measurements
finalResults = results(lastIdx, :);

groups = {'Control', 'IUGR', 'SGA'};
stimType = 'Stimulus';

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
allX = finalResults.MDI;
allY = finalResults.Omega; 
allIDs = string(finalResults.ID); % Ensure you have an 'id' column

figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

for i = 1:3
    subplot(1, 3, i);
    hold on;

    % Filter data for current sex
    isData = strcmp(finalResults.Group, groups{i});
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



%filePath = fullfile('..', 'spreadsheet', 'AvgOmega_Barley_AUD_DEV.csv');
%results = readtable(filePath);

% Calculate the mean for numeric columns
subjSummary = groupsummary(results, 'ID', 'mean', {'Omega', 'MDI'});

% Add the 'Group' back in. 
% Since Group is the same for all entries of an ID, we just take the first occurrence.
[~, uniqueIdx] = unique(results.ID, 'first');
idGroupMap = results(uniqueIdx, {'ID', 'Group'});

% Join them together
subjSummary = innerjoin(subjSummary, idGroupMap, 'Keys', 'ID');

% Rename columns for the rest of your script to work
subjSummary.Properties.VariableNames{'mean_Omega'} = 'Omega';
subjSummary.Properties.VariableNames{'mean_MDI'}   = 'MDI';

groups = {'Control', 'IUGR', 'SGA'};

% Colors: Control (Purple), IUGR (Yellow), SGA (Blue)
plotColors = {
    [0.4940, 0.1840, 0.5560], ... % Purple
    [0.9290, 0.6940, 0.1250], ... % Yellow
    [0.0000, 0.4470, 0.7410]      % Blue
    };

figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

for i = 1:3
    subplot(1, 3, i);
    hold on;
    
    % Filter summary table for current group
    isGroup = strcmp(string(subjSummary.Group), groups{i});
    xData = subjSummary.MDI(isGroup);
    yData = subjSummary.Omega(isGroup);
    
    % Remove NaNs
    validIdx = ~isnan(xData) & ~isnan(yData);
    xData = xData(validIdx);
    yData = yData(validIdx);
    
    if isempty(xData), continue; end
    
    % Scatter plot
    scatter(xData, yData, 60, 'filled', 'MarkerFaceColor', plotColors{i}, 'MarkerFaceAlpha', 0.6);
    
    % Trend Line
    if length(xData) > 1
        p = polyfit(xData, yData, 1);
        x_range = linspace(min(xData), max(xData), 100);
        y_fit = polyval(p, x_range);
        plot(x_range, y_fit, 'k--', 'LineWidth', 1.5);
        
        % Stats
        [r, pval] = corr(xData, yData, 'Type', 'Spearman', 'Rows', 'complete');
        
        if pval < 0.001
            subtitle(sprintf('rho = %.2f, p < 0.001', r), 'FontAngle', 'italic');
        else
            subtitle(sprintf('rho = %.2f, p = %.3f', r, pval), 'FontAngle', 'italic');
        end
    end
    
    xlabel('Mental Development Index (MDI)');
    ylabel('O-Information');
    title(groups{i}, 'FontSize', 14);    
    set(gca, 'FontSize', 11, 'TickDir', 'out');
    grid on;
end

sgtitle([stimType, ' Condition: Subject-Average O-Info vs. MDI'], ...
    'FontSize', 16, 'FontWeight', 'bold');
