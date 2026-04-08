% Settings
% Colors: AUDIO (Blue), LIGHT (Yellow), SPONT (Purple)
colors = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560];

pc_cols = {'Explained1', 'Explained2', 'Explained3'};
titles = {'PC1 variance', 'PC2 variance', 'PC3 variance'};

% Create Figure
figure('Color', 'w', 'Position', [100, 100, 1300, 450]);
% Changed 'relaxed' to 'loose'
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'loose');

for i = 1:3
    nexttile;
    hold on;
    
    % Group data by condition
    conds = categorical(results.Condition);
    yData = results.(pc_cols{i});
    
    % Draw underlying distribution
    boxchart(conds, yData, 'BoxFaceColor', 'k', 'BoxFaceAlpha', 0.1, 'MarkerStyle', 'none');
    
    % Overlay individual data points
    % Note: If results.SSW is still a string, we convert to double for point sizing
%     if iscell(results.SSW) || isstring(results.SSW)
%         sz = str2double(results.SSW);
%     else
%         sz = results.SSW;
%     end
%     
%     % Scale sizes to look good on the plot (e.g., 20 to 100)
%     ptSizes = normalize(sz, 'range', [30, 150]);
    
    s = swarmchart(conds, yData, 30, conds, 'filled', 'MarkerFaceAlpha', 0.6);
    
    colororder(colors);
    title(titles{i});
    ylabel('Variance explained (%)');
    grid on;
    set(gca, 'FontSize', 11);
end

% Add a legend for the sizing (only on the last plot)
%legend('Distribution', 'Subjects (Size = Age)', 'Location', 'northeastoutside');

%sgtitle('PCA explained variance: Condition comparison', 'FontSize', 15, 'FontWeight', 'bold');