clearvars
close all

finalData = readtable(fullfile('spreadsheet','AvgOmega_Model_LIGHT_norm.csv'));
%finalData = readtable(fullfile('spreadsheet','AvgOmega_Model_AUD_DEV_norm.csv'));
%finalData = readtable(fullfile('spreadsheet','AvgOmega_Model_AUD_STD_norm.csv'));
%finalData = readtable(fullfile('spreadsheet','AvgOmega_Model_SPONT_norm.csv'));

finalData.Group = categorical(finalData.Group);
finalData.Group = reordercats(finalData.Group, {'Control', 'SGA', 'IUGR'});

lme = fitlme(finalData, 'Omega ~ 1 + GA_days*Group + GA_days*Sex + (1|ID)');
disp(lme);


% figure('Color', 'w');
% plotResiduals(lme, 'fitted'); % Check for homoscedasticity
% title('Residuals vs Fitted Values');
% 
% figure('Color', 'w');
% % Plot raw data vs model predictions
% scatter(finalData.GA_days, finalData.Omega, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% plot(finalData.GA_days, fitted(lme), 'r.');
% xlabel('Gestational Age (Days)');
% ylabel('O-Information (\Omega)');
% legend('Observed', 'LME Predicted');
% title('Mixed Model Fit of Fetal Brain Integration');