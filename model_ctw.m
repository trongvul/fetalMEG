clearvars
close all

finalData = readtable('spreadsheet\Fetal_Model_IUGR_CTW_Spont.csv');

finalData.Group = categorical(finalData.Group);
finalData.Group = reordercats(finalData.Group, {'Control', 'SGA', 'IUGR'});

lme = fitlme(finalData, 'CTW ~ 1 + GA_days*Group + GA_days*Sex + (1|ID)');
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