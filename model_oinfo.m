clearvars
close all

%finalData = readtable('spreadsheet\AvgOmega_SSW_IUGR_AUD_DEV.csv');
%finalData = readtable('spreadsheet\AvgOmega_SSW_IUGR_AUD_STD.csv');
%finalData = readtable('spreadsheet\AvgOmega_SSW_IUGR_LIGHT.csv');
finalData = readtable('spreadsheet\AvgOmega_Model_SPONT_notrials.csv');

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