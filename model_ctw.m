clearvars
close all

finalData = readtable('spreadsheet/CTW_Model_AUD_DEV_15.csv');

finalData.Group = categorical(finalData.Group);
finalData.Group = reordercats(finalData.Group, {'Control', 'SGA', 'IUGR'});

lme = fitlme(finalData, 'CTW ~ 1 + GA_days*Group + GA_days*Sex + (1|ID)');
disp(lme);