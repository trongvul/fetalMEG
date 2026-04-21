clearvars
close all

finalData = readtable('spreadsheet/LZ_Model_LIGHT_Moser.csv');

finalData.Group = categorical(finalData.Group);
finalData.Group = reordercats(finalData.Group, {'Control', 'SGA', 'IUGR'});

lme = fitlme(finalData, 'LZ ~ 1 + GA_days*Group + GA_days*Sex + (1|ID)');
disp(lme);
