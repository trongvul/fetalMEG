% Choose results file used to merge with the demographics spreadsheet.
clearvars
close all

%load("C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_LIGHT_15.mat");
%load("C:\Users\idm\Desktop\fetalMEG\CTW\CTW_SSW_IUGR_AUD_DEV.mat");

%load("/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_SPONT_notrials.mat");
load("C:\Users\idm\OneDrive\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_LIGHT_15.mat")


tempTable = results;

tempTable.ID = regexp(tempTable.SubjectID, '^[^-]+', 'match', 'once');

% Convert ID to categorical for easier matching
tempTable.ID = categorical(tempTable.ID);

opts = detectImportOptions('spreadsheet/Barley.csv');
opts.VariableNamingRule = 'preserve'; 
spreadsheet = readtable('spreadsheet/Barley.csv', opts);

% Ensure the spreadsheet ID column is also categorical
spreadsheet.ID = categorical(spreadsheet.ID);

%The Merge (Inner Join)
%This matches rows where ID exists in BOTH tables
finalData = innerjoin(tempTable, spreadsheet, 'Keys', 'ID');

% Set categorical groups for the LME
finalData.Group = categorical(finalData.Group);

% Remove row with NaN
%finalData = rmmissing(finalData);

writetable(finalData, 'spreadsheet/AvgOmega_Barley_LIGHT_full.csv');

disp('Merge Complete. First 5 rows:');
disp(head(finalData, 5));