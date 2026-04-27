% Choose results file used to merge with the demographics spreadsheet.
clearvars
close all

%load("C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_LIGHT_15.mat");
%load("C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_AUD_DEV_15.mat");
%load("C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_AUD_STD_15.mat");
%load("C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_SPONT_notrials.mat");

%load("C:\Users\idm\Desktop\fetalMEG\CTW\CTW_SSW_IUGR_LIGHT_AvgTrials_15.mat");

%load("/home/trongvule/Desktop/fetalMEG/CTW/CTW_SSW_IUGR_AUD_DEV_Moser.mat");
%load("/home/trongvule/Desktop/fetalMEG/LZ/LZ_SSW_IUGR_LIGHT_Moser.mat")

load("/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_LIGHT_norm_avgtrials.mat");
%load("/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_SPONT_notrials_norm.mat");



tempTable = results;

tempTable.ID = regexp(tempTable.SubjectID, '^[^-]+', 'match', 'once');

% Convert ID to categorical for easier matching
tempTable.ID = categorical(tempTable.ID);

opts = detectImportOptions('spreadsheet/IUGR_demographics.xls');
opts.VariableNamingRule = 'preserve'; 
spreadsheet = readtable('spreadsheet/IUGR_demographics.xls', opts);

% Ensure the spreadsheet ID column is also categorical
spreadsheet.ID = categorical(spreadsheet.ID);

% Change the column's name
spreadsheet.Properties.VariableNames{'Sex_(1=M,2=F)'} = 'Sex';

idsInMEG = unique(tempTable.ID);
idsInClinical = unique(spreadsheet.ID);

missingInClinical = setdiff(idsInMEG, idsInClinical);
missingInMEG = setdiff(idsInClinical, idsInMEG);

if ~isempty(missingInClinical)
    fprintf('Warning: %d subjects in MEG have no clinical data\n', ...
        length(missingInClinical));
end

%The Merge (Inner Join)
%This matches rows where ID exists in BOTH tables
finalData = innerjoin(tempTable, spreadsheet, 'Keys', 'ID');

% Set categorical groups for the LME
finalData.Group = categorical(finalData.Group);
finalData.Sex = categorical(finalData.Sex);

% Remove row with NaN
finalData = rmmissing(finalData);

writetable(finalData, 'spreadsheet/AvgOmega_Model_LIGHT_norm_avgtrials.csv');

disp('Merge Complete. First 5 rows:');
disp(head(finalData, 5));