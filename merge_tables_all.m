clearvars
close all

% List all the files you want to merge
filePaths = {
    '/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_LIGHT_norm_avgtrials.mat', ...
    '/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_AUD_STD_norm_avgtrials.mat', ...
    '/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_AUD_DEV_norm_avgtrials.mat', ...
    '/home/trongvule/Desktop/fetalMEG/o-information/AvgOmega_SSW_IUGR_SPONT_notrials_norm.mat'

    %'C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_AUD_DEV_15.mat', ...
    %'C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_AUD_STD_15.mat', ...
    %'C:\Users\idm\Desktop\fetalMEG\o-information\AvgOmega_SSW_IUGR_SPONT.mat',
    %"/home/trongvule/Desktop/fetalMEG/CTW/CTW_SSW_IUGR_LIGHT_AvgTrials_15.mat";
    %"/home/trongvule/Desktop/fetalMEG/LZ/LZ_SSW_IUGR_LIGHT_AvgTrials.mat";

};

allResults = table(); % Initialize an empty table to hold everything

% Loop through each file and stack them
for i = 1:length(filePaths)
    if exist(filePaths{i}, 'file')
        dataStruct = load(filePaths{i});
        
        % Assuming the table inside each .mat file is named 'results'
        currentTable = dataStruct.results;
        
        % Stack the table vertically
        allResults = [allResults; currentTable]; 
    else
        fprintf('Warning: File not found: %s\n', filePaths{i});
    end
end

% Process the combined table
tempTable = allResults;
tempTable.ID = regexp(tempTable.SubjectID, '^[^-]+', 'match', 'once');
tempTable.ID = categorical(tempTable.ID);

% Load the Spreadsheet
opts = detectImportOptions('spreadsheet/IUGR_demographics.xls');
opts.VariableNamingRule = 'preserve'; 
spreadsheet = readtable('spreadsheet/IUGR_demographics.xls', opts);

% Clean up Spreadsheet columns
spreadsheet.ID = categorical(spreadsheet.ID);
if ismember('Sex_(1=M,2=F)', spreadsheet.Properties.VariableNames)
    spreadsheet.Properties.VariableNames{'Sex_(1=M,2=F)'} = 'Sex';
end

% The Merge (Inner Join)
% This now matches ALL subjects from ALL mat files with the clinical data
finalData = innerjoin(tempTable, spreadsheet, 'Keys', 'ID');

% Set categorical groups
finalData.Group = categorical(finalData.Group);
finalData.Sex = categorical(finalData.Sex);

% Clean and Save
finalData = rmmissing(finalData);
writetable(finalData, 'spreadsheet/AvgOmega_full_avgtrials.csv');

disp('Merge Complete. Combined data contains:');
fprintf('%d total observations.\n', height(finalData));
disp(head(finalData, 5));