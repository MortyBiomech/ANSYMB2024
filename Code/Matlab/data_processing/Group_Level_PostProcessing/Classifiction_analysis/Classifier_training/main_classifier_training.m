clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
classification_data_path = [data_path, '8_Classification\'];


%% Creating the input structure for the Classification Learner App

% First: Load the ROI files which contain the subjects-ICs of the multiple
% k-means clustering process

cd([classification_data_path, 'ROIs_features'])
load("ROIs_1_Flexion.mat")
ROIs_Flexion = ROIs;
load("ROIs_0_Extension.mat")
ROIs_Extension = ROIs;
load("ROIs_0_FlextoFlex.mat")
ROIs_FlextoFlex = ROIs;


%% Initialize and fill the tables
selected_ROI = ROIs_Flexion;
type = 'Flexion';

regions = fieldnames(selected_ROI); % Get all region names (e.g., Left_PreMot_SuppMot)
subjectIDs = 5:18; % Subject IDs to process

for subjectID = subjectIDs
    % Step 1: Initialize an empty table with no columns
    columnNames = {};
    data = [];
    ICs_in_regions = zeros(1, length(regions));

    for regionIdx = 1:numel(regions)
        regionName = regions{regionIdx};
        regionData = selected_ROI.(regionName); % Extract data for the current region

        % Find rows corresponding to the current subject
        subjectRows = find(cell2mat(regionData(:, 1)) == subjectID);
        if isempty(subjectRows)
            continue; % Skip if no ICs for this subject in this region
        end

        % Step 2: Process each IC for the subject in the current region
        for rowIdx = subjectRows'
            icID = regionData{rowIdx, 2}; % IC ID (2nd column)
            featureStruct = regionData{rowIdx, 3}; % Structure with classes (P1, P3, etc.)
            data_temp = [featureStruct.P1; featureStruct.P3; featureStruct.P6];
            data = cat(2, data, data_temp);
            
            featureNames = strcat(regionName, '_', {'delta', 'theta', 'alpha', 'beta', 'gamma'}, ...
                '_IC', string(icID));
            columnNames = cat(2, columnNames, featureNames);
            ICs_in_regions(regionIdx) = ICs_in_regions(regionIdx) + 1;
        end
    end

    % Step 3: Add the 'class' column at the end
    columnNames = cat(2, columnNames, 'Pressure');
    
    % Replace numeric classes with string labels ('P1', 'P3', 'P6') and make them categorical
    classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                   repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
                   repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data, 'VariableNames', columnNames(1:end-1));
    subjectTable.Pressure = classLabels; % Add the 'Pressure' column with categorical labels

    % Step 4: Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjectID), '_', type], subjectTable);
    assignin('base', ['Subject_', num2str(subjectID), '_ICs_in_regions'], ICs_in_regions);

end



