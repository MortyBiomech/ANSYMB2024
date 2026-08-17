clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
classification_data_path = [data_path, '8_Classification\ROIs_features\'];


%% Creating the input structure for the Classification Learner App

% First: Load the ROI files which contain the subjects-ICs of the multiple
% k-means clustering process

% % Flexion
% epoch_type = 'ROIs_1_Flexion_all_epochs.mat';
% data = load(fullfile(classification_data_path, epoch_type));
% name = fieldnames(data);
% ROIs_Flexion_all_epochs = data.(name{1});
% 
% epoch_type = 'ROIs_1_Flexion_per_trial.mat';
% data = load(fullfile(classification_data_path, epoch_type));
% name = fieldnames(data);
% ROIs_Flexion_per_trial = data.(name{1});
% 
% % Extension
% epoch_type = 'ROIs_1_Extension_all_epochs.mat';
% data = load(fullfile(classification_data_path, epoch_type));
% name = fieldnames(data);
% ROIs_Extension_all_epochs = data.(name{1});
% 
% epoch_type = 'ROIs_1_Extension_per_trial.mat';
% data = load(fullfile(classification_data_path, epoch_type));
% name = fieldnames(data);
% ROIs_Extension_per_trial = data.(name{1});

% FlextoFlex
% epoch_type = 'ROIs_1_FlextoFlex_all_epochs.mat';
% data = load(fullfile(classification_data_path, epoch_type));
% name = fieldnames(data);
% ROIs_FlextoFlex_all_epochs = data.(name{1});

epoch_type = 'ROIs_0_FlextoFlex_per_trial_pressure_Based_Flx_Ext_separate.mat';
data = load(fullfile(classification_data_path, epoch_type));
name = fieldnames(data);
ROIs_FlextoFlex_per_trial = data.(name{1});



%% Initialize and fill the tables
per_trial_or_all_epochs = 'per_trial';
selected_ROI = ROIs_FlextoFlex_per_trial;
type = 'FlextoFlex'; % 'FlextoFlex' 'Extension' 'Flexion'
classes = 'P1P3P6'; % 'P1P3P6' 'P1P6' 'P1P3' 'P3P6'

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
            featureStruct = regionData{rowIdx, 3}; % 3: nonNorm RMS, 4: Norm RMS
            %% Three classes (P1, P3, P6)
            data_temp = [featureStruct.P1; featureStruct.P3; featureStruct.P6]; 
            %% Two Classes 
            % data_temp = [featureStruct.P1; featureStruct.P6]; 
            data = cat(2, data, data_temp);
            
            % for the case we extracted the flexion and extension features
            % together. PSD was calsulated for the whole FlextoFlex epochs.
            % featureNames = strcat(regionName, '_', {'delta', 'theta', 'alpha', 'beta', 'gamma'}, ...
            %     '_IC', string(icID));
            
            % for the case we extracted the flexion and extension features
            % separately and created a 1*10 feature vector instead of 1*5.
            featureNames = strcat(regionName, '_', ...
                {'delta_Flx', 'theta_Flx', 'alpha_Flx', 'beta_Flx', 'gamma_Flx', ...
                 'delta_Ext', 'theta_Ext', 'alpha_Ext', 'beta_Ext', 'gamma_Ext'}, ...
                '_IC', string(icID));

            columnNames = cat(2, columnNames, featureNames);
            ICs_in_regions(regionIdx) = ICs_in_regions(regionIdx) + 1;

        end
    end

    % Z-score normalization
    data_zscore_norm = zscore(data);

    % Step 3: Add the 'class' column at the end
    columnNames = cat(2, columnNames, 'Pressure');
    
    % Replace numeric classes with string labels ('P1', 'P3', 'P6') and make them categorical
    %% Three classes (P1, P3, P6)
    classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                   repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
                   repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    %% Two classes 
    % classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
    %                repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data_zscore_norm, 'VariableNames', columnNames(1:end-1));
    subjectTable.Pressure = classLabels; % Add the 'Pressure' column with categorical labels

    % Step 4: Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjectID), '_', type, '_classes_', classes, '_', per_trial_or_all_epochs], subjectTable);
    assignin('base', ['Subject_', num2str(subjectID), '_ICs_in_regions'], ICs_in_regions);

end


