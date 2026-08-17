clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
ROIs_features_path = [data_path, ...
    '8_Classification\ROIs_features\regardless_of_clustering\'];
subjects_classification_results_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\Group_Level_PostProcessing\', ...
    'Classifiction_analysis\regardless_of_GroupLevel_clustering\',...
    'subjects_classification_results\'];


%% load ROIs containing the PSD integrals for each Brain ICs per subject
subject_list = 5:18;

fileName = ['ROIs_0_FlextoFlex_per_trial_pressure_Based_Flx_Ext', ...
    '_separate_NoBrainClustering_PSD_integral.mat'];
load(fullfile(ROIs_features_path, fileName))
selected_ROI = ROIs;


%% Initialize and fill the tables
type = 'FlextoFlex'; % 'FlextoFlex' 'Extension' 'Flexion'
classes = 'P1P3P6'; % 'P1P3P6' 'P1P6' 'P1P3' 'P3P6'
per_trial_or_all_epochs = 'per_trial';

regions = fieldnames(selected_ROI); 

for subjectID = subject_list
    % Step 1: Initialize an empty table with no columns
    columnNames = {};
    data = [];
    % ICs_in_regions = zeros(1, length(regions));

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
            % %% Two Classes 
            % data_temp = [featureStruct.P1; featureStruct.P6]; 
            data = cat(2, data, data_temp);
            
            featureNames = strcat(regionName, '_', ...
                {'delta_Flx', 'theta_Flx', 'alpha_Flx', 'beta_Flx', 'gamma_Flx', ...
                 'delta_Ext', 'theta_Ext', 'alpha_Ext', 'beta_Ext', 'gamma_Ext'}, ...
                '_IC', string(icID));
            columnNames = cat(2, columnNames, featureNames);
            % ICs_in_regions(regionIdx) = ICs_in_regions(regionIdx) + 1;
        end
    end

    % Step 3: Add the 'class' column at the end
    columnNames = cat(2, columnNames, 'Label');
    
    % Replace numeric classes with string labels ('P1', 'P3', 'P6') and make them categorical
    %% Three classes (P1, P3, P6)
    classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                   repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
                   repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    % %% Two classes 
    % classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
    %                repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data, 'VariableNames', columnNames(1:end-1));
    subjectTable.Label = classLabels; % Add the 'Pressure' column with categorical labels

    % Step 4: Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjectID), '_', type, '_classes_', classes, '_', per_trial_or_all_epochs], subjectTable);
    % assignin('base', ['Subject_', num2str(subjectID), '_ICs_in_regions'], ICs_in_regions);

end

data = struct();
for i = 1:length(subject_list)
    data.(['subject_', num2str(subject_list(i))]) = ...
        eval(['Subject_', num2str(subject_list(i)), ...
        '_FlextoFlex_classes_', classes, '_per_trial']);
end

struct_Name = 'subjects_psdIntegral_noBrainCluster_classificationTables';
folder_path = subjects_classification_results_path;
if ~exist(folder_path, "dir")
    mkdir(folder_path)
end
save(fullfile(folder_path, [struct_Name, '.mat']), 'data');




%% Main Loop over subjects to run the classification
features_from_epochs = 'FlextoFlex';
classes_pressure_based = 'P1P3P6';
% classes_score_based = 'S1S3';

for i = 1:length(subject_list)

    %% Pressure based EEG features
    % Create the input table for the classification 
    subjectTable_EEG = data.(['subject_', num2str(subject_list(i))]);
    % [subjectTable_EEG, ICs_in_regions] = ...
    %     createInputTableforClassification(ROIs_score_Based, ...
    %     subject_list(i), classes_score_based);

    % Save the table in the subject's folder
    tableName_EEG = ['Subject_', num2str(subject_list(i)), '_', ...
        features_from_epochs, '_classes_', classes_pressure_based, '_', ...
        per_trial_or_all_epochs];
    folder_path = [subjects_classification_results_path, ...
        'sub-', num2str(subject_list(i))];
    if ~exist(folder_path, "dir")
        mkdir(folder_path)
    end
    save(fullfile(folder_path, [tableName_EEG, '.mat']), 'subjectTable_EEG');


    % train classifiers and export the best one
    X = subjectTable_EEG(:, 1:end-1);
    Y = subjectTable_EEG.Label;
    [bestModel, meanAccuracy, stdAccuracy, acc_full, acc_selected, featureTable] = ...
        runClassificationWithoutSMOTE(X, Y, subject_list(i), folder_path);

    final_struct = struct('bestModel', bestModel, ...
        'Mean_Accuracy', meanAccuracy, ...
        'STD_Accuracy', stdAccuracy, ...
        'Repeated_5CV_Full_Features_Acc', acc_full, ...
        'Repeated_5CV_Selected_Features_Acc', acc_selected, ...
        'Sorted_Feature_Table', featureTable);
    save([folder_path, '\final_results_5CV_', 'sub-', num2str(subject_list(i)), '.mat'], ...
        'final_struct')


end