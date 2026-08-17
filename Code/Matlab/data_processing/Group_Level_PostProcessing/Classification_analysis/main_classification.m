clc
clear


%% Add and Define Necessary Paths
main_project_folder = 'D:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
ROIs_RMS_features_path = [data_path, '8_Classification\ROIs_features\'];
study_folder = [data_path, '7_STUDY'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
subjects_classification_results_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\Group_Level_PostProcessing\', ...
    'Classification_analysis\subjects_classification_results\'];
EMG_features_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\EMG_processing\structured_EMG_data'];



%% Load ALLEEG 
% if ~exist('ALLCOM','var')
%     eeglab;
% end
% 
% % Define subjects
% subject_list = 5:18;  % List of subject IDs
% 
% for i = 1:length(subject_list)
%     file_name = ['sub-', num2str(subject_list(i)), '_cleaned_with_ICA.set'];
%     dataset_path = fullfile([processed_data_path, 'sub-', ...
%         num2str(subject_list(i)), filesep]);
%     EEG = pop_loadset('filename', file_name, 'filepath', dataset_path);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
% end


%% Create meta STUDY files
% all_STUDY_names = {'Left_PreMot_SuppMot', 'Left_Paracentral_Lobule', ...
%     'Left_Dorsal_ACC', 'Left_VisMotor', 'Left_PrimVisual', ...
%     'Right_PreMot_SuppMot', 'Right_VisMotor', 'Right_PrimVisual'};
% all_STUDY_files = cell(size(all_STUDY_names));
% 
% for study = 1:length(all_STUDY_files)
%     filepath = [study_folder, '\multiple_clustering\', ...
%         all_STUDY_names{1, study}];
%     filename = ['main_study_potential_brain_ICs_RV-15_', ...
%         all_STUDY_names{1, study}, '.study'];
%     [all_STUDY_files{1, study}, ~] = pop_loadstudy('filepath', filepath, 'filename', filename);
% end


%% Load meta STUDY files
% load(fullfile(study_folder, "all_STUDY_names.mat"));
% load(fullfile(study_folder, "all_STUDY_files.mat"));


%% Creat a metafile containg all ROIs and Subjects/ICs RMS features
epoch_type = 'Epochs_FlextoFlex_based.mat'; % 'Epochs_FlextoFlex_based', 'Epochs_Flexion_based', 'Epochs_Extension_based'
features_from_epochs = 'FlextoFlex'; % 'FlextoFlex', 'Extension', 'Flexion'
per_trial_or_all_epochs = 'per_trial'; % 'per_trial', 'all_epochs'
subject_list = 5:18;  % List of subject IDs

% RMS feature of flexion and extension parts are calculated separately as
% two features based on pressure conditions
% ROIs_pressure_Based = ROIs_with_features_Flx_Ext_separate(all_STUDY_names, all_STUDY_files, subject_list, ...
%     epoch_type, features_from_epochs, data_path, per_trial_or_all_epochs, 'pressure'); 
ROIs_pressure_Based = load(fullfile(ROIs_RMS_features_path, 'ROIs_1_FlextoFlex_per_trial_pressure_Based_Flx_Ext_separate.mat'));
name = fieldnames(ROIs_pressure_Based);
ROIs_pressure_Based = ROIs_pressure_Based.(name{1});


% RMS feature of flexion and extension parts are calculated separately as
% two features based on score clustering results
% ROIs_score_Based    = ROIs_with_features_Flx_Ext_separate(all_STUDY_names, all_STUDY_files, subject_list, ...
%     epoch_type, features_from_epochs, data_path, per_trial_or_all_epochs, 'score');
ROIs_score_Based = load(fullfile(ROIs_RMS_features_path, 'ROIs_0_FlextoFlex_per_trial_score_Based_Flx_Ext_separate.mat'));
name = fieldnames(ROIs_score_Based);
ROIs_score_Based = ROIs_score_Based.(name{1});


%% Main Loop over subjects to create the classifiers input and run the classification
classes_pressure_based = 'P1P3P6';
classes_score_based = 'S1S2S3';

for i = 2:length(subject_list)

    %% Score based EEG features
    % Create the input table for the classification 
    [subjectTable_EEG, ICs_in_regions] = ...
        createInputTableforClassification(ROIs_score_Based, ...
        subject_list(i), classes_score_based);

    % Store the table in the workspace and save it in the subject's folder
    tableName_EEG = ['Subject_', num2str(subject_list(i)), '_', ...
        features_from_epochs, '_classes_', classes_score_based, '_', ...
        per_trial_or_all_epochs];
    assignin('base', tableName_EEG, subjectTable_EEG);
    assignin('base', ['Subject_', num2str(subject_list(i)), ...
        '_ICs_in_regions'], ICs_in_regions);

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
        runClassificationWithSMOTE(X, Y, subject_list(i), folder_path);

    final_struct = struct('bestModel', bestModel, ...
        'Mean_Accuracy', meanAccuracy, ...
        'STD_Accuracy', stdAccuracy, ...
        'Repeated_5CV_Full_Features_Acc', acc_full, ...
        'Repeated_5CV_Selected_Features_Acc', acc_selected, ...
        'Sorted_Feature_Table', featureTable);
    save([folder_path, '\final_results_5CV_', 'sub-', num2str(subject_list(i)), '.mat'], ...
        'final_struct')



    %% Score based EMG features
    % Create the input table for the classification 
    subjectTable_EMG = ...
        createInputTableforClassification_EMG(subject_list(i), EMG_features_path, ...
        'score', epoched_data_path, classes_score_based);

    % Find rows where all 8 columns are NaN in EMG table
    nanRows = all(ismissing(subjectTable_EMG(:, 1:8)), 2);  
    % Remove these rows from EMG table
    subjectTable_EMG = subjectTable_EMG(~nanRows, :);

    % Z-score normalization
    temp_EMG = subjectTable_EMG(:, 1:end-1);
    temp_EMG = zscore(table2array(temp_EMG));
    subjectTable_EMG(:, 1:end-1) = array2table(temp_EMG);


    % Store the table in the workspace and save it in the subject's folder
    tableName_EMG = ['Subject_', num2str(subject_list(i)), '_', ...
        features_from_epochs, '_classes_', classes_score_based, '_', ...
        per_trial_or_all_epochs, '_EMG'];
    assignin('base', tableName_EMG, subjectTable_EMG);

    folder_path = [subjects_classification_results_path, 'Score_based_EMG\',...
        'sub-', num2str(subject_list(i))];
    if ~exist(folder_path, "dir")
        mkdir(folder_path)
    end
    save(fullfile(folder_path, [tableName_EMG, '.mat']), 'subjectTable_EMG');


    % train classifiers and export the best one
    X = subjectTable_EMG(:, 1:end-1);
    Y = subjectTable_EMG.Label;
    [bestModel, meanAccuracy, stdAccuracy, acc_full, acc_selected, featureTable] = ...
        runClassificationWithSMOTE(X, Y, subject_list(i), folder_path);

    final_struct = struct('bestModel', bestModel, ...
        'Mean_Accuracy', meanAccuracy, ...
        'STD_Accuracy', stdAccuracy, ...
        'Repeated_5CV_Full_Features_Acc', acc_full, ...
        'Repeated_5CV_Selected_Features_Acc', acc_selected, ...
        'Sorted_Feature_Table', featureTable);
    save([folder_path, '\final_results_5CV_', 'sub-', num2str(subject_list(i)), '.mat'], ...
        'final_struct')


    %% Score based - EEG & EMG Features

    % EEG feature table
    [subjectTable_EEG, ICs_in_regions] = ...
        createInputTableforClassification(ROIs_score_Based, ...
        subject_list(i), classes_score_based);

    % EMG feature table
    subjectTable_EMG = ...
        createInputTableforClassification_EMG(subject_list(i), EMG_features_path, ...
        'score', epoched_data_path, classes_score_based);

    % Combining EEG and EMG
    test_labels_in_EEG_and_EMG = ...
        subjectTable_EEG.Label == subjectTable_EMG.Label;
    if sum(test_labels_in_EEG_and_EMG) ~= size(subjectTable_EEG, 1)
        disp(['Something is wrong with the EEG or EMG table in subject ', num2str(subject_list(i))])
        continue
    end

    subjectTable_EEG_EMG = [subjectTable_EEG(:, 1:end-1), subjectTable_EMG(:, 1:end-1), subjectTable_EEG(:, end)];

    % EMG part in the combined table needs to be z-normalized
    % Find rows where all 8 columns are NaN in EMG table
    nanRows = all(ismissing(subjectTable_EMG(:, 1:8)), 2);  
    % Remove these rows from EMG table
    subjectTable_EEG_EMG = subjectTable_EEG_EMG(~nanRows, :);

    % Z-score normalization
    temp_EMG = subjectTable_EEG_EMG(:, size(subjectTable_EEG,2):end-1);
    temp_EMG = zscore(table2array(temp_EMG));
    subjectTable_EEG_EMG(:, size(subjectTable_EEG,2):end-1) = array2table(temp_EMG);


    % Store the table in the workspace and save it in the subject's folder
    tableName_EEG_EMG = ['Subject_', num2str(subject_list(i)), '_', ...
        features_from_epochs, '_classes_', classes_score_based, '_', ...
        per_trial_or_all_epochs, '_EEG_EMG'];
    assignin('base', tableName_EEG_EMG, subjectTable_EEG_EMG);

    folder_path = [subjects_classification_results_path, 'Score_based_EEG_EMG\',...
        'sub-', num2str(subject_list(i))];
    if ~exist(folder_path, "dir")
        mkdir(folder_path)
    end
    save(fullfile(folder_path, [tableName_EEG_EMG, '.mat']), 'subjectTable_EEG_EMG');


    % train classifiers and export the best one
    X = subjectTable_EEG_EMG(:, 1:end-1);
    Y = subjectTable_EEG_EMG.Label;
    [bestModel, meanAccuracy, stdAccuracy, acc_full, acc_selected, featureTable] = ...
        runClassificationWithSMOTE(X, Y, subject_list(i), folder_path);

    final_struct = struct('bestModel', bestModel, ...
        'Mean_Accuracy', meanAccuracy, ...
        'STD_Accuracy', stdAccuracy, ...
        'Repeated_5CV_Full_Features_Acc', acc_full, ...
        'Repeated_5CV_Selected_Features_Acc', acc_selected, ...
        'Sorted_Feature_Table', featureTable);
    save([folder_path, '\final_results_5CV_', 'sub-', num2str(subject_list(i)), '.mat'], ...
        'final_struct')






end