clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
ROIs_RMS_features_path = [data_path, '8_Classification\ROIs_features\'];
study_folder = [data_path, '7_STUDY'];


%% Load meta STUDY files
load(fullfile(study_folder, "all_STUDY_names.mat"));
load(fullfile(study_folder, "all_STUDY_files.mat"));


%% Creat a metafile containg all ROIs and Subjects/ICs RMS features
epoch_type = 'Epochs_FlextoFlex_based.mat'; % 'Epochs_FlextoFlex_based', 'Epochs_Flexion_based', 'Epochs_Extension_based'
features_from_epochs = 'FlextoFlex'; % 'FlextoFlex', 'Extension', 'Flexion'
per_trial_or_all_epochs = 'per_trial'; % 'per_trial', 'all_epochs'
subject_list = 5:18;  % List of subject IDs
pressure_score = 'score'; % 'pressure', 'score' : define the conditions for later classification
ROIs = ROIs_with_features_ScoreBased(all_STUDY_names, all_STUDY_files, subject_list, ...
    epoch_type, features_from_epochs, data_path, main_project_folder, per_trial_or_all_epochs, pressure_score);


%% Load ROIs with RMS features
% epoch_type = 'FlextoFlex.mat';
% file_name = ['ROIs_0_', epoch_type];
% 
% data = load(fullfile(ROIs_RMS_features_path, file_name));
% name = fieldnames(data);
% data = data.(name{1});

