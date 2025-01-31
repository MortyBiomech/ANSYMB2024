clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
ROIs_RMS_features_path = [data_path, '8_Classification\ROIs_features\'];
study_folder = [data_path, '7_STUDY'];


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
epoch_type = 'Epochs_Extension_based.mat'; % 'Epochs_FlextoFlex_based', 'Epochs_Flexion_based', 'Epochs_Extension_based'
features_from_epochs = 'Extension'; % 'FlextoFlex', 'Extension', 'Flexion'
per_trial_or_all_epochs = 'all_epochs'; % 'per_trial', 'all_epochs'
subject_list = 5:18;  % List of subject IDs
ROIs = ROIs_with_features(all_STUDY_names, all_STUDY_files, subject_list, ...
    epoch_type, features_from_epochs, data_path, main_project_folder, per_trial_or_all_epochs);


%% Load ROIs with RMS features
% epoch_type = 'FlextoFlex.mat';
% file_name = ['ROIs_0_', epoch_type];
% 
% data = load(fullfile(ROIs_RMS_features_path, file_name));
% name = fieldnames(data);
% data = data.(name{1});

