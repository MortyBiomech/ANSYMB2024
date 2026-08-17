clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

grouplevel_postprocessing_path = fullfile([main_project_folder, ...
    '\Code\Matlab\data_processing\Group_Level_PostProcessing']);
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processed_data_path = [data_path, '5_single-subject-EEG-analysis\'];
trialInfo_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
ROIs_RMS_features_path = [data_path, '8_Classification\ROIs_features\'];
study_folder = [data_path, '7_STUDY'];


%% load data containing the finalized Brain ICs per subject
subject_list = 5:18;
load(fullfile(grouplevel_postprocessing_path, ...
    'Brain_PotentialBrain_AcceptedPotentialBrain.mat'));

Brain_ICs = cell(size(subject_list));
for i = 1:numel(subject_list)
    Brain_ICs{i} = [cell2mat(ICs(i, 1)); cell2mat(ICs(i, 3))];
end


%% Create meta STUDY files
all_STUDY_names = cell(size(subject_list));
for i = 1:length(all_STUDY_names)
    all_STUDY_names{i} = ['TotalBrain_sub_', num2str(subject_list(i))];
end

all_STUDY_files = cell(size(all_STUDY_names));

for i = 1:length(all_STUDY_files)
    all_STUDY_files{1, i} = Brain_ICs{1, i};
end

ROIs = struct();
for i = 1:length(all_STUDY_names)
    N = length(all_STUDY_files{1, i});
    ROIs.(all_STUDY_names{1, i}) = cell(N, 3);
end

for i = 1:length(all_STUDY_names)
    sets = length(all_STUDY_files{1, i});
    comps = all_STUDY_files{1, i};
    ROIs.(all_STUDY_names{1, i})(1:sets, 1) = num2cell(subject_list(i));
    ROIs.(all_STUDY_names{1, i})(1:sets, 2) = num2cell(comps);
end


%% Create the table containing RMS features for each subject

epoch_type = 'Epochs_FlextoFlex_based.mat'; % 'Epochs_FlextoFlex_based', 
                                            % 'Epochs_Flexion_based', 
                                            % 'Epochs_Extension_based'
features_from_epochs = 'FlextoFlex';         % 'FlextoFlex', 
                                            % 'Extension', 
                                            % 'Flexion'
per_trial_or_all_epochs = 'per_trial';     % 'per_trial', 
                                            % 'all_epochs'

pressure_score = 'pressure';                % 'pressure' based
                                            % 'score' based

for i = 1:length(subject_list)
    
    disp(['Subject ', num2str(subject_list(i))])
    
    data = load(fullfile(trialInfo_path, ['sub-', num2str(subject_list(i))], epoch_type));
    [~, baseName, ~] = fileparts(epoch_type);
    data = data.(baseName);
    
    Trials_Info = load(fullfile(trialInfo_path, ['sub-', num2str(subject_list(i))], 'Trials_Info.mat'));
    [~, baseName, ~] = fileparts('Trials_Info.mat');
    Trials_Info = Trials_Info.(baseName);


    % Find all conditions indices in trials 
    % P1, P3, and P6
    condition_indices = condition_indices_identifier(Trials_Info, subject_list(i));


    % frequencies of the Power Spectral Density
    [data, frequencies] = calculate_PSD_Flx_Ext_separate(...
        data, Trials_Info, subject_list(i));

     
        
    subject_sets_indx = ...
        find(cell2mat(ROIs.(all_STUDY_names{1, i})(:,1)) == subject_list(i));
    subject_comps = cell2mat(ROIs.(all_STUDY_names{1, i})(subject_sets_indx, 2));
        
    for k = 1:length(subject_sets_indx)
        IC = subject_comps(k);
        % [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
        %     ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
        %     RMS_features_generator(condition_indices, data, IC, frequencies, per_trial_or_all_epochs);
        
        ROIs.(all_STUDY_names{1, i})(subject_sets_indx(k), 3) =  ...
             PSD_integral_feature_generator_Flx_Ext_separate(condition_indices, data, IC, frequencies, per_trial_or_all_epochs, pressure_score);
    end


    clear data
    clear Trials_Info
        

end



%% Save ROI (check if there is a older version, save a new one)
ROIs_features_path = [ROIs_RMS_features_path, 'regardless_of_clustering'];

% Define the base name and extension
baseName = 'ROIs';
suffix = features_from_epochs;
folder = ROIs_features_path; 
extension = ['_', per_trial_or_all_epochs, '_', ...
        pressure_score, '_Based','_Flx_Ext_separate_NoBrainClustering_PSD_integral.mat']; 
    

% Initialize version number
version = 0;
fileName = sprintf('%s_%d_%s%s', baseName, version, suffix, extension);
fullFilePath = fullfile(folder, fileName);

% Increment version until a unique file name is found
while exist(fullFilePath, 'file') == 2
    version = version + 1;
    fileName = sprintf('%s_%d_%s%s', baseName, version, suffix, extension);
    fullFilePath = fullfile(folder, fileName);
end

% Save your file
save(fullFilePath, 'ROIs'); 
disp(['File saved as: ' fileName]);

