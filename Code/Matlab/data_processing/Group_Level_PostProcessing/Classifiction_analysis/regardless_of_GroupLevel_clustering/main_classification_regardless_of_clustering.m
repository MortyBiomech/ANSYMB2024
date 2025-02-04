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
    ROIs.(all_STUDY_names{1, i}) = cell(N, 4);
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


    if subject_list(i) < 10 % data structure for subject 10 and above was saved differently
            
        % frequencies of the Power Spectral Density
        frequencies = data{1, 1}.EEG_stream.Preprocessed.Freq_Domain.Freqs;
        
        for j = 1:length(all_STUDY_files)
            
            subject_sets_indx = ...
                find(cell2mat(ROIs.(all_STUDY_names{1, j})(:,1)) == subject_list(i));
            subject_comps = cell2mat(ROIs.(all_STUDY_names{1, j})(subject_sets_indx, 2));
            if ~isempty(subject_comps)
                
                for k = 1:length(subject_sets_indx)
                    IC = subject_comps(k);
                    [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
                        ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
                        RMS_features_generator(condition_indices, data, IC, frequencies, per_trial_or_all_epochs);
                end
    
            end
    
        end

    else

        % frequencies of the Power Spectral Density
        [data, frequencies] = calculate_PSD(data, Trials_Info);
        for j = 1:length(all_STUDY_files)
            
            subject_sets_indx = ...
                find(cell2mat(ROIs.(all_STUDY_names{1, j})(:,1)) == subject_list(i));
            subject_comps = cell2mat(ROIs.(all_STUDY_names{1, j})(subject_sets_indx, 2));
            if ~isempty(subject_comps)
                
                for k = 1:length(subject_sets_indx)
                    IC = subject_comps(k);
                    [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
                        ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
                        RMS_features_generator(condition_indices, data, IC, frequencies, per_trial_or_all_epochs);
                end
    
            end
    
        end

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
extension = ['_', per_trial_or_all_epochs,'.mat']; 

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


%% Initialize and fill the tables
selected_ROI = ROIs;
type = 'FlextoFlex'; % 'FlextoFlex' 'Extension' 'Flexion'
classes = 'P1P6'; % 'P1P3P6' 'P1P6' 'P1P3' 'P3P6'

regions = fieldnames(selected_ROI); % Get all region names (e.g., Left_PreMot_SuppMot)

for subjectID = subject_list
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
            % data_temp = [featureStruct.P1; featureStruct.P3; featureStruct.P6]; 
            %% Two Classes 
            data_temp = [featureStruct.P1; featureStruct.P6]; 
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
    %% Three classes (P1, P3, P6)
    % classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
    %                repmat({'P3'}, size(featureStruct.P3, 1), 1); ...
    %                repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    %% Two classes 
    classLabels = [repmat({'P1'}, size(featureStruct.P1, 1), 1); ...
                   repmat({'P6'}, size(featureStruct.P6, 1), 1)];
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data, 'VariableNames', columnNames(1:end-1));
    subjectTable.Pressure = classLabels; % Add the 'Pressure' column with categorical labels

    % Step 4: Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjectID), '_', type, '_classes_', classes, '_', per_trial_or_all_epochs], subjectTable);
    assignin('base', ['Subject_', num2str(subjectID), '_ICs_in_regions'], ICs_in_regions);

end

