clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
EMG_features_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\EMG_processing\structured_EMG_data'];
EEG_ROIs_path = [data_path, '8_Classification\ROIs_features'];


%% Create the EEG feature table (Score Based)

% FlextoFlex
epoch_type = 'ROIs_2_FlextoFlex_per_trial_ScoreBased.mat';
data = load(fullfile(EEG_ROIs_path, epoch_type));
name = fieldnames(data);
ROIs_FlextoFlex_per_trial_ScoreBased = data.(name{1});

% Initialize and fill the tables
per_trial_or_all_epochs = 'per_trial';
selected_ROI = ROIs_FlextoFlex_per_trial_ScoreBased;
type = 'FlextoFlex'; % 'FlextoFlex' 'Extension' 'Flexion'
classes = 'S1S2S3'; % 'S1S2S3' 'S1S2' 'S1S3' 'S2S3'

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
            %% Three classes (S1, S2, S3)
            data_temp = [featureStruct.S1; featureStruct.S2; featureStruct.S3]; 
            %% Two Classes 
            % data_temp = [featureStruct.S1; featureStruct.S3]; 
            data = cat(2, data, data_temp);
            
            featureNames = strcat(regionName, '_', {'delta', 'theta', 'alpha', 'beta', 'gamma'}, ...
                '_IC', string(icID));
            columnNames = cat(2, columnNames, featureNames);
            ICs_in_regions(regionIdx) = ICs_in_regions(regionIdx) + 1;
        end
    end

    % Step 3: Add the 'class' column at the end
    columnNames = cat(2, columnNames, 'Pressure');
    
    % Replace numeric classes with string labels ('S1', 'S2', 'S3') and make them categorical
    %% Three classes (S1, S2, S3)
    classLabels = [repmat({'S1'}, size(featureStruct.S1, 1), 1); ...
                   repmat({'S2'}, size(featureStruct.S2, 1), 1); ...
                   repmat({'S3'}, size(featureStruct.S3, 1), 1)];
    %% Two classes 
    % classLabels = [repmat({'S1'}, size(featureStruct.S1, 1), 1); ...
    %                repmat({'S3'}, size(featureStruct.S3, 1), 1)];
    
    classLabels = categorical(classLabels); % Convert to categorical
    
    % Create the table for the current subject
    subjectTable = array2table(data, 'VariableNames', columnNames(1:end-1));
    subjectTable.Score = classLabels; % Add the 'Pressure' column with categorical labels

    % Step 4: Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjectID), '_', type, '_classes_', classes, '_', per_trial_or_all_epochs, '_EEG'], subjectTable);
    assignin('base', ['Subject_', num2str(subjectID), '_ICs_in_regions'], ICs_in_regions);

end

for subjectID = subjectIDs
    table_Name = ['Subject_', num2str(subjectID), '_', type, ...
        '_classes_', classes, '_', per_trial_or_all_epochs, '_EEG'];
    temp_table = eval(table_Name);

    featureCols = temp_table.Properties.VariableNames(1:end-1); % Feature column names
    for j = 1:length(featureCols)
        temp_table.(featureCols{j}) = normalize(temp_table.(featureCols{j}), 'zscore');
    end

    assignin('base', ['Subject_', num2str(subjectID), '_', type, ...
        '_classes_', classes, '_', per_trial_or_all_epochs, '_EEG_normalized'], ...
        temp_table);
end



%% Load EMG Features and create the table for classification
subjects_list = [5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18];
for i = 1:numel(subjects_list)
    disp(['subject ', num2str(subjects_list(i))])

    filename = ['sub-', num2str(subjects_list(i)), '\sub-', num2str(subjects_list(i)),'_structured_EMG_data.mat'];
    data = load(fullfile(EMG_features_path, filename));
    name = fieldnames(data);
    data = data.(name{1});


    %% identify the conditions
    condition_per = 2; % 1 = pressure_based; 2 = score_based
    if condition_per == 1
        condition_indices = condition_indices_identifier_EMG(data);
    elseif condition_per == 2
        % other function that labels the data based on score clustering
        condition_indices = condition_indices_identifier_EMG_ScoreBased(data, subjects_list(i));
    end

    %% Creating the tables for classification

    per_trial_or_all_epochs = 'per_trial';
    type = 'FlextoFlex'; % 'FlextoFlex' 'Extension' 'Flexion'
    classes = 'S1S2S3'; % 'P1P3P6' 'P1P6' 'P1P3' 'P3P6'; 'S1S2S3' 'S1S2' 'S1S3' 'S2S3'

    % P1_count = 0;
    % P3_count = 0;
    % P6_count = 0;
    % P1_features = [];
    % P3_features = [];
    % P6_features = [];

    S1_count = 0;
    S2_count = 0;
    S3_count = 0;
    S1_features = [];
    S2_features = [];
    S3_features = [];

    for j = 1:length(data)

        if isempty(data{1, j}.Classification_Features.per_trial)
            
            if ismember(j, condition_indices.S1) % condition_indices.P1
                % P1_features = cat(1, P1_features, data{1, j}.Classification_Features.per_trial);
                % P1_count = P1_count + 1;
                S1_features = cat(1, S1_features, NaN(1, 8));
                S1_count = S1_count + 1;
            elseif ismember(j, condition_indices.S2) % condition_indices.P3
                % P3_features = cat(1, P3_features, data{1, j}.Classification_Features.per_trial);
                % P3_count = P3_count + 1;
                S2_features = cat(1, S2_features, NaN(1, 8));
                S2_count = S2_count + 1;
            elseif ismember(j, condition_indices.S3) % condition_indices.P6
                % P6_features = cat(1, P6_features, data{1, j}.Classification_Features.per_trial);
                % P6_count = P6_count + 1;
                S3_features = cat(1, S3_features, NaN(1, 8));
                S3_count = S3_count + 1;
            end

            continue
        end

        if ismember(j, condition_indices.S1) % condition_indices.P1
            % P1_features = cat(1, P1_features, data{1, j}.Classification_Features.per_trial);
            % P1_count = P1_count + 1;
            S1_features = cat(1, S1_features, data{1, j}.Classification_Features.per_trial);
            S1_count = S1_count + 1;
        elseif ismember(j, condition_indices.S2) % condition_indices.P3
            % P3_features = cat(1, P3_features, data{1, j}.Classification_Features.per_trial);
            % P3_count = P3_count + 1;
            S2_features = cat(1, S2_features, data{1, j}.Classification_Features.per_trial);
            S2_count = S2_count + 1;
        elseif ismember(j, condition_indices.S3) % condition_indices.P6
            % P6_features = cat(1, P6_features, data{1, j}.Classification_Features.per_trial);
            % P6_count = P6_count + 1;
            S3_features = cat(1, S3_features, data{1, j}.Classification_Features.per_trial);
            S3_count = S3_count + 1;
        end
        
    end

    % all_features = [P1_features; P3_features; P6_features];
    all_features = [S1_features; S2_features; S3_features];
    % Three-class labels (P1, P3, P6)
    % classLabels = [repmat({'P1'}, P1_count, 1); ...
    %                repmat({'P3'}, P3_count, 1); ...
    %                repmat({'P6'}, P6_count, 1)];
    classLabels = [repmat({'S1'}, S1_count, 1); ...
                   repmat({'S2'}, S2_count, 1); ...
                   repmat({'S3'}, S3_count, 1)];
    classLabels = categorical(classLabels); % Convert to categorical

    % table columns
    columnNames = {'Vastus_med_R_Flexsion_RMS', 'Vastus_med_R_Extension_RMS', ...
                   'Rectus_femoris_R_Flexsion_RMS', 'Rectus_femoris_R_Extension_RMS', ...
                   'Gastrocnemius_R_Flexsion_RMS', 'Gastrocnemius_R_Extension_RMS', ...
                   'Biceps_femoris_R_Flexsion_RMS', 'Biceps_femoris_R_Extension_RMS', ...
                   'Label'};

    % Create the table for the current subject
    subjectTable = array2table(all_features, 'VariableNames', columnNames(1:end-1));
    subjectTable.Label = classLabels; % Add the 'Pressure' column with categorical labels

    % Store the table in the workspace dynamically
    assignin('base', ['Subject_', num2str(subjects_list(i)), '_', type, '_classes_', classes, '_', per_trial_or_all_epochs, '_EMG'], subjectTable);
    % assignin('base', ['Subject_', num2str(subjects_list(i)), '_classes_count'], [P1_count, P3_count, P6_count]);
    assignin('base', ['Subject_', num2str(subjects_list(i)), '_classes_count'], [S1_count, S2_count, S3_count]);


end




%% Combine two sets of EEG and EMG features 
% find Nan rows in EMG tables
for i = 1:length(subjects_list)
    
    % Construct the variable name as a string
    varName_EMG = ['Subject_', num2str(subjects_list(i)), '_', type, ...
        '_classes_', classes, '_', per_trial_or_all_epochs, '_EMG'];
    % Use EVAL to read the variable
    tmp_EMG = eval(varName_EMG);

    % Construct the variable name as a string
    varName_EEG = ['Subject_', num2str(subjects_list(i)), '_', type, ...
        '_classes_', classes, '_', per_trial_or_all_epochs, '_EEG'];
    % Use EVAL to read the variable
    tmp_EEG = eval(varName_EEG);

    % Find rows where all first 8 columns are NaN in EMG table
    nanRows = all(ismissing(tmp_EMG(:, 1:8)), 2);  

    % Remove these rows from E tables
    new_EMG_table = tmp_EMG(~nanRows, :);
    new_EEG_table = tmp_EEG(~nanRows, :);

    % Check if the last columns are identical for each row and add to new column
    identicalLabels = new_EMG_table.Label == new_EEG_table.Score;

    % Make the new table with both EEG and EMG features. 
    final_table = [new_EMG_table(:, 1:end-1), new_EEG_table(:, 1:end-1)];

    % Add the last column of new_EMG_table for rows with identical labels
    final_table.Label = repmat(categorical(missing), height(final_table), 1); 
    final_table.Label(identicalLabels) = new_EMG_table.Label(identicalLabels);


    % Store the table in the workspace dynamically
    varName_EMG_EEG = ['Subject_', num2str(subjects_list(i)), '_', ...
        type, '_classes_', classes, '_', per_trial_or_all_epochs, '_EMG_EEG'];
    assignin('base', varName_EMG_EEG, final_table);
    
    
    
    % Count and display the number of missing values in the last column
    numMissing = sum(ismissing(final_table.Label));
    disp(['subject ', num2str(subjects_list(i)), ...
        ', Number of missing labels: ', num2str(numMissing)]);
    

end



