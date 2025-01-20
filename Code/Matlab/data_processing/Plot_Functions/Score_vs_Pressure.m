clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];


%% Load Trials Info for all subjects
cd(epoched_data_path)
subjects_id = 5:18;

for i = subjects_id
    cd(['sub-', num2str(i)])
    load("Trials_Info.mat")
    assignin('base', ['Trials_Info_sub_', num2str(i)], Trials_Info)
    cd(epoched_data_path)
end


%% Step 1: Organize the Data
% Dynamically find variables with the pattern 'Trials_Info_sub_*'
vars = who('Trials_Info_sub_*');

% Initialize a cell array to store all data
all_data = cell(1, numel(vars));

% Loop through the variables and assign their values to all_data
for i = 1:numel(vars)
    all_data{i} = eval(vars{i});
end

SubjectID = [];
Condition = [];
Score = [];
for i = 1:numel(subjects_id)

    condition_indices = condition_indices_identifier(all_data{i}, str2double(vars{i}(17:end)));
    
    % P1
    score_id_count_P1 = [];
    for j = 1:length(condition_indices.P1)
        score_id_count_P1 = cat(2, score_id_count_P1, all_data{1, i}{1, condition_indices.P1(j)}.General.Score);
    end
    % P3
    score_id_count_P3 = [];
    for j = 1:length(condition_indices.P3)
        score_id_count_P3 = cat(2, score_id_count_P3, all_data{1, i}{1, condition_indices.P3(j)}.General.Score);
    end
    % P6
    score_id_count_P6 = [];
    for j = 1:length(condition_indices.P6)
        score_id_count_P6 = cat(2, score_id_count_P6, all_data{1, i}{1, condition_indices.P6(j)}.General.Score);
    end

    score_id_count = [score_id_count_P1, score_id_count_P3, score_id_count_P6];
    condition_id_count = [ones(1, length(condition_indices.P1)), ...
        3*ones(1, length(condition_indices.P3)), ...
        6*ones(1, length(condition_indices.P6))];
    subject_id_count = str2double(vars{i}(17:end))*ones(1, length(condition_id_count));


    SubjectID = cat(2, SubjectID, subject_id_count);
    Condition = cat(2, Condition, condition_id_count);
    Score = cat(2, Score, score_id_count);

end

SubjectID = SubjectID';
Condition = Condition';
Score = Score';
dataTable = table(SubjectID, Condition, Score);
dataTable.SubjectID = categorical(dataTable.SubjectID);
dataTable.Condition = categorical(dataTable.Condition);

%% Step 2: Fit the Linear Mixed-Effects Model
% Define the model
model1 = fitlme(dataTable, 'Score ~ Condition + (1|SubjectID)'); % Random intercept only
model2 = fitlme(dataTable, 'Score ~ Condition + (Condition|SubjectID)'); % Random slope + intercept

compare(model1, model2)

