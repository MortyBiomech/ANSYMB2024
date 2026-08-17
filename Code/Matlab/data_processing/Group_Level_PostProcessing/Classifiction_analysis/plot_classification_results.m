clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
subjects_classification_results_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\Group_Level_PostProcessing\', ...
    'Classifiction_analysis\subjects_classification_results\'];


%% Read the accuracy values from EEG, EMG, and EEG+EMG cases
subjects_list = 5:18;

folderName_EEG = 'Score_based_EEG';
folderName_EMG = 'Score_based_EMG';
folderName_EEG_EMG = 'Score_based_EEG_EMG';
EEG_mean_Acc = zeros(size(subjects_list));
EEG_std_Acc = zeros(size(subjects_list));
EMG_mean_Acc = zeros(size(subjects_list));
EMG_std_Acc = zeros(size(subjects_list));
EEG_EMG_mean_Acc = zeros(size(subjects_list));
EEG_EMG_std_Acc = zeros(size(subjects_list));
for i = 1:length(subjects_list)
    % EEG
    data = load([subjects_classification_results_path, folderName_EEG, ...
        '\sub-', num2str(subjects_list(i)), ...
        '\final_results_5CV_sub-', num2str(subjects_list(i)), '.mat']);
    EEG_mean_Acc(i) = data.final_struct.Mean_Accuracy;
    EEG_std_Acc(i) = data.final_struct.STD_Accuracy;
    % EMG
    data = load([subjects_classification_results_path, folderName_EMG, ...
        '\sub-', num2str(subjects_list(i)), ...
        '\final_results_5CV_sub-', num2str(subjects_list(i)), '.mat']);
    EMG_mean_Acc(i) = data.final_struct.Mean_Accuracy;
    EMG_std_Acc(i) = data.final_struct.STD_Accuracy;
    % EEG + EMG
    data = load([subjects_classification_results_path, folderName_EEG_EMG, ...
        '\sub-', num2str(subjects_list(i)), ...
        '\final_results_5CV_sub-', num2str(subjects_list(i)), '.mat']);
    EEG_EMG_mean_Acc(i) = data.final_struct.Mean_Accuracy;
    EEG_EMG_std_Acc(i) = data.final_struct.STD_Accuracy;
end


%% Calculate the random chance values
random_chance = zeros(size(subjects_list));
for i = 1:length(subjects_list)
    data = load([subjects_classification_results_path, folderName_EEG, ...
        '\sub-', num2str(subjects_list(i)), ...
        '\Subject_', num2str(subjects_list(i)),'_FlextoFlex_classes_S1S2S3_per_trial.mat']);
    C1 = sum(data.subjectTable.Label == 'S1');
    C2 = sum(data.subjectTable.Label == 'S2');
    C3 = sum(data.subjectTable.Label == 'S3');
    random_chance(i) = max([C1, C2, C3])/(C1+C2+C3);
end


%% plot the accuracies and random chance levels
figure()
hold on;
if length(subjects_list) == 14
    subjects_list([4, 6]) = [];
    random_chance([4, 6]) = [];
    EEG_mean_Acc([4, 6]) = [];
    EMG_mean_Acc([4, 6]) = [];
    EEG_EMG_mean_Acc([4, 6]) = [];
end
X = 1:length(subjects_list);
plot(X, random_chance*100, 'Marker', 'o', ...
    'LineStyle', '-', 'MarkerFaceColor', 'none', 'MarkerSize', 8,...
    'MarkerEdgeColor', 'k', 'LineWidth', 1, 'Color', 'k');
plot(X, EEG_mean_Acc*100, 'Marker', 'o', ...
    'LineStyle', '-', 'MarkerFaceColor', 'none', 'MarkerSize', 8,...
    'MarkerEdgeColor', 'r', 'LineWidth', 1, 'Color', 'r');
plot(X, EMG_mean_Acc*100, 'Marker', 'o', ...
    'LineStyle', '-', 'MarkerFaceColor', 'none', 'MarkerSize', 8,...
    'MarkerEdgeColor', 'b', 'LineWidth', 1, 'Color', 'b');
plot(X, EEG_EMG_mean_Acc*100, 'Marker', 'o', ...
    'LineStyle', '-', 'MarkerFaceColor', 'none', 'MarkerSize', 8,...
    'MarkerEdgeColor', 'g', 'LineWidth', 1, 'Color', 'g');
grid on

set(gca, 'XTick', X);
set(gca, 'XTickLabel', {'Sub-5', 'Sub-6', 'Sub-7', 'Sub-9', 'Sub-11', ...
    'Sub-12', 'Sub-13', 'Sub-14', 'Sub-15', 'Sub-16', 'Sub-17', 'Sub-18'});
ylim([0 100])
xlim([0, 13])
title('Score-Based Classification with EEG, EMG, and EEG+EMG Features')
ylabel('Accuracy [%]')
legend({'Random Chance', 'EEG Features', 'EMG Features', 'EEG+EMG Features'}, ...
    "Location", "southeast")
set(gca, 'FontSize', 12)
set(gcf, 'Position', [200, 200, 800, 400])