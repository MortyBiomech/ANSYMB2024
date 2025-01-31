clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
source_data_path = [data_path, '0_source_data\'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];

subject = 15;

%% Transform and Save the Calibrated Force (not length normalized)
% epoch_type = 'Epochs_FlextoFlex_based.mat';
% subject = 11; % [11, 12, 15, 16, 17, 18]; These 6 datasets have force sensor data
% 
% calibrated_Force = save_calibrated_force(epoch_type, epoched_data_path, subject, ...
%     source_data_path, data_path);


%% Length normalization (considering extension_start event)
% before taking the average, make sure that the extension_start events are
% aligning together as well as flextoflex start/end events. Use median
% length as a reference to decrease interpolation. 
EXP_Analysis_path = [data_path, '9_EXP_Analysis\'];
force_path = [EXP_Analysis_path, 'sub-', num2str(subject)];

data = load(fullfile(force_path, 'calibrated_Force.mat'));
name = fieldnames(data);
data = data.(name{1});

trials_info_path = [epoched_data_path, 'sub-', num2str(subject)];
trials_info = load(fullfile(trials_info_path, 'Trials_Info.mat'));
name = fieldnames(trials_info);
trials_info = trials_info.(name{1});


save_calibrated_time_normalized_force(trials_info, data, subject) 

