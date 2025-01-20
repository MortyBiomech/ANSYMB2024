clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
source_data_path = [data_path, '0_source_data\'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];


%% Transform and Save the Calibrated Force 
epoch_type = 'Epochs_FlextoFlex_based.mat';
subject = 11; % [11, 14, 15, 16, 17, 18]; % These datasets have force sensor data
save_calibrated_force(epoch_type, epoched_data_path, subject, ...
    source_data_path, data_path)