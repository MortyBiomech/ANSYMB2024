clc
clear

%% add paths
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master'))

%% Subject 4 
% there is no data in the computer. They might moved to the flash memory.

%% Subject 3
filepath = 'C:\Morteza\MyProjects\ANSYMB2024\data\P001\sub-P001\ses-S001\eeg';
filename = 'sub-P001_ses-S001_task-Default_run-002_eeg.xdf';

streams = load_xdf(fullfile(filepath, filename));

%% Extract EXP stream
All_Experiment_time = streams{1, 3}.time_stamps;
All_Experiment = streams{1, 3}.time_series;

%% 
figure()
plot(All_Experiment_time, All_Experiment(5, :))

% hold on
% xline(All_Experiment_time( diff(All_Experiment(6,:)) == 1 ), 'Color', 'b', 'LineStyle', '--')
% xline(All_Experiment_time( diff(All_Experiment(6,:)) == -1 ), 'Color', 'r', 'LineStyle', '--')

ylabel('Force Sensor Voltage [a.u.]')
xlabel('Time [s]')