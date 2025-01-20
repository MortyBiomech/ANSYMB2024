clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
rawdata_path = [data_path, '0_source_data\'];


%% All signals from all sessions concatenated (it takes time!)
subject_id = 17;
output = runs_concatenated(subject_id, rawdata_path);


%% Extract data
All_EEG = output.All_EEG;
All_EEG_time = output.All_EEG_time;
All_EMG = output.All_EMG;
All_EMG_time = output.All_EMG_time;
All_Experiment = output.All_Exp;
All_Experiment_time = output.All_Exp_time;


%% load Preprocessed EEG (Cleaned with ICA)
if ~exist('ALLCOM','var')
	eeglab;
end
% load cleaned dataset (.set)
filename = ['sub-', num2str(subject_id), '_cleaned_with_ICA.set'];
filepath = [data_path, '5_single-subject-EEG-analysis\', 'sub-', ...
    num2str(subject_id), filesep];
EEG = pop_loadset('filename', filename, 'filepath', filepath);


%% Create the Trial_Info.mat file
filepath = [data_path, '6_0_Trials_Info_and_Events', filesep, 'sub-', ...
    num2str(subject_id)];
filename = ['sub-', num2str(subject_id),'_Trials_encoder_events.mat'];
load(fullfile(filepath, filename))

sessions_trial_id = [42, 78, 114, 150];  % subject 18

Trials_Info = Main_event_selection(output, ...
                                   EEG, ...
                                   Trials_encoder_events, ...
                                   subject_id, ...
                                   data_path, ...
                                   sessions_trial_id);


%% EMG sensors id 
% sensors which were used for measuring muscles activity (Delsys System)
EMG_sensor_id = [2, 3, 4, 5, 8, 9];


%% Split dataset based on events
Main_epoch_selection(output, EEG, Trials_Info, EMG_sensor_id, ...
    subject_id, data_path)



