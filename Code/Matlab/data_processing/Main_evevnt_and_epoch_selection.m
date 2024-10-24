clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
rawdata_path = [data_path, '0_source_data\'];


%% All signals from all sessions concatenated (it takes time!)
subject_id = 10;
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


%% Finding peaks of Encoder data for using in App to remove undesired peaks
% [start_beep, ...
%     finish_beep, XLimits, ...
%     pks_high_peaks, locs_high_peaks, ...
%     pks_low_peaks, locs_low_peaks, ...
%     trial_pks_high_peaks, trial_locs_high_peaks, ...
%     trial_pks_low_peaks, trial_locs_low_peaks, ...
%     N_Trials, ...
%     Trials_encoder_events] = find_peaks_and_select_events(output);


%% Matlab App for deselecting undesired peaks of Encoder data
% find_flexion_extension_events
filepath = [data_path, '6_0_Trials_Info_and_Events', filesep, 'sub-', ...
    num2str(subject_id)];
filename = ['subj_', num2str(subject_id),'_Trials_encoder_events.mat'];
load(fullfile(filepath, filename))


% %% Matlab App for marking the bad trials to exclude from post-processing
% % bad_trials_EEG_based = [];
% % mark_bad_trials_of_EEG_data
% % %%
% % filepath = [data_path, '6_0_Trials_Info_and_Events\', 'sub-', ...
% %     num2str(subject_id), filesep];
% % filename = 'bad_trials_EEG_based.mat';
% % save(fullfile(filepath, filename), 'bad_trials_EEG_based')
% filepath = [data_path, '6_0_Trials_Info_and_Events', filesep, 'sub-', ...
%     num2str(subject_id)];
% filename = 'bad_trials_EEG_based.mat';
% load(fullfile(filepath, filename))

%% Find all events based on entire trials, flexions, extension, and
% flextoflex  epochs and store in a big structure (Trials_Info.mat)

% for subject 10:
sessions_trial_id = [57, 108, 159, 210]; 

Trials_Info = Main_event_selection(output, ...
                                   EEG, ...
                                   Trials_encoder_events, ...
                                   subject_id, ...
                                   data_path, ...
                                   sessions_trial_id);


%% EMG sensors id 
% sensors which were used for measuring muscles activity (Delsys System)
EMG_sensor_id = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];


%% Split dataset based on events
% Note:
%      Epoch selection was performed on subject 7 and datasets are stored.

Main_epoch_selection(output, EEG, Trials_Info, EMG_sensor_id, ...
    subject_id, data_path)



