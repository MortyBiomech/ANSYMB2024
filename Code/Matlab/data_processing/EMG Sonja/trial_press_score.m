
clc
clear

%%
% Author: Sonja Hanek
% This script extracts and saves the pressure levels and scorings for each
% trial

%% load data
subject_id = 7;
sessions   = 4;

[~,~, ~,~, All_Experiment,All_Experiment_time] ...
    = runs_concatenated(subject_id, sessions);

path=['/home/sonja/Documents/Studium/ANSYMB2/Analysis/data/subj_',num2str(subject_id),filesep];

load([path,'subj_',num2str(subject_id),'_Trials_encoder_events.mat'])

%%
% Beeps mark the beginning and end of each trial
start_beep = find(diff(All_Experiment(6, :)) == 1);
start_beep(1:6) = [];
finish_beep = find(diff(All_Experiment(6, :)) == -1);
finish_beep(1:6) = [];


trial_pressure_score=cell(1,length(start_beep));


for trial_nr=1:length(start_beep)-1
    trial_pressure_score{1,trial_nr}=struct( ...
        'pressure', All_Experiment(3,start_beep(trial_nr)), ...
        'score', All_Experiment(4,start_beep(trial_nr+1))); % trial_nr+1 because the scoring is given after the end of the trial
end
trial_pressure_score{1,end}=struct( ...
        'pressure', All_Experiment(3,start_beep(end)), ...
        'score', All_Experiment(4,end)); % for the last trial, the scoring is the last given score

%% save data
filename=sprintf([path,'subj_',num2str(subject_id),'_trial_pressure_score.mat']);
disp(filename)
save(filename, 'trial_pressure_score')