
clc
clear

%%
% Author: Sonja Hanek
% Identify and save the timestamps of the beginning and end of epochs, i.e.
% flexions and extensions for eacht trial. Also saves the pressure level
% and scoring of the respective trials.

%% load data
path='/home/sonja/Documents/Studium/ANSYMB2/Analysis/data/';
subj=7;
load([path,'subj_',num2str(subj),'/subj_',num2str(subj),'_Trials_encoder_events.mat'])
load([path,'subj_',num2str(subj),'/subj_',num2str(subj),'_trial_pressure_score.mat'])
%% find epochs

epochs=cell(1,length(Trials_encoder_events)-1); %last trial in encoder_events is empty!!

for trial_nr = 1:length(epochs)

    high_time=Trials_encoder_events{1,trial_nr}.high_peaks.time;
    low_time=Trials_encoder_events{trial_nr}.low_peaks.time;

    if high_time(1)<low_time(1)
        first_peak_high=true;
    else
        first_peak_high=false;
    end
    if high_time(end)<low_time(end)
        nr_flexions=length(high_time);
        nr_extensions=length(low_time)-1;

    else
        nr_flexions=length(high_time)-1;
        nr_extensions=length(low_time);

    end

    flexion_start =zeros(1,nr_flexions);
    extension_start=zeros(1,nr_extensions);
    
    flexion_end =zeros(1,nr_flexions);
    extension_end=zeros(1,nr_extensions);
    

    for flexion_start_index=1:nr_flexions
        if(first_peak_high)
            flexion_end_index=flexion_start_index;
        else
            flexion_end_index=flexion_start_index+1;
        end
        flexion_start(flexion_start_index)= high_time(flexion_start_index);
        flexion_end(flexion_start_index)=low_time(flexion_end_index);

    end

    for extension_start_index=1:nr_extensions
        if(first_peak_high)
            extension_end_index=extension_start_index+1;
        else
            extension_end_index=extension_start_index;
        end
        extension_start(extension_start_index)= low_time(extension_start_index);
        extension_end(extension_start_index)= high_time(extension_end_index);

    end
    epochs{trial_nr}=struct( ...
        'flexion_start',flexion_start, ...
        'flexion_end',flexion_end, ...
        'extension_start',extension_start, ...
        'extension_end',extension_end, ...
        'pressure', trial_pressure_score{1,trial_nr}.pressure, ...
        'score',trial_pressure_score{1,trial_nr}.score);
end

%% check if result makes sense
for trial=epochs
    flexion_start=trial{1,1}.flexion_start;
    flexion_end=trial{1,1}.flexion_end;
    if length(flexion_start)~=length(flexion_end)
        disp('ERROR: flexion start and end have not the same length')
    end
    for i=1:length(flexion_start)-1
        if (flexion_start(i)>=flexion_end(i))
            disp('ERROR: flexion')
        end
    end
    extension_start=trial{1,1}.extension_start;
    extension_end=trial{1,1}.extension_end;
    if length(extension_start)~=length(extension_end)
        disp('ERROR: extension start and end have not the same length')
    end
    for i=1:length(extension_start)-1
        if (extension_start(i)>=extension_end(i))
            disp('ERROR:extension')
        end
    end
end

%% save epoch timestamps
filename=sprintf([path,'/subj_',num2str(subj),'/subj_',num2str(subj),'_epoch_timestamps.mat']);
disp(filename)
save(filename, 'epochs')



