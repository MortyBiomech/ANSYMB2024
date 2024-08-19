% Main EMG processing
clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024')); % main folder containing all codes and data
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
rawdata_path = [data_path, '0_source_data\'];


%% Colors for plotting the signals
All_colours = struct('dark_blue', [0, 0.4470, 0.7410], 'light_blue', [0.3010, 0.7450, 0.9330], ...
        'dark_orange', [0.8500, 0.3250, 0.0980], 'light_orange', [0.9290, 0.6940, 0.1250], ...
        'dark_green', [0.4660, 0.6740, 0.1880], 'light_green', [0.5960, 0.8740, 0.5410]);




%% Load data
load([data_path, '6_Trials_Info_and_Epoched_data', filesep, 'sub-', ...
    num2str(7), filesep,'Trials_Info.mat']);

load([data_path, '6_Trials_Info_and_Epoched_data', filesep, 'sub-', ...
    num2str(7), filesep,'Epochs_Flexion_based.mat']);

data = Epochs_Flexion_based;


%% Initialize a cell for saving length-normalized EMG data
EMG_Preprocessed = repmat({struct('without_outlier_removal', [], ...
    'with_outlier_removal', [])}, 1, length(Trials_Info));



%% Select the trials with the same Pressure
P1 = struct();
P1.trials = [];
P1.scores = [];

P3 = struct();
P3.trials = [];
P3.scores = [];

P6 = struct();
P6.trials = [];
P6.scores = [];

for i = 1:length(Trials_Info)
    p = Trials_Info{1, i}.General.Pressure;
    switch p
        case 1
            P1.trials(1, end+1) = i;
            P1.scores(1, end+1) = Trials_Info{1, i}.General.Score;
        case 3
            P3.trials(1, end+1) = i;
            P3.scores(1, end+1) = Trials_Info{1, i}.General.Score;
        case 6
            P6.trials(1, end+1) = i;
            P6.scores(1, end+1) = Trials_Info{1, i}.General.Score;
    end
            
end



%% Lenght Normalization - without outlier removal
%%% P1
L_P1 = zeros(1, length(P1.trials));
for i = 1:length(P1.trials)
    l = cellfun("length", data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P1(1, i) = max(l);
end
[max_L_P1, ~] = max(L_P1);

for i = 1:length(P1.trials)
    Q = length(data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P1.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P1);
        y_old = data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal(:, :, 1) = [];
end

%%% P3
L_P3 = zeros(1, length(P3.trials));
for i = 1:length(P3.trials)
    l = cellfun("length", data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P3(1, i) = max(l);
end
[max_L_P3, ~] = max(L_P3);

for i = 1:length(P3.trials)
    Q = length(data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P3.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P3);
        y_old = data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal(:, :, 1) = [];
end

%%% P6
L_P6 = zeros(1, length(P6.trials));
for i = 1:length(P6.trials)
    l = cellfun("length", data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed);
    L_P6(1, i) = max(l);
end
[max_L_P6, ~] = max(L_P6);

for i = 1:length(P6.trials)
    Q = length(data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed);
    for j = 1:Q
        x_old = data{1, P6.trials(i)}.EMG_stream.Times{1, j};
        x_new = linspace(x_old(1), x_old(end), max_L_P6);
        y_old = data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed{1, j};
        y_new = interp1(x_old', y_old', x_new', "spline");
    
        EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal(:, :, end + 1) = y_new';
    end
    EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal(:, :, 1) = [];
end


%% All the EMGs for one pressure condition together
EMG_P1 = [];
for i = 1:length(P1.trials)
    EMG_P1 = cat(3, EMG_P1, EMG_Preprocessed{1, P1.trials(i)}.without_outlier_removal);
end

EMG_P3 = [];
for i = 1:length(P3.trials)
    EMG_P3 = cat(3, EMG_P3, EMG_Preprocessed{1, P3.trials(i)}.without_outlier_removal);
end

EMG_P6 = [];
for i = 1:length(P6.trials)
    EMG_P6 = cat(3, EMG_P6, EMG_Preprocessed{1, P6.trials(i)}.without_outlier_removal);
end



%% Plot EMG signals before outlier removal
X_P1 = linspace(0, 100, size(EMG_P1, 2));
X_P3 = linspace(0, 100, size(EMG_P3, 2));
X_P6 = linspace(0, 100, size(EMG_P6, 2));

%%% P1
figure();
tiledlayout(3,4)
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;

    plot(X_P1, squeeze(EMG_P1(i, :, :))', 'Color', All_colours.light_blue, 'LineWidth', 0.5);
    plot(X_P1, mean(EMG_P1(i, :, :), 3), 'Color', All_colours.dark_blue, 'LineWidth', 2);
    % plot(X_P1, median(EMG_P1(i,:,:), 3), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
    
    hold off;
end


%%% P3
figure();
tiledlayout(3,4)
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;

    plot(X_P3, squeeze(EMG_P3(i, :, :))', 'Color', All_colours.light_orange, 'LineWidth', 0.5);
    plot(X_P3, mean(EMG_P3(i, :, :), 3), 'Color', All_colours.dark_orange, 'LineWidth', 2);
    
    hold off;
end

%%% P6
figure();
tiledlayout(3,4)
for i = 1:size(EMG_P1, 1)
    nexttile; hold on;
    
    plot(X_P6, squeeze(EMG_P6(i, :, :))', 'Color', All_colours.light_green, 'LineWidth', 0.5);
    plot(X_P6, mean(EMG_P6(i, :, :), 3), 'Color', All_colours.dark_green, 'LineWidth', 2);
    
    hold off;
end


%% Compute medians and our custom-error to find most similar signals

%%% P1
err_P1_flx = zeros(size(EMG_P1, 1), size(EMG_P1, 3));
for m = 1:size(EMG_P1, 1)
    EMG_P1_median = median(EMG_P1(m,:,:), 3);
    for i = 1:length(err_P1_flx)
        err_P1_flx(m, i) = sum((1 + abs(EMG_P1(m, :, i) - EMG_P1_median) ).^6); 
    end
end
figure()
tiledlayout(3,4)
for i = 1:size(EMG_P1, 1)
    nexttile
    bar(sort(err_P1_flx(m, :)));
end


%%% P3
err_P3_flx = zeros(size(EMG_P3, 1), size(EMG_P3, 3));
for m = 1:size(EMG_P3, 1)
    EMG_P3_median = median(EMG_P3(m,:,:), 3);
    for i = 1:length(err_P3_flx)
        err_P3_flx(m, i) = ...
            sum((1 + abs(EMG_P3(m, :, i) - EMG_P3_median) ).^6); 
    end
end
figure()
tiledlayout(3,4)
for i = 1:size(EMG_P3, 1)
    nexttile
    bar(sort(err_P3_flx(m, :)));
end


%%% P6
err_P6_flx = zeros(size(EMG_P6, 1), size(EMG_P6, 3));
for m = 1:size(EMG_P6, 1)
    EMG_P6_median = median(EMG_P6(m,:,:), 3);
    for i = 1:length(err_P6_flx)
        err_P6_flx(m, i) = ...
            sum((1 + abs(EMG_P6(m, :, i) - EMG_P6_median) ).^6); 
    end
end
figure()
tiledlayout(3,4)
for i = 1:size(EMG_P6, 1)
    nexttile
    bar(sort(err_P6_flx(m, :)));
end


%% EMG preprocessing
% Author: Sonja Hanek
% This script applies the preprocessing to the EMG data and structures the
% EMG data.
% Preprocessing consists of filtering, rectification, removal of outliers,
% amplitude normalization and event extraction.
% The resulting data is structured as follows:
% EMG_data --- muscle names 
%          |-- time stamps for the raw data
%          |-- raw data of the EMG sensors
%          |-- preprocessed data of each trial  
%                   |--- trial_nr
%                   |-- pressure level of the trial
%                   |-- scoring of the trial
%                   |-- flexions            
%                   |        |-- preprocessed data of individual flexions
%                   |            within the trial
%                   |-- extensions            
%                   |        |-- preprocessed data of individual extensions
%                   |            within the trial
%                   |-- epoch packs
%                   |        |-- preprocessed data of individual epoch
%                   |        |    packs, i.e. flexion-extension events
%                   |        |-- index of the start of the flexion within
%                   |             the epoch pack
%                   |-- epoch timing
%                           |-- time stamps of beginnings of flexions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of beginnings of flexions of
%                           |    the trial within EMG data
%                           |-- time stamps of endings of flexions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of endings of flexions of
%                           |    the trial within EMG data
%                           |-- time stamps of beginnings of extensions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of beginnings of extensions of
%                           |    the trial within EMG data
%                           |-- time stamps of endings of extensions of
%                           |    the trial in the timestamps of EMG data
%                           |-- index of endings of extensions of
%                               the trial within EMG data
%                                               
% Authors: 
% Sonja Hanek (abc@gmail.com) 
% Morteza Khosrotabar (mkhosrotabar@gmail.com)

%% Load the data
% Before running this section you need to modify the paths in the 
% xdf_load_matlab.m function based on your directory
% you need to run trial_press_score first, then find_epochs.m  !

% All signals from all sessions concatenated (it takes time!)
subject_id = 7;

output = runs_concatenated(subject_id, rawdata_path);
All_EMG = output.All_EMG;
All_EMG_time = output.All_EMG_time;

%%
% % for faster loading, save data
%  load('All_EMG_raw.mat')
%  load('All_EMG_time.mat')

load([path, 'subj_', num2str(subject_id),'_epoch_timestamps.mat'])
load([path,'subj_',num2str(subject_id),'_trial_pressure_score.mat'])



%% filter data

% measurement at 2 kHz 
fs=2e3;
%first muscle is just zeros, since we used the sensors 2-11
EMG_filtered=zeros(10,length(All_EMG));
for muscle=2:11

    % get data
    Raw=All_EMG(muscle,:);



    %bandpass filter
    
    %BANDPASS FILTER 
    %A bandpass filter is used to filter the EMG.
    
    fn=fs/2;    %Hz Nyquist Frequency is 1/2 Sampling Frequency
    low_freq=20;    %Hz
    high_freq=450;  %Hz
    order=4;
    [b,a]=butter(order,([low_freq high_freq]/fn)); % determine filter coefficients
    EMG_r_filt=filtfilt(b,a,Raw); %filtfilt provides zero-lag

    % rectification
    rect=abs(EMG_r_filt);

    %LOWPASS
    fc = 4; % Cut-off frequency (Hz)
    order = 2; % Filter order
    [b,a] = butter(order,fc/fn);
    EMG_r_filt=filtfilt(b,a,rect); %filtfilt provides zero-lag

    EMG_filtered(muscle-1,:)=EMG_r_filt;
end


% plot(All_EMG_time, EMG_filtered(1, :))
% zoom on
% hold on
% plot(All_EMG_time, EMG_filtered(3, :))

%% separate events
% find and save timestamps and indices of the beginning and end of each
% flexion and extension

Emg_trial_epochs=cell(1,length(epochs));

for trial_nr=1:length(epochs)
    trial=epochs{1,trial_nr};
    nr_flexions=length(trial.flexion_start);
    emg_flexion_start_time=zeros(1,nr_flexions);
    emg_flexion_start_index=zeros(1,nr_flexions);

    emg_flexion_end_time=zeros(1,nr_flexions);
    emg_flexion_end_index=zeros(1,nr_flexions);

    nr_extensions=length(trial.extension_start);
    emg_extension_start_time=zeros(1,nr_extensions);
    emg_extension_start_index=zeros(1,nr_extensions);

    emg_extension_end_time=zeros(1,nr_extensions);
    emg_extension_end_index=zeros(1,nr_extensions);

    for flexion=1:nr_flexions
        [~,start_index]=min(abs(All_EMG_time-trial.flexion_start(flexion)));
        [~, end_index]=min(abs(All_EMG_time-trial.flexion_end(flexion)));
        emg_flexion_start_time(flexion)=All_EMG_time(start_index);
        emg_flexion_start_index(flexion)=start_index;
        emg_flexion_end_time(flexion)=All_EMG_time(end_index);
        emg_flexion_end_index(flexion)=end_index;

    end
    for extension=1:nr_extensions
        [~,start_index]=min(abs(All_EMG_time-trial.extension_start(extension)));
        [~, end_index]=min(abs(All_EMG_time-trial.extension_end(extension)));
        emg_extension_start_time(extension)=All_EMG_time(start_index);
        emg_extension_start_index(extension)=start_index;
        emg_extension_end_time(extension)=All_EMG_time(end_index);
        emg_extension_end_index(extension)=end_index;

    end
    Emg_trial_epochs{trial_nr}=struct('emg_flexion_start_time', emg_flexion_start_time, ...
        'emg_flexion_start_index',emg_flexion_start_index, ...
        'emg_flexion_end_time',emg_flexion_end_time, ...
        'emg_flexion_end_index',emg_flexion_end_index, ...
        'emg_extension_start_time',emg_extension_start_time, ...
        'emg_extension_start_index',emg_extension_start_index, ...
        'emg_extension_end_time',emg_extension_end_time, ...
        'emg_extension_end_index',emg_extension_end_index);
end



%% Remove outliers
% For each muscle and each pressure level, epoch packs are identified.
% Epoch packs consist of one flexion and one extension, i.e. they last from
% the beginning of one flexion to the beginning of the next flexion within 
% the trial.
% The epoch packs are time normalized, i.e. interpolated to the sample
% count of the longest epoch pack within the muscle-pressure-level
% category. The error to the median of the epoch packs in the category is 
% calculated:
% The sum of [absolute distances of each sample to the median plus one, to
% the power of 6]
% The 300 epoch packs with the lowest error are kept, the EMG of all other
% epochs is set to NaN.
% This ensures that the evaluated data is of epochs where the muscle
% activity expresses a somewhat typical behaviour and outliers are removed.
% This does not mean that we don't exclude epochs that are actually no
% outliers.
% 

trial_epoch_packs=cell(1,length(Emg_trial_epochs));
for muscle=1:10
    for pressure=[1,3,6]
        epoch_pack_indices=[];
        for trial_nr=1:length(Emg_trial_epochs)
           if trial_pressure_score{1,trial_nr}.pressure==pressure
                flexion_starts=Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index;
                epoch_pack=cell(1,length(flexion_starts)-1);
                for flexion = 1:length(flexion_starts)-1
                    epoch_pack_indices=[epoch_pack_indices;
                        flexion_starts(flexion),flexion_starts(flexion+1),trial_nr];
                end
           end
        end
        epoch_pack_indices(:,4)=epoch_pack_indices(:,2)-epoch_pack_indices(:,1);
        time_norm_length=max(epoch_pack_indices(:,4));
        time_norm_epochs=zeros(length(epoch_pack_indices),time_norm_length);
        for pack_nr=1:length(epoch_pack_indices)
            pack=epoch_pack_indices(pack_nr,:);
            epoch_pack=EMG_filtered(muscle,pack(1):pack(2));
            epoch_time_norm=interp1(epoch_pack,(1:time_norm_length)/(time_norm_length/length(epoch_pack)));
            time_norm_epochs(pack_nr,:)=epoch_time_norm;
        end
        median_pack=median(time_norm_epochs,1,"omitnan");
        errors=abs(time_norm_epochs-median_pack);
        errors=(errors+1).^6;
        errors=sum(errors,2, "omitnan");
        outlier_detection_matrix=[epoch_pack_indices, errors,time_norm_epochs];
        outlier_detection_matrix=sortrows(outlier_detection_matrix,5);
%        % uncomment to examine plots of errors
%         figure;plot(outlier_detection_matrix(:,5))
%         title(join(['muscle: ',num2str(muscle),' pressure: ', num2str(pressure)]))
        for outs=301:length(errors) % keep only the 300 epoch packs (around half) that are most consistent
            EMG_filtered(muscle,outlier_detection_matrix(outs,1):outlier_detection_matrix(outs,2))=NaN;
        end
    end
end

%% amplitude normalization
% For each muscle, the maximum value of the muscle activity in the first 
% trial with a pressure level of 1 is taken as reference for amplitude 
% normalization

trial_press=[trial_pressure_score{:}];
trial_press=[trial_press(:).pressure];
ref_trial_nr=find(trial_press==1,1);
ref_trial=Emg_trial_epochs{1,ref_trial_nr};
trial_start_index=min(ref_trial.emg_flexion_start_index(1),ref_trial.emg_extension_start_index(1));
trial_end_index=max(ref_trial.emg_flexion_end_index(end), ref_trial.emg_extension_end_index(end));
ref_val=max(EMG_filtered(:,trial_start_index:trial_end_index),[],2,"omitnan");
for muscle=1:10
    EMG_filtered(muscle,:)=EMG_filtered(muscle,:)/ref_val(muscle);
end

%% structure data

EMG_filtered_epochs=cell(1,length(epochs));
EMG_filtered_epoch_packs=cell(1,length(epoch_pack_indices));
for trial_nr=1:length(epochs)

    flexion=cell(1,length(Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index));
    for flex_nr=1:length(Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index)
        flexion{1,flex_nr}=EMG_filtered(:, ...
            Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index(flex_nr): ...
            Emg_trial_epochs{1,trial_nr}.emg_flexion_end_index(flex_nr));
    end

    extension=cell(1,length(Emg_trial_epochs{1,trial_nr}.emg_extension_start_index));
    
    for ext_nr=1:length(Emg_trial_epochs{1,trial_nr}.emg_extension_start_index)
        extension{1,ext_nr}=EMG_filtered(:, ...
            Emg_trial_epochs{1,trial_nr}.emg_extension_start_index(ext_nr): ...
            Emg_trial_epochs{1,trial_nr}.emg_extension_end_index(ext_nr));
    end

    flexion_starts=Emg_trial_epochs{1,trial_nr}.emg_flexion_start_index;
    extension_starts=Emg_trial_epochs{1,trial_nr}.emg_extension_start_index;
    epoch_pack=cell(1,length(flexion_starts)-1);
    for pack_nr = 1:length(flexion_starts)-1
        ext_start=extension_starts(flexion_starts(pack_nr)< extension_starts & ...
            extension_starts<flexion_starts(pack_nr+1))-flexion_starts(pack_nr);
        epoch_pack{1,pack_nr}=struct(...
            'filtered_data', EMG_filtered(:,...
            flexion_starts(pack_nr):flexion_starts(pack_nr+1)),...
            'relative_extension_start_index', ext_start);
    end

        
    EMG_filtered_epochs{trial_nr}=struct( ...
        'trial_nr', trial_nr, ...
        'pressure', epochs{1,trial_nr}.pressure, ...
        'score',epochs{1,trial_nr}.score, ...
        'flexions', {flexion}, ...    
        'extensions',{extension}, ...
        'epoch_packs', {epoch_pack},...
        'epoch_timing', Emg_trial_epochs{1,trial_nr});
end

names={'Vastus_med_R', 'Rectus_femoris_R', 'Gastrocnemius_R', 'Biceps_femoris_R', ...
    'Vastus_med_L', 'Rectus_femoris_L', 'Gastrocnemius_L', 'Biceps_femoris_L', ...
    'Trapezius_R', 'Trapezius_L'};

EMG_data=struct( ...
    'Muscle_names', {names},...
    'Timestamps', All_EMG_time, ...
    'Raw_data', All_EMG, ...
    'Epochs_preprocessed', {EMG_filtered_epochs});

%% save data
filename=sprintf([path,'subj_',num2str(subject_id),'_emg.mat']);
disp(filename)
save(filename, 'EMG_data')