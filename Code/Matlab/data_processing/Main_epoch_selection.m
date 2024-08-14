clc
clear

%% Add and Define Necessary Paths
addpath(genpath('C:\Morteza\Analysis\ANSYMB2024')); % main folder containing all codes and data
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master')); % required for loading the XDF files.

% Change path to the directory on your PC which raw XDF files are stored:
data_path = 'C:\Morteza\Analysis\ANSYMB2024\data\';
rawdata_path = [data_path, '0_source_data\'];

%% All signals from all sessions concatenated (it takes time!)
subject_id = 7;
output = runs_concatenated(subject_id, rawdata_path);

%% Extract data
All_EEG = output.All_EEG;
All_EEG_time = output.All_EEG_time;
All_EMG = output.All_EMG;
All_EMG_time = output.All_EMG_time;
All_Experiment = output.All_Exp;
All_Experiment_time = output.All_Exp_time;


%% load Preprocessed EEG
if ~exist('ALLCOM','var')
	eeglab;
end
% load cleaned dataset (.set)
filename = ['sub-', num2str(subject_id), '_cleaned_with_ICA.set'];
filepath = [data_path, '5_single-subject-EEG-analysis\', 'sub-', ...
    num2str(subject_id), filesep];
EEG = pop_loadset('filename', filename, 'filepath', filepath);


%% load Trials_Info
Trials_Info_path = [data_path, '6_Trials_Info_and_Epoched_data\', ...
    'sub-', num2str(subject_id)];
load(fullfile(Trials_Info_path, 'Trials_Info.mat'));

%% EMG sensors id 
% sensors which were used for measuring muscles activity
EMG_sensor_id = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

%% Trials_Based epoch selection (takes time to save the mat file!)
clc;
Epochs_Trial_based = main_epoch_selection_TrialsBased(output, ...
                                                      EEG, ...
                                                      Trials_Info, ...
                                                      EMG_sensor_id, ...
                                                      subject_id, ...
                                                      data_path);



%% Cut EEG signal (sensor and source level)





%% plot the scores vs pressures

% addpath('C:\Morteza\Analysis\ANSYMB2024\Code\data_processing\Plot functions\Violin Plots\distributionPlot')
% figure;
% distributionPlot({P1.scores, P3.scores, P6.scores}, 'showMM', 2);
% set(gca, 'XTickLabel', {'1 bar', '3 bar', '6 bar'});
% xlabel('PAM Pressure');
% ylabel('Scores');
% title('Scores Comparison Across PAM Pressures');


% Calculate means and standard errors
means = [mean(P1.scores), mean(P3.scores), mean(P6.scores)];
stdErrors = [std(P1.scores)/sqrt(length(P1.scores)), ...
             std(P3.scores)/sqrt(length(P3.scores)), ...
             std(P6.scores)/sqrt(length(P6.scores))];
stdDevs = [std(P1.scores), std(P3.scores), std(P6.scores)];

figure;
bar(means);
hold on;
h = errorbar(1:3, means, stdDevs, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
h.CapSize = 20; % Increase the size of the horizontal caps
set(gca, 'XTick', 1:3, 'XTickLabel', {'1', '3', '6'});
xlabel('PAM Pressure [bar]');
ylabel('Scores');
title('Scores Comparison Across PAM Pressures');
set(gca, 'FontSize', 10); 


% Perform pairwise unpaired t-tests
[~, p12] = ttest2(P1.scores, P3.scores);
[~, p13] = ttest2(P1.scores, P6.scores);
[~, p23] = ttest2(P3.scores, P6.scores);

% Bonferroni correction
alpha = 0.05;
alpha_corrected = alpha / 3;

% Add significance markers
yMax = max(means + stdDevs) + 1;
if p12 < alpha_corrected
    plot([1 2], [yMax yMax], '-k', 'LineWidth', 1);
    text(1.5, yMax, '$\ast$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14, 'BackgroundColor', 'white', 'Rotation', 90, 'Margin', 0.1,'Interpreter', 'latex');
    yMax = yMax + 0.5;
end
if p13 < alpha_corrected
    plot([1 3], [yMax yMax], '-k', 'LineWidth', 1);
    text(2, yMax, '$\ast$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14, 'BackgroundColor', 'white', 'Rotation', 90, 'Margin', 0.1,'Interpreter', 'latex');
    yMax = yMax + 0.5;
end
if p23 < alpha_corrected
    plot([2 3], [yMax yMax], '-k', 'LineWidth', 1);
    text(2.5, yMax, '$\ast$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14, 'BackgroundColor', 'white', 'Rotation', 90, 'Margin', 0.1,'Interpreter', 'latex');
end

hold off;



% % Create a jittered scatter plot
% means = [mean(P1.scores), mean(P3.scores), mean(P6.scores)];
% stdErrors = [std(P1.scores)/sqrt(length(P1.scores)), ...
%              std(P3.scores)/sqrt(length(P3.scores)), ...
%              std(P6.scores)/sqrt(length(P6.scores))];
% jitterAmount = 0.1; % Adjust as necessary
% figure;
% hold on;
% scatter(1 + jitterAmount*(rand(size(P1.scores))-0.5), P1.scores, 'r');
% scatter(2 + jitterAmount*(rand(size(P3.scores))-0.5), P3.scores, 'g');
% scatter(3 + jitterAmount*(rand(size(P6.scores))-0.5), P6.scores, 'b');
% % Plot mean and standard deviation
% errorbar(1:3, means, stdErrors, 'k', 'LineStyle', 'none');
% set(gca, 'XTick', 1:3, 'XTickLabel', {'1 bar', '2 bar', '3 bar'});
% xlabel('PAM Pressure');
% ylabel('Scores');
% title('Scores Comparison Across PAM Pressures');
% hold off;



%% Cut and Store EEG sensor and source activity based on trials pressure

% Predefinen the nested stucture

length_structure = struct('not_normalized', [], 'length_normalized', []);
sensor_source_structure = struct('channels', length_structure, 'sources', length_structure, 'EEG_times', []);
freq_sensor_source_structure = struct('channels', [], 'sources', []);
time_freq_struct = struct('Time_Domain', sensor_source_structure, 'Frequency_Domain', freq_sensor_source_structure);
EEG_epoched = struct('EEG_P1_trials', time_freq_struct, 'EEG_P3_trials', time_freq_struct, 'EEG_P6_trials', time_freq_struct, ...
    'EEG_P1_flexion', time_freq_struct, 'EEG_P3_flexion', time_freq_struct, 'EEG_P6_flexion', time_freq_struct, ...
    'EEG_P1_extension', time_freq_struct, 'EEG_P3_extension', time_freq_struct, 'EEG_P6_extension', time_freq_struct, ...
    'EEG_P1_flextoflex', time_freq_struct, 'EEG_P3_flextoflex', time_freq_struct, 'EEG_P6_flextoflex', time_freq_struct);

%% loop over the trials to fill the big structure

% EEG_P1_trials
L_P1_trials = zeros(1, numel(P1.trials));
for i = 1:numel(P1.trials)
    trial_i = P1.trials(i);
    start_indx = Trials_Info{1, P1.trials(i)}.Trial_start_indx;
    end_indx   = Trials_Info{1, P1.trials(i)}.Trial_end_indx;
    t = EEG.times(start_indx:end_indx);
    EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times{1, 1, i} = t;
    L_P1_trials(1, i) = length(t);

    EEG_epoched.EEG_P1_trials.Time_Domain.channels.not_normalized{1, 1, i} = ...
        channel_data(:,start_indx:end_indx);

    EEG_epoched.EEG_P1_trials.Time_Domain.sources.not_normalized{1, 1, i} = ...
        source_data(:,start_indx:end_indx);

end

% % check with a plot
% figure(); plot(L_P1_trials)
% set(gca, 'YLim', [10000 12000])
% ylabel(P1_trials_length)
% 
% figure(); hold on
% for i = 1:length(L_P1_trials)
%     plot([EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times{1, 1, i}(1), ...
%         EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times{1, 1, i}(end)], ...
%         [1 1], 'b');
% end
% set(gca, 'YLim', [0.8 1.2], 'YTick', 1, 'YTicklabel', '')
% hold off


%% EEG_P3_trials
L_P3_trials = zeros(1, numel(P3.trials));
for i = 1:numel(P3.trials)
    trial_i = P3.trials(i);
    start_indx = Trials_Info{1, P3.trials(i)}.Trial_start_indx;
    end_indx   = Trials_Info{1, P3.trials(i)}.Trial_end_indx;
    t = EEG.times(start_indx:end_indx);
    EEG_epoched.EEG_P3_trials.Time_Domain.EEG_times{1, 1, i} = t;
    L_P3_trials(1, i) = length(t);

    EEG_epoched.EEG_P3_trials.Time_Domain.channels.not_normalized{1, 1, i} = ...
        channel_data(:,start_indx:end_indx);

    EEG_epoched.EEG_P3_trials.Time_Domain.sources.not_normalized{1, 1, i} = ...
        source_data(:,start_indx:end_indx);

end


%% EEG_P6_trials
L_P6_trials = zeros(1, numel(P6.trials));
for i = 1:numel(P6.trials)
    trial_i = P6.trials(i);
    start_indx = Trials_Info{1, P6.trials(i)}.Trial_start_indx;
    end_indx   = Trials_Info{1, P6.trials(i)}.Trial_end_indx;
    t = EEG.times(start_indx:end_indx);
    EEG_epoched.EEG_P6_trials.Time_Domain.EEG_times{1, 1, i} = t;
    L_P6_trials(1, i) = length(t);

    EEG_epoched.EEG_P6_trials.Time_Domain.channels.not_normalized{1, 1, i} = ...
        channel_data(:,start_indx:end_indx);

    EEG_epoched.EEG_P6_trials.Time_Domain.sources.not_normalized{1, 1, i} = ...
        source_data(:,start_indx:end_indx);

end



%% EEG_P1_flexion

flexion_epochs_number = zeros(1, numel(P1.trials));
for i = 1:numel(P1.trials)
    flexion_epochs_number(1, i) = length(Trials_Info{1, P1.trials(i)}.flexion_start_indx);
end

L_P1_flexion = [];
for i = 1:numel(P1.trials)
    trial_i = P1.trials(i);
    for j = 1:numel(Trials_Info{1, P1.trials(i)}.flexion_start_indx)
        start_indx = Trials_Info{1, P1.trials(i)}.flexion_start_indx(j);
        end_indx   = Trials_Info{1, P1.trials(i)}.flexion_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P1_flexion.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P1_flexion(1, end+1) = length(t);

        EEG_epoched.EEG_P1_flexion.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P1_flexion.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P1_flexion.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P1_flexion.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P1_flexion.Time_Domain.sources.not_normalized(:, :, 1) = [];



%% EEG_P3_flexion

flexion_epochs_number = zeros(1, numel(P3.trials));
for i = 1:numel(P3.trials)
    flexion_epochs_number(1, i) = length(Trials_Info{1, P3.trials(i)}.flexion_start_indx);
end

L_P3_flexion = [];
for i = 1:numel(P3.trials)
    trial_i = P3.trials(i);
    for j = 1:numel(Trials_Info{1, P3.trials(i)}.flexion_start_indx)
        start_indx = Trials_Info{1, P3.trials(i)}.flexion_start_indx(j);
        end_indx   = Trials_Info{1, P3.trials(i)}.flexion_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P3_flexion.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P3_flexion(1, end+1) = length(t);

        EEG_epoched.EEG_P3_flexion.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P3_flexion.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P3_flexion.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P3_flexion.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P3_flexion.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P6_flexion

flexion_epochs_number = zeros(1, numel(P6.trials));
for i = 1:numel(P6.trials)
    flexion_epochs_number(1, i) = length(Trials_Info{1, P6.trials(i)}.flexion_start_indx);
end

L_P6_flexion = [];
for i = 1:numel(P6.trials)
    trial_i = P6.trials(i);
    for j = 1:numel(Trials_Info{1, P6.trials(i)}.flexion_start_indx)
        start_indx = Trials_Info{1, P6.trials(i)}.flexion_start_indx(j);
        end_indx   = Trials_Info{1, P6.trials(i)}.flexion_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P6_flexion.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P6_flexion(1, end+1) = length(t);

        EEG_epoched.EEG_P6_flexion.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P6_flexion.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P6_flexion.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P6_flexion.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P6_flexion.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P1_extension

extension_epochs_number = zeros(1, numel(P1.trials));
for i = 1:numel(P1.trials)
    extension_epochs_number(1, i) = length(Trials_Info{1, P1.trials(i)}.extension_start_indx);
end

L_P1_extension = [];
for i = 1:numel(P1.trials)
    trial_i = P1.trials(i);
    for j = 1:numel(Trials_Info{1, P1.trials(i)}.extension_start_indx)
        start_indx = Trials_Info{1, P1.trials(i)}.extension_start_indx(j);
        end_indx   = Trials_Info{1, P1.trials(i)}.extension_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P1_extension.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P1_extension(1, end+1) = length(t);

        EEG_epoched.EEG_P1_extension.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P1_extension.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P1_extension.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P1_extension.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P1_extension.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P3_extension

extension_epochs_number = zeros(1, numel(P3.trials));
for i = 1:numel(P3.trials)
    extension_epochs_number(1, i) = length(Trials_Info{1, P3.trials(i)}.extension_start_indx);
end

L_P3_extension = [];
for i = 1:numel(P3.trials)
    trial_i = P3.trials(i);
    for j = 1:numel(Trials_Info{1, P3.trials(i)}.extension_start_indx)
        start_indx = Trials_Info{1, P3.trials(i)}.extension_start_indx(j);
        end_indx   = Trials_Info{1, P3.trials(i)}.extension_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P3_extension.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P3_extension(1, end+1) = length(t);

        EEG_epoched.EEG_P3_extension.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P3_extension.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P3_extension.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P3_extension.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P3_extension.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P6_extension

extension_epochs_number = zeros(1, numel(P6.trials));
for i = 1:numel(P6.trials)
    extension_epochs_number(1, i) = length(Trials_Info{1, P6.trials(i)}.extension_start_indx);
end

L_P6_extension = [];
for i = 1:numel(P6.trials)
    trial_i = P6.trials(i);
    for j = 1:numel(Trials_Info{1, P6.trials(i)}.extension_start_indx)
        start_indx = Trials_Info{1, P6.trials(i)}.extension_start_indx(j);
        end_indx   = Trials_Info{1, P6.trials(i)}.extension_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P6_extension.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P6_extension(1, end+1) = length(t);

        EEG_epoched.EEG_P6_extension.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P6_extension.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P6_extension.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P6_extension.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P6_extension.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P1_flextoflex

flextoflex_epochs_number = zeros(1, numel(P1.trials));
for i = 1:numel(P1.trials)
    flextoflex_epochs_number(1, i) = length(Trials_Info{1, P1.trials(i)}.flextoflex_start_indx);
end

L_P1_flextoflex = [];
for i = 1:numel(P1.trials)
    trial_i = P1.trials(i);
    for j = 1:numel(Trials_Info{1, P1.trials(i)}.flextoflex_start_indx)
        start_indx = Trials_Info{1, P1.trials(i)}.flextoflex_start_indx(j);
        end_indx   = Trials_Info{1, P1.trials(i)}.flextoflex_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P1_flextoflex.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P1_flextoflex(1, end+1) = length(t);

        EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P1_flextoflex.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P3_flextoflex

flextoflex_epochs_number = zeros(1, numel(P3.trials));
for i = 1:numel(P3.trials)
    flextoflex_epochs_number(1, i) = length(Trials_Info{1, P3.trials(i)}.flextoflex_start_indx);
end

L_P3_flextoflex = [];
for i = 1:numel(P3.trials)
    trial_i = P3.trials(i);
    for j = 1:numel(Trials_Info{1, P3.trials(i)}.flextoflex_start_indx)
        start_indx = Trials_Info{1, P3.trials(i)}.flextoflex_start_indx(j);
        end_indx   = Trials_Info{1, P3.trials(i)}.flextoflex_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P3_flextoflex.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P3_flextoflex(1, end+1) = length(t);

        EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P3_flextoflex.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% EEG_P6_flextoflex

flextoflex_epochs_number = zeros(1, numel(P6.trials));
for i = 1:numel(P6.trials)
    flextoflex_epochs_number(1, i) = length(Trials_Info{1, P6.trials(i)}.flextoflex_start_indx);
end

L_P6_flextoflex = [];
for i = 1:numel(P6.trials)
    trial_i = P6.trials(i);
    for j = 1:numel(Trials_Info{1, P6.trials(i)}.flextoflex_start_indx)
        start_indx = Trials_Info{1, P6.trials(i)}.flextoflex_start_indx(j);
        end_indx   = Trials_Info{1, P6.trials(i)}.flextoflex_end_indx(j);
        t = EEG.times(start_indx:end_indx);
        EEG_epoched.EEG_P6_flextoflex.Time_Domain.EEG_times{1, 1, end+1} = t;
        L_P6_flextoflex(1, end+1) = length(t);

        EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.not_normalized{1, 1, end+1} = ...
            channel_data(:,start_indx:end_indx);
        EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.not_normalized{1, 1, end+1} = ...
            source_data(:,start_indx:end_indx);
    end
end

EEG_epoched.EEG_P6_flextoflex.Time_Domain.EEG_times(:, :, 1) = [];
EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.not_normalized(:, :, 1) = [];
EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.not_normalized(:, :, 1) = [];




%% Now let's do the length normalization
% Each condition is normalized to the maximum lenght of its group. For
% example: P1_flexion epochs are normalized to the maximum lenght of
% P1_flexion.


%% EEG_P1_trials_lenght_normalization
[max_L_P1_trials, max_L_P1_trials_indx] = max(L_P1_trials);

for i = 1:size(EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times, 3)
    if i ~= max_L_P1_trials_indx
        x_old = EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P1_trials);

        % channels data
        data = EEG_epoched.EEG_P1_trials.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_trials.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_trials.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P1_trials.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_trials.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_trials.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P3_trials_lenght_normalization
[max_L_P3_trials, max_L_P3_trials_indx] = max(L_P3_trials);

for i = 1:size(EEG_epoched.EEG_P3_trials.Time_Domain.EEG_times, 3)
    if i ~= max_L_P3_trials_indx
        x_old = EEG_epoched.EEG_P3_trials.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P3_trials);

        % channels data
        data = EEG_epoched.EEG_P3_trials.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_trials.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_trials.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P3_trials.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_trials.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_trials.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P6_trials_lenght_normalization
[max_L_P6_trials, max_L_P6_trials_indx] = max(L_P6_trials);

for i = 1:size(EEG_epoched.EEG_P6_trials.Time_Domain.EEG_times, 3)
    if i ~= max_L_P6_trials_indx
        x_old = EEG_epoched.EEG_P6_trials.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P6_trials);

        % channels data
        data = EEG_epoched.EEG_P6_trials.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_trials.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_trials.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P6_trials.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_trials.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_trials.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P1_flexion_lenght_normalization
[max_L_P1_flexion, max_L_P1_flexion_indx] = max(L_P1_flexion);

for i = 1:size(EEG_epoched.EEG_P1_flexion.Time_Domain.EEG_times, 3)
    if i ~= max_L_P1_flexion_indx
        x_old = EEG_epoched.EEG_P1_flexion.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P1_flexion);

        % channels data
        data = EEG_epoched.EEG_P1_flexion.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_flexion.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_flexion.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P1_flexion.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_flexion.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_flexion.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P3_flexion_lenght_normalization
[max_L_P3_flexion, max_L_P3_flexion_indx] = max(L_P3_flexion);

for i = 1:size(EEG_epoched.EEG_P3_flexion.Time_Domain.EEG_times, 3)
    if i ~= max_L_P3_flexion_indx
        x_old = EEG_epoched.EEG_P3_flexion.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P3_flexion);

        % channels data
        data = EEG_epoched.EEG_P3_flexion.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_flexion.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_flexion.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P3_flexion.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_flexion.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_flexion.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P6_flexion_lenght_normalization
[max_L_P6_flexion, max_L_P6_flexion_indx] = max(L_P6_flexion);

for i = 1:size(EEG_epoched.EEG_P6_flexion.Time_Domain.EEG_times, 3)
    if i ~= max_L_P6_flexion_indx
        x_old = EEG_epoched.EEG_P6_flexion.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P6_flexion);

        % channels data
        data = EEG_epoched.EEG_P6_flexion.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_flexion.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_flexion.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P6_flexion.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_flexion.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_flexion.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end



%% EEG_P1_extension_lenght_normalization
[max_L_P1_extension, max_L_P1_extension_indx] = max(L_P1_extension);

for i = 1:size(EEG_epoched.EEG_P1_extension.Time_Domain.EEG_times, 3)
    if i ~= max_L_P1_extension_indx
        x_old = EEG_epoched.EEG_P1_extension.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P1_extension);

        % channels data
        data = EEG_epoched.EEG_P1_extension.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_extension.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_extension.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P1_extension.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_extension.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_extension.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P3_extension_lenght_normalization
[max_L_P3_extension, max_L_P3_extension_indx] = max(L_P3_extension);

for i = 1:size(EEG_epoched.EEG_P3_extension.Time_Domain.EEG_times, 3)
    if i ~= max_L_P3_extension_indx
        x_old = EEG_epoched.EEG_P3_extension.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P3_extension);

        % channels data
        data = EEG_epoched.EEG_P3_extension.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_extension.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_extension.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P3_extension.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_extension.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_extension.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P6_extension_lenght_normalization
[max_L_P6_extension, max_L_P6_extension_indx] = max(L_P6_extension);

for i = 1:size(EEG_epoched.EEG_P6_extension.Time_Domain.EEG_times, 3)
    if i ~= max_L_P6_extension_indx
        x_old = EEG_epoched.EEG_P6_extension.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P6_extension);

        % channels data
        data = EEG_epoched.EEG_P6_extension.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_extension.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_extension.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P6_extension.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_extension.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_extension.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P1_flextoflex_lenght_normalization
[max_L_P1_flextoflex, max_L_P1_flextoflex_indx] = max(L_P1_flextoflex);

for i = 1:size(EEG_epoched.EEG_P1_flextoflex.Time_Domain.EEG_times, 3)
    if i ~= max_L_P1_flextoflex_indx
        x_old = EEG_epoched.EEG_P1_flextoflex.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P1_flextoflex);

        % channels data
        data = EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P3_flextoflex_lenght_normalization
[max_L_P3_flextoflex, max_L_P3_flextoflex_indx] = max(L_P3_flextoflex);

for i = 1:size(EEG_epoched.EEG_P3_flextoflex.Time_Domain.EEG_times, 3)
    if i ~= max_L_P3_flextoflex_indx
        x_old = EEG_epoched.EEG_P3_flextoflex.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P3_flextoflex);

        % channels data
        data = EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end




%% EEG_P6_flextoflex_lenght_normalization
[max_L_P6_flextoflex, max_L_P6_flextoflex_indx] = max(L_P6_flextoflex);

for i = 1:size(EEG_epoched.EEG_P6_flextoflex.Time_Domain.EEG_times, 3)
    if i ~= max_L_P6_flextoflex_indx
        x_old = EEG_epoched.EEG_P6_flextoflex.Time_Domain.EEG_times{1, 1, i};
        x_new = linspace(x_old(1), x_old(end), max_L_P6_flextoflex);

        % channels data
        data = EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);
        

        % source data
        data = EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
        S = size(EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.not_normalized{1, 1, 1}, 1);
        interpolated_data = arrayfun(@(row) interp1(x_old, data(row, :), ...
            x_new, 'spline'), (1:S)', 'UniformOutput', false);

        EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.length_normalized{1, 1, i} = ...
            cell2mat(interpolated_data);

    end
end







%% Calculate PSD (frequency domain)

Fs = 500; % Sampling frequency
nfft = 2048; % number of FFT points


%% EEG_P1_trials
for i = 1:size(EEG_epoched.EEG_P1_trials.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P1_trials.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_trials.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P1_trials.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_trials.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P1_trials.Frequency_Domain.freqs = F;




%% EEG_P3_trials
for i = 1:size(EEG_epoched.EEG_P3_trials.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P3_trials.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_trials.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P3_trials.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_trials.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P3_trials.Frequency_Domain.freqs = F;




%% EEG_P6_trials
for i = 1:size(EEG_epoched.EEG_P6_trials.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P6_trials.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_trials.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P6_trials.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/10); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_trials.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P6_trials.Frequency_Domain.freqs = F;




%% EEG_P1_flexion
for i = 1:size(EEG_epoched.EEG_P1_flexion.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P1_flexion.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_flexion.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P1_flexion.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_flexion.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P1_flexion.Frequency_Domain.freqs = F;




%% EEG_P3_flexion
for i = 1:size(EEG_epoched.EEG_P3_flexion.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P3_flexion.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_flexion.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P3_flexion.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_flexion.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P3_flexion.Frequency_Domain.freqs = F;




%% EEG_P6_flexion
for i = 1:size(EEG_epoched.EEG_P6_flexion.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P6_flexion.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_flexion.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P6_flexion.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_flexion.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P6_flexion.Frequency_Domain.freqs = F;




%% EEG_P1_extension
for i = 1:size(EEG_epoched.EEG_P1_extension.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P1_extension.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_extension.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P1_extension.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_extension.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P1_extension.Frequency_Domain.freqs = F;




%% EEG_P3_extension
for i = 1:size(EEG_epoched.EEG_P3_extension.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P3_extension.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_extension.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P3_extension.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_extension.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P3_extension.Frequency_Domain.freqs = F;




%% EEG_P6_extension
for i = 1:size(EEG_epoched.EEG_P6_extension.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P6_extension.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_extension.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P6_extension.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_extension.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P6_extension.Frequency_Domain.freqs = F;




%% EEG_P1_flextoflex
for i = 1:size(EEG_epoched.EEG_P1_flextoflex.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P1_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_flextoflex.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P1_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P1_flextoflex.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P1_flextoflex.Frequency_Domain.freqs = F;




%% EEG_P3_flextoflex
for i = 1:size(EEG_epoched.EEG_P3_flextoflex.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P3_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_flextoflex.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P3_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P3_flextoflex.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P3_flextoflex.Frequency_Domain.freqs = F;




%% EEG_P6_flextoflex
for i = 1:size(EEG_epoched.EEG_P6_flextoflex.Time_Domain.EEG_times, 3)
    
    % channels
    signal = EEG_epoched.EEG_P6_flextoflex.Time_Domain.channels.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, ~] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_flextoflex.Frequency_Domain.channels{1, 1, i} = Pxx';

    % sources
    signal = EEG_epoched.EEG_P6_flextoflex.Time_Domain.sources.not_normalized{1, 1, i};
    % Define parameters for the PSD calculation
    window = floor(size(signal,2)/2); % length of each segment
    noverlap = floor(0.75*window); % number of samples to overlap between segments

    [Pxx, F] = pwelch(signal', window, noverlap, nfft, Fs);
    EEG_epoched.EEG_P6_flextoflex.Frequency_Domain.sources{1, 1, i} = Pxx';

end
F = F';
EEG_epoched.EEG_P6_flextoflex.Frequency_Domain.freqs = F;



%% save data separately

filepath = ['C:\Morteza\Analysis\ANSYMB2024\data\', ...
    '5_single-subject-EEG-analysis_with_epochs\', ...
    'sub-', num2str(subject_id)];


%% trials
EEG_epoched_trials = struct('EEG_P1_trials', EEG_epoched.EEG_P1_trials, ...
    'EEG_P3_trials', EEG_epoched.EEG_P3_trials, ...
    'EEG_P6_trials', EEG_epoched.EEG_P6_trials);

filename = 'EEG_epoched_trials.mat';
save(fullfile(filepath, filename), 'EEG_epoched_trials', '-v7.3');

%% flexion
EEG_epoched_flexion = struct('EEG_P1_flexion', EEG_epoched.EEG_P1_flexion, ...
    'EEG_P3_flexion', EEG_epoched.EEG_P3_flexion, ...
    'EEG_P6_flexion', EEG_epoched.EEG_P6_flexion);

filename = 'EEG_epoched_flexion.mat';
save(fullfile(filepath, filename), 'EEG_epoched_flexion', '-v7.3');

%% extension
EEG_epoched_extension = struct('EEG_P1_extension', EEG_epoched.EEG_P1_extension, ...
    'EEG_P3_extension', EEG_epoched.EEG_P3_extension, ...
    'EEG_P6_extension', EEG_epoched.EEG_P6_extension);

filename = 'EEG_epoched_extension.mat';
save(fullfile(filepath, filename), 'EEG_epoched_extension', '-v7.3');

%% flextoflex
EEG_epoched_flextoflex = struct('EEG_P1_flextoflex', EEG_epoched.EEG_P1_flextoflex, ...
    'EEG_P3_flextoflex', EEG_epoched.EEG_P3_flextoflex, ...
    'EEG_P6_flextoflex', EEG_epoched.EEG_P6_flextoflex);

filename = 'EEG_epoched_flextoflex.mat';
save(fullfile(filepath, filename), 'EEG_epoched_flextoflex', '-v7.3');


%% clear the big data to save memory
clear EEG_epoched

%% Plot the results
load('Channels_Names.mat');
Ch_Names = Channels_Names;
All_colours = struct('dark_blue', [0, 0.4470, 0.7410], 'light_blue', [0.3010, 0.7450, 0.9330], ...
    'dark_orange', [0.8500, 0.3250, 0.0980], 'light_orange', [0.9290, 0.6940, 0.1250], ...
    'dark_green', [0.4660, 0.6740, 0.1880], 'light_green', [0.5960, 0.8740, 0.5410]);


%% EEG_trials - Time-Domain: channels
%%% making 3d matrix 
numChannels = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.channels.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.channels.length_normalized , 3);
signal_P1_3d = zeros(numChannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P1_trials.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_trials.EEG_P3_trials.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_trials.EEG_P3_trials.Time_Domain.channels.length_normalized , 3);
signal_P3_3d = zeros(numChannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P3_trials.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_trials.EEG_P6_trials.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_trials.EEG_P6_trials.Time_Domain.channels.length_normalized , 3);
signal_P6_3d = zeros(numChannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P6_trials.Time_Domain.channels.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Entire Trial): %d trials P1, %d trials P3, %d trials P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_trials - Time-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.sources.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_trials.EEG_P1_trials.Time_Domain.sources.length_normalized , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P1_trials.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_trials.EEG_P3_trials.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_trials.EEG_P3_trials.Time_Domain.sources.length_normalized , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P3_trials.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_trials.EEG_P6_trials.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_trials.EEG_P6_trials.Time_Domain.sources.length_normalized , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_trials.EEG_P6_trials.Time_Domain.sources.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Entire Trial): %d trials P1, %d trials P3, %d trials P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_trials - Frequency-Domain: channels
%%% making 3d matrix 
numchannels = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.channels{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.channels , 3);
signal_P1_3d = zeros(numchannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.channels{1,1, i});
end

numPointsP3   = size(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.channels , 3);
signal_P3_3d = zeros(numchannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.channels{1,1, i});
end

numPointsP6   = size(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.channels , 3);
signal_P6_3d = zeros(numchannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.channels{1,1, i});
end


%
XP1 = EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.freqs;
XP3 = EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.freqs;
XP6 = EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Entire Trial): %d trials P1, %d trials P3, %d trials P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['Frequency-Domain: ', Ch_Names{i}], All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_trials - Frequency-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.sources{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.sources , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.sources{1,1, i});
end

numPointsP3   = size(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.sources , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.sources{1,1, i});
end

numPointsP6   = size(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.sources , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.sources{1,1, i});
end


%
XP1 = EEG_epoched_trials.EEG_P1_trials.Frequency_Domain.freqs;
XP3 = EEG_epoched_trials.EEG_P3_trials.Frequency_Domain.freqs;
XP6 = EEG_epoched_trials.EEG_P6_trials.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Entire Trial): %d trials P1, %d trials P3, %d trials P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['Frequency-Domain: IC ', num2str(i)], All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flexion - Time-Domain: channels
%%% making 3d matrix 
numChannels = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.channels.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.channels.length_normalized , 3);
signal_P1_3d = zeros(numChannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.channels.length_normalized , 3);
signal_P3_3d = zeros(numChannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.channels.length_normalized , 3);
signal_P6_3d = zeros(numChannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.channels.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Flexion): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flexion - Time-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.sources.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.sources.length_normalized , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P1_flexion.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.sources.length_normalized , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P3_flexion.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.sources.length_normalized , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_flexion.EEG_P6_flexion.Time_Domain.sources.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Entire Trial): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flexion - Frequency-Domain: channels
%%% making 3d matrix 
numchannels = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.channels{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.channels , 3);
signal_P1_3d = zeros(numchannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.channels{1,1, i});
end

numPointsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.channels , 3);
signal_P3_3d = zeros(numchannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.channels{1,1, i});
end

numPointsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.channels , 3);
signal_P6_3d = zeros(numchannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.channels{1,1, i});
end


%
XP1 = EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.freqs;
XP3 = EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.freqs;
XP6 = EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Entire Trial): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flexion - Frequency-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.sources{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.sources , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.sources{1,1, i});
end

numPointsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.sources , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.sources{1,1, i});
end

numPointsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.sources , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.sources{1,1, i});
end


%
XP1 = EEG_epoched_flexion.EEG_P1_flexion.Frequency_Domain.freqs;
XP3 = EEG_epoched_flexion.EEG_P3_flexion.Frequency_Domain.freqs;
XP6 = EEG_epoched_flexion.EEG_P6_flexion.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Entire Trial): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_extension - Time-Domain: channels
%%% making 3d matrix 
numChannels = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.channels.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.channels.length_normalized , 3);
signal_P1_3d = zeros(numChannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P1_extension.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_extension.EEG_P3_extension.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_extension.EEG_P3_extension.Time_Domain.channels.length_normalized , 3);
signal_P3_3d = zeros(numChannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P3_extension.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_extension.EEG_P6_extension.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_extension.EEG_P6_extension.Time_Domain.channels.length_normalized , 3);
signal_P6_3d = zeros(numChannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P6_extension.Time_Domain.channels.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Extension): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_extension - Time-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.sources.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_extension.EEG_P1_extension.Time_Domain.sources.length_normalized , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P1_extension.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_extension.EEG_P3_extension.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_extension.EEG_P3_extension.Time_Domain.sources.length_normalized , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P3_extension.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_extension.EEG_P6_extension.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_extension.EEG_P6_extension.Time_Domain.sources.length_normalized , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_extension.EEG_P6_extension.Time_Domain.sources.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (Extension): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_extension - Frequency-Domain: channels
%%% making 3d matrix 
numchannels = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.channels{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.channels , 3);
signal_P1_3d = zeros(numchannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.channels{1,1, i});
end

numPointsP3   = size(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.channels , 3);
signal_P3_3d = zeros(numchannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.channels{1,1, i});
end

numPointsP6   = size(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.channels , 3);
signal_P6_3d = zeros(numchannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.channels{1,1, i});
end


%
XP1 = EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.freqs;
XP3 = EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.freqs;
XP6 = EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Extension): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_extension - Frequency-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.sources{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.sources , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.sources{1,1, i});
end

numPointsP3   = size(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.sources , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.sources{1,1, i});
end

numPointsP6   = size(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.sources , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.sources{1,1, i});
end


%
XP1 = EEG_epoched_extension.EEG_P1_extension.Frequency_Domain.freqs;
XP3 = EEG_epoched_extension.EEG_P3_extension.Frequency_Domain.freqs;
XP6 = EEG_epoched_extension.EEG_P6_extension.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (Extension): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom





%% EEG_flextoflex - Time-Domain: channels
%%% making 3d matrix 
numChannels = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.channels.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.channels.length_normalized , 3);
signal_P1_3d = zeros(numChannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.channels.length_normalized , 3);
signal_P3_3d = zeros(numChannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.channels.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.channels.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.channels.length_normalized , 3);
signal_P6_3d = zeros(numChannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.channels.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (flextoflex): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flextoflex - Time-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.sources.length_normalized{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.sources.length_normalized , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P1_flextoflex.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.sources.length_normalized , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P3_flextoflex.Time_Domain.sources.length_normalized{1,1, i};
end

numPointsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.sources.length_normalized{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.sources.length_normalized , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        EEG_epoched_flextoflex.EEG_P6_flextoflex.Time_Domain.sources.length_normalized{1,1, i};
end


%
XP1 = linspace(0, 100, numPointsP1);
XP3 = linspace(0, 100, numPointsP3);
XP6 = linspace(0, 100, numPointsP6);

mean_signalP1 = mean(signal_P1_3d, 3);
mean_signalP3 = mean(signal_P3_3d, 3);
mean_signalP6 = mean(signal_P6_3d, 3);

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Cycle [%]';
YLabel_n = 'Amplitude [\muV]';
XLim_n = [0 100];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Time-Domain (flextoflex): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flextoflex - Frequency-Domain: channels
%%% making 3d matrix 
numchannels = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.channels{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.channels , 3);
signal_P1_3d = zeros(numchannels, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.channels{1,1, i});
end

numPointsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.channels , 3);
signal_P3_3d = zeros(numchannels, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.channels{1,1, i});
end

numPointsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.channels{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.channels , 3);
signal_P6_3d = zeros(numchannels, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.channels{1,1, i});
end


%
XP1 = EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.freqs;
XP3 = EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.freqs;
XP6 = EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (flextoflex): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:64
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(Ch_Names{i});
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, Ch_Names{i}, XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom




%% EEG_flextoflex - Frequency-Domain: sources
%%% making 3d matrix 
numsources = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.sources{1,1,1} , 1);

numPointsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP1   = size(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.sources , 3);
signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
for i = numEpochsP1
    signal_P1_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.sources{1,1, i});
end

numPointsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP3   = size(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.sources , 3);
signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
for i = numEpochsP3
    signal_P3_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.sources{1,1, i});
end

numPointsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.sources{1,1,1} , 2);
numEpochsP6   = size(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.sources , 3);
signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
for i = numEpochsP6
    signal_P6_3d(:,:, i) = ...
        10*log10(EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.sources{1,1, i});
end


%
XP1 = EEG_epoched_flextoflex.EEG_P1_flextoflex.Frequency_Domain.freqs;
XP3 = EEG_epoched_flextoflex.EEG_P3_flextoflex.Frequency_Domain.freqs;
XP6 = EEG_epoched_flextoflex.EEG_P6_flextoflex.Frequency_Domain.freqs;

mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');

std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

% plot the results
h1 = []; h2 = []; h3 = [];
XLabel_n = 'Frequency [Hz]';
YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
XLim_n = [0.5 50];
figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
tiledlayout(8,8)
sgtitle(sprintf('Frequency-Domain (flextoflex): %d epochs P1, %d epochs P3, %d epochs P6', numEpochsP1, numEpochsP3, numEpochsP6));
for i = 1:61
    nexttile; hold on

    % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
    %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
    %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);

    % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
    %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
    %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);

    % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
    %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
    %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
    
    set(gca, 'XLim', [0.5 50])

    title(['IC ', num2str(i)]);
    ax = gca;
    set(ax, 'ButtonDownFcn', ...
        @(src, event)showDetails({XP1, XP3, XP6}, ...
        {mean_signalP1, mean_signalP3, mean_signalP6}, ...
        {std_signalP1, std_signalP3, std_signalP6}, ...
        i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
end

% Create a legend for the entire tiled layout
lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south'; % Position the legend at the bottom






%% Callback function to display detailed plot
function showDetails(X, mean_signals, std_signals, ch, title_n, XLabel_n, YLabel_n, XLim_n, All_colours)
    figure; hold on;
    alpha = 0.2;

    fill([X{1,1} fliplr(X{1,1})], [mean_signals{1,1}(ch, :) + std_signals{1,1}(ch, :), ...
        fliplr(mean_signals{1,1}(ch, :) - std_signals{1,1}(ch, :))], ...
        All_colours.light_blue, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    fill([X{1,2} fliplr(X{1,2})], [mean_signals{1,2}(ch, :) + std_signals{1,2}(ch, :), ...
        fliplr(mean_signals{1,2}(ch, :) - std_signals{1,2}(ch, :))], ...
        All_colours.light_orange, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    fill([X{1,3} fliplr(X{1,3})], [mean_signals{1,3}(ch, :) + std_signals{1,3}(ch, :), ...
        fliplr(mean_signals{1,3}(ch, :) - std_signals{1,3}(ch, :))], ...
        All_colours.light_green, 'FaceAlpha', alpha, 'EdgeColor', 'none');


    h1 = plot(X{1,1}, mean_signals{1,1}(ch, :), 'Color', All_colours.dark_blue, 'LineWidth', 2); 

    h2 = plot(X{1,2}, mean_signals{1,2}(ch, :), 'Color', All_colours.dark_orange, 'LineWidth', 2); 
    
    h3 = plot(X{1,3}, mean_signals{1,3}(ch, :), 'Color', All_colours.dark_green, 'LineWidth', 2); 
    
    set(gca, 'XLim', XLim_n)

    title(title_n)
    xlabel(XLabel_n);
    ylabel(YLabel_n);
    grid on;

    legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal')

    zoom on
end




