clc
clear

%% load XDF file
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master'));
addpath(genpath('C:\Morteza\Analysis\ANSYMB2024'));

subject = 1;
files_path = ['C:\Morteza\Analysis\ANSYMB2024\Subjects_Dataset\P00', ...
    num2str(subject), '\sub-P00', num2str(subject), '\ses-S001\eeg'];
filename1 = ['\sub-P00', num2str(subject), '_ses-S001_task-Default_run-001_eeg.xdf'];
filename2 = ['\sub-P00', num2str(subject), '_ses-S001_task-Default_run-002_eeg.xdf'];
filename3 = ['\sub-P00', num2str(subject), '_ses-S001_task-Default_run-003_eeg.xdf'];
filename4 = ['\sub-P00', num2str(subject), '_ses-S001_task-Default_run-004_eeg.xdf'];

%% diagnosis 
% files_path = ['C:\Morteza\Analysis\ANSYMB2024\Subjects_Dataset\P00000', ...
%      '\sub-P00000diag\ses-S001\eeg\'];
% filename1 = ['sub-P00000diag_ses-S001_task-Default_run-011_eeg.xdf'];

files_path = 'C:\Users\morte\OneDrive\Documents\CurrentStudy\sub-P005\ses-S001\eeg\';
filename1 = 'sub-P005_ses-S001_task-Default_run-002_eeg.xdf';

%%
streams1 = load_xdf([files_path, filename1]);
streams2 = load_xdf([files_path, filename2]);
streams3 = load_xdf([files_path, filename3]);
streams4 = load_xdf([files_path, filename4]);

%%
figure()
plot(streams1{1, 2}.time_stamps - streams1{1, 2}.time_stamps(1), ...
    streams1{1, 2}.time_series(6, :))
hold on
plot(streams{1, 1}.time_stamps - streams{1, 1}.time_stamps(1), ...
    streams{1, 1}.time_series(12, :))
set(gca, 'YLim', [-5 5])
xlabel('Time [s]')
ylabel('mV')


%% Check EEG signals
h = figure();
tiledlayout(3,1)
nexttile
plot(streams{1,2}.time_stamps - streams{1,2}.time_stamps(1), ...
    streams{1,2}.time_series(11,:))
nexttile
plot(streams{1,2}.time_stamps - streams{1,2}.time_stamps(1), ...
    streams{1,2}.time_series(12,:))
nexttile
plot(streams{1,2}.time_stamps - streams{1,2}.time_stamps(1), ...
    streams{1,2}.time_series(13,:))



%% check EMG signals
figure()
tiledlayout(4,1)
nexttile
plot(streams{1,4}.time_stamps - streams{1,4}.time_stamps(1), ...
    streams{1,4}.time_series(1,:))
nexttile
plot(streams{1,4}.time_stamps - streams{1,4}.time_stamps(1), ...
    streams{1,4}.time_series(2,:))
nexttile
plot(streams{1,4}.time_stamps - streams{1,4}.time_stamps(1), ...
    streams{1,4}.time_series(3,:))
nexttile
plot(streams{1,4}.time_stamps - streams{1,4}.time_stamps(1), ...
    streams{1,4}.time_series(4,:))





%% test plots
figure()
subplot(2,1,1)
plot(streams{1, 3}.time_stamps, streams{1, 3}.time_series(end,:))
set(gca, 'XLim', [1.078 1.079]*10000)
subplot(2,1,2)
plot(streams{1, 3}.time_stamps, streams{1, 3}.time_series(1,:))
set(gca, 'XLim', [1.078 1.079]*10000)
