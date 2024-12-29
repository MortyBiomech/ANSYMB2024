
%% Load dataset

filepath = 'C:\Morteza\MyProjects\ANSYMB2024\data\0_source_data\sub-18\ses-S004\eeg';
filename = 'sub-18_ses-S004_task-Default_run-001_eeg.xdf';

% filepath = 'C:\Users\morte\OneDrive\Documents\CurrentStudy\sub-100\ses-S001\eeg';
% filename = 'sub-100_ses-S001_task-Default_run-001_eeg.xdf';

data = load_xdf(fullfile(filepath, filename));


%% Plot experiment data stream
figure();
i1 = 3; % cell id
X = 5; % force sensor
plot(data{1,i1}.time_stamps - data{1,i1}.time_stamps(1), data{1,i1}.time_series(X,:))
title('Force Sensor')


%% Muscles EMG
figure();
i2 = 4; % cell id 

tiledlayout(3, 2);
nexttile
X = 2; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])

nexttile
X = 3; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])

nexttile
X = 4; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])

nexttile
X = 5; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])

nexttile
X = 8; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])

nexttile
X = 9; % EMG id
plot(data{1,i2}.time_stamps - data{1,i2}.time_stamps(1), data{1,i2}.time_series(X,:))
title(['EMG sensor ', num2str(X)])


