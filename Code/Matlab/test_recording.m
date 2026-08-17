%% Load dataset

filepath = 'C:\Users\morte\OneDrive\Documents\CurrentStudy\sub-AMAMtest\ses-S001\eeg';
filename = 'sub-AMAMtest_ses-S001_task-Default_run-001_eeg.xdf';

data = load_xdf(fullfile(filepath, filename));


%% Apply a band-pass filter to demonstrate the EEG signals
eeg_id = 1;
fs = 500;           % Sampling frequency in Hz
EEG = data{1, eeg_id}.time_series;
EEG = EEG';

% 4th-order Butterworth bandpass filter from 1 to 30 Hz
[b, a] = butter(4, [1 30]/(fs/2), 'bandpass');

% Apply the filter
EEG_filt = filtfilt(b, a, EEG); % zero-phase filtering (recommended)
EEG_filt = EEG_filt';


