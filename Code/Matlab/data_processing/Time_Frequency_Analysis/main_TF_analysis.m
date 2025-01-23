clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];


%% load one data for generating sample code
subject = 18;
cd([epoched_data_path, 'sub-', num2str(subject)])

data = load("Epochs_FlextoFlex_based.mat");
name = fieldnames(data);
data = data.(name{1});


%% Extract the TF content using Morlet Wavelet 
fs = 500; % sampling frequency
trial = 55;
epoch = 1;
ic = 1;
signal = data{1, trial}.EEG_stream.Preprocessed.Sources{1, epoch};
time = data{1, trial}.EEG_stream.Preprocessed.Times{1, epoch};

%% Compute the Continuous Wavelet Transform using Morlet Wavelets
[cwt_coeffs, frequencies] = ...
    cwt(signal, 'amor', fs, 'FrequencyLimits', [1 50], 'VoicesPerOctave', 20); 

% Compute power (magnitude squared of coefficients)
power = abs(cwt_coeffs).^2;

% Calculate mean power across time for each frequency
mean_power_per_frequency = mean(power, 2); % Average across time (rows)

% Normalize power using decibel scale
normalized_power_db = 10 * log10(bsxfun(@rdivide, power, mean_power_per_frequency));


%%
figure;
surf(time, frequencies, normalized_power_db, 'EdgeColor', 'none');

axis tight;
view(2);

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Normalized Time-Frequency Plot (dB)');
colorbar;
colormap jet;
clim([-20 20]); % Adjust color limits if needed


%%
% Parameters
Fs = 500; % Sampling rate in Hz
t = -2:1/Fs:2; % Time vector for Morlet wavelet
f_min = 1; % Minimum frequency in Hz
f_max = 50; % Maximum frequency in Hz
num_frequencies = 250; % Number of frequency bands
frequencies = logspace(log10(f_min), log10(f_max), num_frequencies); % Log-spaced frequencies
width = linspace(4, 10, num_frequencies); % Width of the Morlet wavelet (number of cycles)

% Example EEG signal
EEG_signal = signal; 

% Preallocate output
power = zeros(num_frequencies, length(EEG_signal));

% Morlet wavelet transform
for fi = 1:num_frequencies
    freq = frequencies(fi);
    w = width(fi);
    
    % Create Morlet wavelet
    sigma_t = w / (2 * pi * freq); % Temporal standard deviation
    wavelet = exp(2 * 1i * pi * freq * t) .* exp(-t.^2 / (2 * sigma_t^2)); % Complex Morlet wavelet
    wavelet = wavelet / sqrt(sum(abs(wavelet).^2)); % Normalize wavelet
    
    % Convolve wavelet with EEG signal
    convolution = conv(EEG_signal, wavelet, 'same');
    
    % Compute power
    power(fi, :) = abs(convolution).^2;
end

% Calculate mean power across time for each frequency
mean_power_per_frequency = mean(power, 2); % Average across time (rows)

% Normalize power using decibel scale
normalized_power_db = 10 * log10(bsxfun(@rdivide, power, mean_power_per_frequency));


% Subtract baseline from power
ERSP = power - mean_power_per_frequency;

%% Time-frequency plot
figure;
imagesc(time, frequencies, ERSP);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Time-Frequency Representation');


% % Define extreme colors
% synch_color = [214, 40, 40]/255; 
% desynch_color = [58, 134, 255]/255;  
% 
% % Define the number of colors for the colormap
% num_colors = 256;
% 
% % Create a gradient for the negative side (red to white)
% neg_colors = [linspace(desynch_color(1), 1, num_colors/2)', ...
%               linspace(desynch_color(2), 1, num_colors/2)', ...
%               linspace(desynch_color(3), 1, num_colors/2)'];
% 
% % Create a gradient for the positive side (white to blue)
% pos_colors = [linspace(1, synch_color(1), num_colors/2)', ...
%               linspace(1, synch_color(2), num_colors/2)', ...
%               linspace(1, synch_color(3), num_colors/2)'];
% 
% % Combine the gradients to create the full colormap
% custom_cmap = [neg_colors; pos_colors];

% Define extreme colors
synch_color = [214, 40, 40]/255; % Red (positive side)
desynch_color = [58, 134, 255]/255; % Blue (negative side)

% Define the number of colors for the colormap
num_colors = 256;

% Create a gradient for the negative side (blue to white) using logspace
neg_colors = [logspace(log10(desynch_color(1)), log10(1), num_colors/2)', ...
              logspace(log10(desynch_color(2)), log10(1), num_colors/2)', ...
              logspace(log10(desynch_color(3)), log10(1), num_colors/2)'];

% Create a gradient for the positive side (white to red) using logspace
pos_colors = [logspace(log10(1), log10(synch_color(1)), num_colors/2)', ...
              logspace(log10(1), log10(synch_color(2)), num_colors/2)', ...
              logspace(log10(1), log10(synch_color(3)), num_colors/2)'];

% Combine the gradients to create the full colormap
custom_cmap = [neg_colors; pos_colors];

colormap(custom_cmap)
% Adjust the color axis to center around zero
clim([-max(abs(ERSP(:))), max(abs(ERSP(:)))]);

colorbar;

% figure;
% imagesc(time, frequencies, ERSP);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Time-Frequency Representation (Morlet Wavelet) - subtraction');
% colorbar;
% 
% figure;
% imagesc(time, frequencies, normalized_power_db);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Time-Frequency Representation (Morlet Wavelet) - dB');
% colorbar;
% clim([-1 1])
