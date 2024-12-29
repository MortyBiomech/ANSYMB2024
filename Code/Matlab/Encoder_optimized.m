clc; clear; 

%% Add LSL library
addpath(genpath('C:\Morteza\LSL\liblsl-Matlab'))

%% Setting the variables
Num_trials = 80; % How many times should the (1, 3, 6 bar) conditions repeat; 
Num_Conditions = 3; % Number of conditions (1, 3, 6 bar)

up_lim = -5; % Upper limit of the movement
lo_lim = -85; % Lower limit of the movement
freq_avg = 0.5; % Frequency of the reference movement

%% defining the Arduino
disp('Defining the Arduino...');
a = arduino('COM8','Mega2560','Libraries','rotaryEncoder')
disp('Defining the Encoder...');
encoder = rotaryEncoder(a,'D2','D3',1000)

%% Send data to LSL
% instantiate the library
disp('Loading library...');
lib = lsl_loadlib();
% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Encoder_Pressure_Preference_Force', ...
    'Exp_data',7,0,'cf_float32', 'myuid123456');
disp('Opening an outlet...');
outlet = lsl_outlet(info);
% send data into the outlet, sample by sample
disp('Now transmitting data...');

%% Setting the initial position
resetCount(encoder); Pressure = 0; Preference = 0; mappedPressure = 0;
auditory_input = 0;
score_press = 0;
phaseValue = 0;
sessionNumber = 0;

%% Read the data from the encoder
% Parameters
windowSize = 500; % Number of samples in the sliding window
currentIdx = 0;   % Initialize circular buffer index

% Preallocate buffers
Angle = nan(1, windowSize); % Circular buffer for angle
Time_enc = nan(1, windowSize); % Circular buffer for time
Ref_trj = nan(1, windowSize); % Circular buffer for reference trajectory
Up_lim = up_lim * ones(1, windowSize); % Fixed upper limit
Lo_lim = lo_lim * ones(1, windowSize); % Fixed lower limit

figure();
xlabel('Time (s)')
ylabel('Knee Angle (deg)')

hhh = zeros(1, 2000);
i = 1;

while true
    tic
    % Read encoder data
    [count, time] = readCount(encoder);
    angle = (360 * count) / 4000;

    % Map pressure to bar values
    switch Pressure
        case 1, mappedPressure = 1; % 1 bar
        case 2, mappedPressure = 3; % 3 bar
        case 3, mappedPressure = 6; % 6 bar
    end
    writePWMDutyCycle(a, 'D13', mappedPressure / 6);

    % Compute reference trajectory
    ref_trj = abs((lo_lim - up_lim) / 2) * sin(2 * pi * freq_avg * time) + ...
        (((lo_lim - up_lim) / 2) + up_lim);


    % Update circular buffer
    currentIdx = mod(currentIdx, windowSize) + 1; % Increment circular index
    Angle(currentIdx) = angle;
    Time_enc(currentIdx) = time;
    Ref_trj(currentIdx) = ref_trj;

    % Normalize time for plotting
    normalizedTime = Time_enc - Time_enc(currentIdx); % Subtract current time for alignment

    % Plot sliding window
    plot(normalizedTime, Ref_trj, ...
         normalizedTime, Up_lim, 'r', ...
         normalizedTime, Lo_lim, 'r', ...
         normalizedTime, Angle, 'b', ...
         0, angle, 'bo', 'MarkerSize', 6, 'LineWidth', 2);

    grid on;
    ylim([lo_lim - 20, up_lim + 20]);
    % xlim([-max(normalizedTime), 0]); % Display only the latest `windowSize` samples



    % Force sensor
    force = readVoltage(a, 'A2');
    % Force = [Force; force];



    Encoder_Pressure_Preference_Force = ...
        [angle, ref_trj, mappedPressure, Preference, force, auditory_input, score_press];
    outlet.push_sample(Encoder_Pressure_Preference_Force);

    hhh(i) = toc;
    i = i+1;
end
