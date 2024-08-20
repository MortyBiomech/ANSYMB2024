clc; clear; close all;

%% Add LSL library
addpath(genpath('C:\Morteza\LSL\liblsl-Matlab'))

%% Setting the variables
Num_trials = 62; % How many times should the (1, 3, 6 bar) conditions repeat; 
Num_Conditions = 3; % Number of conditions (1, 3, 6 bar)

up_lim = 15; % Upper limit of the movement
lo_lim = -60; % Lower limit of the movement
freq_avg = 0.6; % Frequency of the reference movement

%% defining the Arduino
disp('Defining the Arduino...');
a = arduino('COM3','Mega2560','Libraries','rotaryEncoder')
disp('Defining the Encoder...');
encoder = rotaryEncoder(a,'D2','D3',1000)

%% Send data to LSL
% instantiate the library
disp('Loading library...');
lib = lsl_loadlib();
% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Encoder_Pressure_Preference_Force', ...
    'Exp_data',6,0,'cf_float32', 'myuid123456');
disp('Opening an outlet...');
outlet = lsl_outlet(info);
% send data into the outlet, sample by sample
disp('Now transmitting data...');

%% Setting the initial position
resetCount(encoder); Pressure = 0; Preference = 0; mappedPressure = 0;
auditory_input = 0;

%% Read the data from the encoder
Angle = []; Time_enc = []; time = 0; StartSpot = 0; Up_lim = []; Lo_lim = []; 
Ref_trj = []; ref_trj = 0; Force = [];

figure1 = figure();
% figure2 = figure();

while true
    
    [count,time] = readCount(encoder);
    angle = (360*count)/4000;

    % Map the incoming Pressure to the desired pressure values
    switch Pressure
        case 1
            mappedPressure = 1; % Set to 1 bar
        case 2
            mappedPressure = 3; % Set to 3.5 bar
        case 3
            mappedPressure = 6; % Set to 6 bar
    end

    % Write the PWM duty cycle based on the mapped pressure
    writePWMDutyCycle(a, 'D13', mappedPressure / 6);
    % writePWMDutyCycle(a, 'D13', Pressure / 6);

    

    % Limited Movement
    % figure(figure1) % plot tracking path for the participant
    ylabel('Knee Angle (deg)')
    xlabel('Time (s)')
    ref_trj = abs((lo_lim-up_lim)/2) * sin(2*pi*freq_avg*time) + ...
        (((lo_lim-up_lim)/2) + up_lim);
    
    Angle = [Angle; angle]; Time_enc = [Time_enc; time];
    Up_lim = [Up_lim; up_lim]; Lo_lim = [Lo_lim; lo_lim]; Ref_trj = [Ref_trj; ref_trj];
    
    p = plot(Time_enc+5, Ref_trj, ...
        Time_enc+5, Up_lim,'r', Time_enc+5, Lo_lim,'r', ...
        Time_enc, Angle,'b', time, angle,'b');
    
    if (time-20 < 0)
        StartSpot = 0;
    else
        StartSpot = time-20;
    end
    
    p(1).LineWidth = 5; 
    p(2).LineWidth = 3; p(3).LineWidth = 3; 
    p(4).LineWidth = 3; 
    p(5).LineWidth = 3; p(5).Marker = 'o'; p(5).MarkerSize = 6; 
    p(1).Color = [0.8 0.8 0.8];
    
    grid on

    axis([StartSpot, (time+5), lo_lim-20, up_lim+20])


    % figure(figure2) % plot tracking path for control 
    % subplot(2,1,1)
    % p1 = plot(Time_enc+5, Ref_trj, ...
    %     Time_enc+5, Up_lim,'r', Time_enc+5, Lo_lim,'r', ...
    %     Time_enc, Angle,'b', time, angle,'b');
    % p1(1).LineWidth = 5; 
    % p1(2).LineWidth = 3; p1(3).LineWidth = 3; 
    % p1(4).LineWidth = 3; 
    % p1(5).LineWidth = 3; p1(5).Marker = 'o'; p1(5).MarkerSize = 6; 
    % p1(1).Color = [0.8 0.8 0.8];
    % grid on
    % axis([StartSpot, (time+5), lo_lim-20, up_lim+20])

    % Force sensor
    force = readVoltage(a, 'A1');
    Force = [Force; force];

    % subplot(2,1,2)
    % figure(figure2);
    % p2 = plot(time, force, 'b', Time_enc, Force, 'b');
    % 
    % p2(1).LineWidth = 3; p2(2).LineWidth = 3; 
    % p2(1).Marker = 'o'; p2(1).MarkerSize = 6;
    % grid on
    % axis([StartSpot, (time+5), 0, 3])

    Encoder_Pressure_Preference_Force = ...
        [angle, ref_trj, mappedPressure, Preference, force, auditory_input];
    outlet.push_sample(Encoder_Pressure_Preference_Force);

end
