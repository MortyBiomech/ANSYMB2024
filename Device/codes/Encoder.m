clc; clear; close all;
%%
disp('Defining the Arduino...');
a = arduino('COM4','Mega2560','Libraries','rotaryEncoder')
disp('Defining the Encoder...');
encoder = rotaryEncoder(a,'D2','D3',1000)
%% Send data to LSL
% instantiate the library
disp('Loading library...');
lib = lsl_loadlib();
% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Encoder','encoder',1,0,'cf_float32');
disp('Opening an outlet...');
outlet = lsl_outlet(info);
% send data into the outlet, sample by sample
disp('Now transmitting data...');
% %% resolving a stream
% disp('Resolving an EEG stream...');
% result = {};
% while isempty(result)
%     result = lsl_resolve_byprop(lib,'type','EEG'); 
% end
% % create a new inlet
% disp('Opening an inlet...');
% inlet = lsl_inlet(result{1});
%% Setting the initial position
resetCount(encoder);
%% Read the data from the encoder
Angle = []; Time_enc = [];
% EEG = []; Time_EEG = [];
while 1
    [count,time] = readCount(encoder);
    angle = (360*count)/4000;
    Angle = [Angle;angle]; Time_enc = [Time_enc;time];
    plot(Time_enc,Angle)
    outlet.push_sample(angle);
%     %% Reading the data from LSL
%     % get data from the inlet
%     [vec,ts] = inlet.pull_sample();
%     EEG = [EEG;vec]; Time_EEG = [Time_EEG; ts];    
end