clear; close all; clc;

% get some user settings
ledPin = 'D13';
deltaT_blink = 0.5;

% Use the Matlab support package for Arduino hardware to instentiate an
% Arduino object that will be used to communicate with the Arduino Mega
% board
port = 'COM4';
board = 'Mega2560';

a = arduino(port,board);
for k = 1:10

    a.writeDigitalPin(ledPin, 0);
    pause(deltaT_blink/2);

    a.writeDigitalPin(ledPin, 1);
    pause(deltaT_blink/2);

end