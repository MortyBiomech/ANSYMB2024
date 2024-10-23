clc; clear; close all;
%%
a = arduino()
Black = []; White = [];

while 1
%     black = readDigitalPin(a, "D2");
%     white = readDigitalPin(a, "D3");
black = readVoltage(a, "A6");
white = readVoltage(a, "A7");
    Black = [Black; black]; White = [White; white];
    plot(Black)
    hold on
    plot(White)

end