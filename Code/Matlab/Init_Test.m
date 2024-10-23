clc; clear;

%%
disp('Defining the Arduino...');
a = arduino('COM12','Mega2560','Libraries','rotaryEncoder')
disp('Defining the Encoder...');
encoder = rotaryEncoder(a,'D2','D3',1000)

%% Setting the initial position
resetCount(encoder);

%% Set the variables
time_limit = 60; % Time limit for the experiment

%% Read the data from the encoder
Angle = []; Time_enc = []; time = 0;
figure
ylabel('Knee Angle (deg)')
xlabel('Time (s)')
while (time<time_limit)
    [count,time] = readCount(encoder);
    angle = (360*count)/4000;
    Angle = [Angle;angle]; Time_enc = [Time_enc;time];
    p = plot(Time_enc,Angle,'b',time,angle,'b');
    p(1).LineWidth = 3; p(2).LineWidth = 10;
    p(2).Marker = 'o';
    grid on
end

%% Calculating the Upper and Lower limits for the experiment
[pks_p,locs_p] = findpeaks(Angle); [pks_n,locs_n] = findpeaks(-Angle);
pks_n = -pks_n;
Up_lim = mean(pks_p); Lo_lim = mean(pks_n);

%% Calculating the average frequency
[m,~] = size(locs_p); period = zeros(m-1,1);
for i = 1:m-1
    period(i) = Time_enc(locs_p(i+1,1))-Time_enc(locs_p(i,1));
end
period_avg = mean(period); freq_avg = 1/period_avg;

%% Plotting the data with upper and lower limits
[n,~] = size(Angle);
figure
plot(Angle)
hold on
plot(locs_p,pks_p,'o')
hold on
plot(locs_n,pks_n,'o')
hold on
plot(Up_lim.*ones(n,1))
hold on
plot(Lo_lim.*ones(n,1))

%% Display the results for the upper and lower limits and the frequency
disp(['Upper Limit is ',num2str(Up_lim), '(deg)'])
disp(['Lower Limit is ',num2str(Lo_lim), '(deg)'])
disp(['Average Frequency is ',num2str(freq_avg), '(Hz)'])

%% Saving the workspace variables
save('Initial_Test.mat')