% 1- angle, 
% 2- ref_trj, 
% 3- mappedPressure, 
% 4- Preference, 
% 5- force, 
% 6- auditory_input, 
% 7- score_press

%% checking the force sensor data
All_Experiment = output.All_Exp;
All_Experiment_time = output.All_Exp_time;

% Load Trials_encoder_events manually

% for subject below 10 load Trials_Info and Trials_encoder_events from the
% correct folder

%% collect start_move events
% diff_AUD = diff(All_Experiment(6, :));
start_move = [];
end_move   = [];
flexion_starts = [];
extension_starts = [];
for i = 1:length(Trials_encoder_events)
    
    % a = find(diff(All_Experiment(6, :)) == 1, 2*i-1);
    a = find(diff(All_Experiment(6, :)) == 1, i);
    start_move = cat(2, start_move, a(end));
    % b = find(diff(All_Experiment(6, :)) == -2, i);
    b = find(diff(All_Experiment(6, :)) == -1, i);
    end_move   = cat(2, end_move, b(end));

    % flexion_starts = cat(2, flexion_starts, Trials_encoder_events{1, i}.Flexion_Start);
    % extension_starts = cat(2, extension_starts, Trials_encoder_events{1, i}.Extension_Start);

    % for subjects 5 to 9
    flexion_starts = cat(2, flexion_starts, Trials_Info{1, i}.Events.EXP_stream.flexion_start_indx');
    extension_starts = cat(2, extension_starts, Trials_Info{1, i}.Events.EXP_stream.extension_start_indx');

end

%% plot data to look at the force sensor
figure()
plot(All_Experiment_time, All_Experiment(5,:))
hold on

% plot(All_Experiment_time - All_Experiment_time(1), All_Experiment(5,2870), ...
%     'Marker','o', 'MarkerFaceColor','r')

xline(All_Experiment_time(start_move), 'Color', 'b', 'LineStyle', '--')
xline(All_Experiment_time(end_move), 'Color', 'r', 'LineStyle', '--')

plot(All_Experiment_time(flexion_starts), All_Experiment(5, flexion_starts), ...
    'Marker', 'v', 'MarkerFaceColor', [0,100,0]/255, 'MarkerEdgeColor', 'none','LineStyle', 'none')
plot(All_Experiment_time(extension_starts), All_Experiment(5, extension_starts), ...
    'Marker', '^', 'MarkerFaceColor', [139,0,0]/255, 'MarkerEdgeColor', 'none', 'LineStyle', 'none')


%%
trial = 144;
xlim([All_Experiment_time(start_move(trial))-1 All_Experiment_time(end_move(trial))+1])
min_y = min(All_Experiment(5, start_move(trial):end_move(trial)));
max_y = max(All_Experiment(5, start_move(trial):end_move(trial)));
ylim([0.999*min_y, 1.001*max_y])

ylabel('Force Sensor Voltage [a.u.]')
xlabel('Time [s]')
% title(['Subject ', num2str(subject_id), ...
%     ', Pressure ', num2str(Trials_encoder_events{1, trial}.Pressure), ...
%     ', trial ', num2str(trial)])

% for subjects 5 to 9
title(['Subject ', num2str(subject_id), ...
    ', Pressure ', num2str(Trials_Info{1, trial}.General.Pressure), ...
    ', trial ', num2str(trial)])

%%
find(diff(All_Experiment(6,:)) == 1, 1)