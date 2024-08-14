clc
clear

%% Load the data
% Before running this section you need to modify the paths in the 
% xdf_load_matlab.m function based on your directory

% All signals from all sessions concatenated (it takes time!)
subject_id = 5;
sessions   = 4;

[~, ~, ~, ~, All_Experiment, All_Experiment_time] ...
    = runs_concatenated(subject_id, sessions);



%%


% figure()
% tiledlayout(2,1)
% nexttile
% plot(All_Experiment_time(1,1:end-1), All_Experiment(6, 1:end-1))
% set(gca, 'XLim', [2.710 2.717]*10000)
% nexttile
% plot(All_Experiment_time(1,1:end-1), diff(All_Experiment(6, :)), 'k')
% set(gca, 'XLim', [2.710 2.717]*10000)

%%
start_beep = find(diff(All_Experiment(6, :)) == 1);
start_beep(1:6) = [];
finish_beep = find(diff(All_Experiment(6, :)) == -1);
finish_beep(1:6) = [];



XLimits = [All_Experiment_time(1,start_beep)' ...
    All_Experiment_time(1,finish_beep)'];
[pks_high_peaks, locs_high_peaks] = findpeaks(All_Experiment(1, :));
[pks_low_peaks, locs_low_peaks] = findpeaks((-1)*All_Experiment(1, :));
pks_low_peaks = -pks_low_peaks;

%% find peaks between single and double beep sounds
trial_pks_high_peaks  = [];
trial_locs_high_peaks = [];

trial_pks_low_peaks   = [];
trial_locs_low_peaks  = [];

for i = 1:length(start_beep)
    
    trial_locs_high_peaks = [trial_locs_high_peaks ...
        locs_high_peaks(1, ...
        start_beep(i) <= locs_high_peaks & ...
        locs_high_peaks <= finish_beep(i))];
    trial_pks_high_peaks = [trial_pks_high_peaks ...
        pks_high_peaks(1, ...
        start_beep(i) <= locs_high_peaks & ...
        locs_high_peaks <= finish_beep(i))];

    trial_locs_low_peaks = [trial_locs_low_peaks ...
        locs_low_peaks(1, ...
        start_beep(i) <= locs_low_peaks & ...
        locs_low_peaks <= finish_beep(i))];
    trial_pks_low_peaks = [trial_pks_low_peaks ...
        pks_low_peaks(1, ...
        start_beep(i) <= locs_low_peaks & ...
        locs_low_peaks <= finish_beep(i))];

end



%% Initialize Trials information
Trials_encoder_events = cell(1, size(XLimits,1));
for i = 1:numel(Trials_encoder_events)
    % Trials_encoder_events{1, i}.Pressure = All_Experiment(3, start_beep(1, i));
    % if i ~= numel(Trials_encoder_events)
    %     Trials_encoder_events{1, i}.Score = All_Experiment(4, start_beep(1, i+1));
    % else
    %     Trials_encoder_events{1, i}.Score = All_Experiment(4, end);
    % end
    
    Trials_encoder_events{1, i}.high_peaks = [];
    Trials_encoder_events{1, i}.low_peaks  = [];
end

%%
% figure()
% % plot(All_Experiment_time(1, start_beep(1):finish_beep(1)), ...
% %     All_Experiment(1, start_beep(1):finish_beep(1)))
% plot(All_Experiment_time(1, :)+5, ...
%     All_Experiment(1, :), 'b')
% hold on
% % plot(All_Experiment_time(1, start_beep(1):finish_beep(1)), ...
% %     All_Experiment(2, start_beep(1):finish_beep(1)))
% plot(All_Experiment_time(1, :), ...
%     All_Experiment(2, :), 'k')
% 
% xline(All_Experiment_time(start_beep), 'Color', 'g')
% xline(All_Experiment_time(finish_beep), 'Color', 'r')
% 
% 
% %% 
% plot(All_Experiment_time, All_Experiment(1, :))
% set(gca, 'XLim', XLimits(7,:))
% hold on
% % [pks, locs] = findpeaks(All_Experiment(1, :));
% plot(All_Experiment_time(1, locs), pks, 'Marker', 'v', 'MarkerFaceColor', 'b')
