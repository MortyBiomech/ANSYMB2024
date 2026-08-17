figure()
hold on

plot(All_EMG_time, EMG.data(1, :), 'k');


EMG_epoch_info = EMG_epoched.epoch(matched_vector(i));
eventtypes = EMG_epoch_info.eventtype;
FlxS_type_indx = strcmpi(eventtypes, 'FlxS');
eventlatencies = EMG_epoch_info.eventlatency;
FlxS_latency_indx = cell2mat(eventlatencies) == 0;
FlxS_indx = find(FlxS_type_indx & FlxS_latency_indx);

init_time_lastEMGpoint = ...
    EMG_epoch_info.eventinit_time_lastEMGpoint{FlxS_indx};
init_time_FlxS_EMG = EMG_epoch_info.eventinit_time{FlxS_indx};
delta_t = 1000*(init_time_FlxS_EMG - init_time_lastEMGpoint); % mili-second


ee = 827;

xlim([All_EMG_time(X_full(ee)) - 5, All_EMG_time(X_full(ee)) + 5])

xline(All_EMG_time(X_full(ee)), 'r--')
% xline(all_event_latency_EMG(ee)*0.001 + All_EMG_time(1) + 0.0005, 'k')
xline(all_event_time_AllEEGtime(ee), 'c-.')
% xline(All_EMG_time(floor(3.599407476000000e6)), 'LineWidth', 2, 'Color', 'r')


% all_event_time_AllEEGtime(ee) - All_EMG_time(X_full(ee)) < 0.0005
% All_EMG_time(X_full(ee)+1) - all_event_time_AllEEGtime(ee) < 0.0005











%%
figure()
hold on
plot(0.0005*(1:length(EMG_events.times)), ...
   repmat(1, 1, length(EMG_events.times)), 'b');
plot(0.002*(1:length(EEG_orig.times)), ...
   repmat(2, 1, length(EEG_orig.times)), 'r');
ylim([-1, 4])
