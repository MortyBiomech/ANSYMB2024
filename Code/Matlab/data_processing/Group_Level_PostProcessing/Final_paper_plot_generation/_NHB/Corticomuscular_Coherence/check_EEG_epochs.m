EEG_epochs_all = EEG_epoched_main.epoch;

EEG_epoch_id = 1;
EEG_epoch_info = EEG_epochs_all(EEG_epoch_id);

FlxS_type_indx = strcmpi(EEG_epoch_info.eventtype, 'FlxS');
FlxS_latency_indx = cell2mat(EEG_epoch_info.eventlatency) == 0;
FlxS_indx = find(FlxS_type_indx & FlxS_latency_indx);