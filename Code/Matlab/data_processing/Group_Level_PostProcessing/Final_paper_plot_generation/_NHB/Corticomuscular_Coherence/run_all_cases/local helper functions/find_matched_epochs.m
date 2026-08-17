function [matched_vector, condition_vector] = find_matched_epochs( ...
        EEG_epoched_main, EMG_epoched)
% FIND_MATCHED_EPOCHS  Replicate the epoch-matching logic from the main script.
%   Returns matched_vector [nEEGEpochs x 1] (EMG epoch index per EEG epoch,
%   NaN if no match) and condition_vector [nEEGEpochs x 1].

    % --- collect FlxS init_time for each EMG epoch ---
    eventtypes_EMG   = {EMG_epoched.epoch.eventtype};
    eventlatency_EMG = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventlatency}, 'UniformOutput', false);

    FlxS_type_lgl = cellfun(@(x) strcmpi(x, 'FlxS'), eventtypes_EMG, ...
        'UniformOutput', false);
    FlxS_lat_lgl  = cellfun(@(x) x == 0, eventlatency_EMG, ...
        'UniformOutput', false);
    FlxS_idx_EMG  = cellfun(@(x, y) find(x & y), ...
        FlxS_type_lgl, FlxS_lat_lgl, 'UniformOutput', false);

    more_one = cellfun(@(x) numel(x) > 1, FlxS_idx_EMG);

    init_time_cell = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventinit_time_EEG}, 'UniformOutput', false);
    init_time_EMG  = cellfun(@(x, y) x(y(1)), init_time_cell, FlxS_idx_EMG, ...
        'UniformOutput', false);
    init_time_EMG  = cell2mat(init_time_EMG);

    N_EEG = length(EEG_epoched_main.epoch);
    matched_vector   = NaN(N_EEG, 1);
    condition_vector = zeros(N_EEG, 1);

    for i = 1:N_EEG
        evtypes = EEG_epoched_main.epoch(i).eventtype;
        evlat   = cell2mat(EEG_epoched_main.epoch(i).eventlatency);
        FlxS_i  = find(strcmpi(evtypes, 'FlxS') & evlat == 0, 1);
        if isempty(FlxS_i), continue; end

        t_EEG = EEG_epoched_main.epoch(i).eventinit_time{FlxS_i};
        hit   = find(abs(init_time_EMG - t_EEG) < 1e-10);

        if numel(hit) > 1 || any(more_one(hit)) || isempty(hit)
            matched_vector(i) = NaN;
        else
            matched_vector(i)   = hit;
            condition_vector(i) = str2double( ...
                EEG_epoched_main.epoch(i).eventcond{FlxS_i});
        end
    end
end