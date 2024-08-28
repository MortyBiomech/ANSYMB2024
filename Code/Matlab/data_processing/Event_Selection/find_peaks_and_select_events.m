function [start_beep, finish_beep, XLimits, ...
    pks_high_peaks, locs_high_peaks, pks_low_peaks, locs_low_peaks, ...
    trial_pks_high_peaks, trial_locs_high_peaks, ...
    trial_pks_low_peaks, trial_locs_low_peaks, ...
    N_Trials, ...
    Trials_encoder_events] = ...
    find_peaks_and_select_events(input_streams)


    %% Extract the data
    All_Experiment = input_streams.All_Exp;
    All_Experiment_time = input_streams.All_Exp_time;
    
    
    %% find start and finish beep moments
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

    %% Number of Trials
    N_Trials = size(XLimits, 1);
    
    %% Initialize Trials information
    Trials_encoder_events = cell(1, size(XLimits,1));
    for i = 1:numel(Trials_encoder_events)
        Trials_encoder_events{1, i}.high_peaks = [];
        Trials_encoder_events{1, i}.low_peaks  = [];
    end

end