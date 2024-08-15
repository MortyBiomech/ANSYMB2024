function Trials_Info = Main_event_selection(input_streams, ...
                                            EEG, ...
                                            Trials_encoder_events, ...
                                            subject_id, ...
                                            data_path)


    %% Extract data
    All_EEG_time = input_streams.All_EEG_time;
    All_EMG_time = input_streams.All_EMG_time;
    All_Experiment = input_streams.All_Exp;
    All_Experiment_time = input_streams.All_Exp_time;
    
    %% Extract events on Experimnet streams
    start_beep = find(diff(All_Experiment(6, :)) == 1);
    start_beep(1:6) = []; % first six trials were part of the familiarization
    start_beep_time_Expdata = All_Experiment_time(start_beep);
    start_beep_indx_All_EEG = ...
        knnsearch(All_EEG_time', start_beep_time_Expdata');
    
    
    pressure_change_time_on_Expdata = All_Experiment_time(start_beep) - 2;
    pressure_change_indx_on_All_EEG = ...
        knnsearch(All_EEG_time', pressure_change_time_on_Expdata');
    pressure_change_time_on_All_EEG = All_EEG_time(pressure_change_indx_on_All_EEG);
    
    
    finish_beep = find(diff(All_Experiment(6, :)) == -1);
    finish_beep(1:6) = []; % first six trials were part of the familiarization
    finish_beep_time_on_Expdata = All_Experiment_time(finish_beep);
    finish_beep_indx_on_All_EEG = ...
        knnsearch(All_EEG_time', finish_beep_time_on_Expdata');
    finish_beep_time_on_All_EEG = All_EEG_time(finish_beep_indx_on_All_EEG);
    
    
    %% Removing trials which have empty high or low peaks (due to mistake in peak selction!)
    empty_trial_peak_indx = [];
    for i = 1:length(Trials_encoder_events)
        if isempty(Trials_encoder_events{1, i}.high_peaks) || ...
            isempty(Trials_encoder_events{1, i}.low_peaks)
            empty_trial_peak_indx(1, end+1) = i;
        end
    end
    if ~isempty(empty_trial_peak_indx)
        Trials_encoder_events(empty_trial_peak_indx) = [];
    end
    
    %% Initialize the structure template
    structTemplate = ...
        struct('Trial_start_time',[], 'Trial_start_indx', [], ...
        'Trial_end_time',[], 'Trial_end_indx', [], ...
        'case', [], ...
        'flexion_start_time', [], 'flexion_start_indx', [], ...
        'flexion_end_time', [], 'flexion_end_indx', [], ...
        'extension_start_time', [], 'extension_start_indx', [], ...
        'extension_end_time', [], 'extension_end_indx', [], ...
        'flextoflex_start_time', [], 'flextoflex_start_indx', [], ...
        'flextoflex_end_time', [], 'flextoflex_end_indx', []);
    
    events_on_eeg = struct('Raw', structTemplate, ...
        'Preprocessed', structTemplate);
    
    general_info = struct('Pressure', [], 'Score', []);
    
    epochs_events = struct('EEG_stream', events_on_eeg, ...
        'EMG_stream', structTemplate, 'EXP_stream', structTemplate);
    
    final_events_structTemplate = struct('General', general_info, ...
        'Events', epochs_events);
    
    % Create a cell array of size 1xN_trials to store events on all streams
    Trials_Info = repmat({final_events_structTemplate}, ...
        1, length(Trials_encoder_events));
    
    
    %% Loop over all trials to fill the events timings and indeces
    for i = 1:length(Trials_encoder_events)
    
        Trials_Info{1, i}.General.Pressure = All_Experiment(3, start_beep(1, i));
        if i ~= length(Trials_encoder_events)
            Trials_Info{1, i}.General.Score = All_Experiment(4, start_beep(1, i+1));
        else
            Trials_Info{1, i}.General.Score = All_Experiment(4, end);
        end
        
        
        if Trials_encoder_events{1, i}.high_peaks.index(1) > ...
                Trials_encoder_events{1, i}.low_peaks.index(1)
            flag1 = 1;
        else
            flag1 = 0;
        end
    
        if Trials_encoder_events{1, i}.high_peaks.index(end) > ...
                Trials_encoder_events{1, i}.low_peaks.index(end)
            flag2 = 1;
        else
            flag2 = 0;
        end
        
    
    
        all_types = {EEG.event.type};
    
    
        % Find the indices where the type matches 'PC_Pressure_Change'
        pressure_change_indx_on_EEGevent = ...
            find(strcmp(all_types, 'PC_Pressure_Change'), i);
        pressure_change_latency_on_EEGevent = ...
            EEG.event(pressure_change_indx_on_EEGevent(i)).latency;
        pressure_change_time_on_EEGtimes = EEG.times(pressure_change_latency_on_EEGevent);
    
        % Trial start time & index on Preprocessed EEG stream
        Trials_Info{1, i}.Events.EEG_stream.Preprocessed.Trial_start_time = ...
            pressure_change_time_on_EEGtimes;
        Trials_Info{1, i}.Events.EEG_stream.Preprocessed.Trial_start_indx = ...
            pressure_change_latency_on_EEGevent;
    
        % Trial start time & index on Raw EEG stream
        Trials_Info{1, i}.Events.EEG_stream.Raw.Trial_start_time = ...
            pressure_change_time_on_All_EEG(i);
        Trials_Info{1, i}.Events.EEG_stream.Raw.Trial_start_indx = ...
            pressure_change_indx_on_All_EEG(i);
    
        % Trial start time & index on EMG stream
        [~, tmp_indx] = min(abs(All_EMG_time - pressure_change_time_on_All_EEG(i)));
        Trials_Info{1, i}.Events.EMG_stream.Trial_start_indx = tmp_indx;
        Trials_Info{1, i}.Events.EMG_stream.Trial_start_time = All_EMG_time(tmp_indx);
    
        % Trial start time & index on EXP stream
        Trials_Info{1, i}.Events.EXP_stream.Trial_start_indx = start_beep(i);
        Trials_Info{1, i}.Events.EXP_stream.Trial_start_time = All_Experiment_time(start_beep(i));
    
    
    
    
    
        % Find the indices where the type matches 'FB_Finish_Beep'
        finish_beep_indx_on_EEGevent = ...
            find(strcmp(all_types, 'FB_Finish_Beep'), i);
        finish_beep_latency_on_EEGevent = ...
            EEG.event(finish_beep_indx_on_EEGevent(i)).latency;
        finish_beep_time_on_EEGtimes = EEG.times(finish_beep_latency_on_EEGevent);
    
        % Trial end time & index on Preprocessed EEG stream
        Trials_Info{1, i}.Events.EEG_stream.Preprocessed.Trial_end_time = ...
            finish_beep_time_on_EEGtimes;
        Trials_Info{1, i}.Events.EEG_stream.Preprocessed.Trial_end_indx = ...
            finish_beep_latency_on_EEGevent;
    
        % Trial end time & index on Raw EEG stream
        Trials_Info{1, i}.Events.EEG_stream.Raw.Trial_end_time = ...
            finish_beep_time_on_All_EEG(i);
        Trials_Info{1, i}.Events.EEG_stream.Raw.Trial_end_indx = ...
            finish_beep_indx_on_All_EEG(i);
    
        % Trial end time & index on EMG stream
        [~, tmp_indx] = min(abs(All_EMG_time - finish_beep_time_on_All_EEG(i)));
        Trials_Info{1, i}.Events.EMG_stream.Trial_end_indx = tmp_indx;
        Trials_Info{1, i}.Events.EMG_stream.Trial_end_time = All_EMG_time(tmp_indx);
    
        % Trial end time & index on EXP stream
        Trials_Info{1, i}.Events.EXP_stream.Trial_end_indx = finish_beep(i);
        Trials_Info{1, i}.Events.EXP_stream.Trial_end_time = All_Experiment_time(finish_beep(i));
    
    
    
    
        
        % Find the indices where the type matches 'SB_Start_Beep'
        start_beep_indx_on_EEGevent = ...
            find(strcmp(all_types, 'SB_Start_Beep'), i);
        start_beep_latency_on_EEGevent = ...
            EEG.event(start_beep_indx_on_EEGevent(i)).latency;
        start_beep_time_on_EEGtimes = EEG.times(start_beep_latency_on_EEGevent);
    
    
    
    
        if flag1 == 1 && flag2 == 1 % case 1
            Case = 1;
            % filling the case parameter
            Trials_Info{1, i}.Events.EEG_stream.Raw.case = 1;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.case = 1;
            Trials_Info{1, i}.Events.EMG_stream.case = 1;
            Trials_Info{1, i}.Events.EXP_stream.case = 1;
    
            epoch_count = numel(Trials_encoder_events{1, i}.high_peaks.time)-1;
    
        elseif flag1 == 1 && flag2 == 0 % case 2
            Case = 2;
            % filling the case parameter
            Trials_Info{1, i}.Events.EEG_stream.Raw.case = 2;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.case = 2;
            Trials_Info{1, i}.Events.EMG_stream.case = 2;
            Trials_Info{1, i}.Events.EXP_stream.case = 2;
    
            epoch_count = numel(Trials_encoder_events{1, i}.high_peaks.time);
    
        elseif flag1 == 0 && flag2 == 1 % case 3
            Case = 3;
            % filling the case parameter
            Trials_Info{1, i}.Events.EEG_stream.Raw.case = 3;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.case = 3;
            Trials_Info{1, i}.Events.EMG_stream.case = 3;
            Trials_Info{1, i}.Events.EXP_stream.case = 3;
    
            epoch_count = numel(Trials_encoder_events{1, i}.high_peaks.time)-1;
    
        elseif flag1 == 0 && flag2 == 0 % case 4
            Case = 4;
            % filling the case parameter
            Trials_Info{1, i}.Events.EEG_stream.Raw.case = 4;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.case = 4;
            Trials_Info{1, i}.Events.EMG_stream.case = 4;
            Trials_Info{1, i}.Events.EXP_stream.case = 4;
    
            epoch_count = numel(Trials_encoder_events{1, i}.high_peaks.time);
    
        end
    
    
    
        for j = 1:epoch_count
    
            if Case == 1
                % flexion start
                [~, FlxS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flexion end
                [~, FlxE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j+1)));
                % extension start
                [~, ExtS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j+1)));
                % extension end
                [~, ExtE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                % flextoflex start
                [~, flextoflexS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flextoflex end
                [~, flextoflexE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                
    
            elseif Case == 2
                % flexion start
                [~, FlxS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flexion end
                [~, FlxE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j+1)));
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    % extension start
                    [~, ExtS_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j+1)));
                    % extension end
                    [~, ExtE_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                    % flextoflex start
                    [~, flextoflexS_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                    % flextoflex end
                    [~, flextoflexE_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                end
               
    
            elseif Case == 3
                % flexion start
                [~, FlxS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flexion end
                [~, FlxE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j)));
                % extension start
                [~, ExtS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j)));
                % extension end
                [~, ExtE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                % flextoflex start
                [~, flextoflexS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flextoflex end
                [~, flextoflexE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
    
    
            elseif Case == 4
                % flexion start
                [~, FlxS_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                % flexion end
                [~, FlxE_indx_in_All_EEG_time] = ...
                    min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j)));
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    % extension start
                    [~, ExtS_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.low_peaks.time(j)));
                    % extension end
                    [~, ExtE_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                    % flextoflex start
                    [~, flextoflexS_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j)));
                    % flexto flex  end
                    [~, flextoflexE_indx_in_All_EEG_time] = ...
                        min(abs(All_EEG_time - Trials_encoder_events{1, i}.high_peaks.time(j+1)));
                end
            end
    
    
            %% Events on the EXP stream
            % flexion start
            Trials_Info{1, i}.Events.EXP_stream.flexion_start_indx(end+1, 1) = ...
                Trials_encoder_events{1, i}.high_peaks.index(j);
            Trials_Info{1, i}.Events.EXP_stream.flexion_start_time(end+1, 1) = ...
                Trials_encoder_events{1, i}.high_peaks.time(j);
            % flexion end
            if Case == 1 || Case == 2
                Trials_Info{1, i}.Events.EXP_stream.flexion_end_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.index(j+1);
                Trials_Info{1, i}.Events.EXP_stream.flexion_end_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.time(j+1);
            else
                Trials_Info{1, i}.Events.EXP_stream.flexion_end_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.index(j);
                Trials_Info{1, i}.Events.EXP_stream.flexion_end_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.time(j);
            end
            % extension start
            if Case == 1
                Trials_Info{1, i}.Events.EXP_stream.extension_start_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.index(j+1);
                Trials_Info{1, i}.Events.EXP_stream.extension_start_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.time(j+1);
            elseif Case == 2
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    Trials_Info{1, i}.Events.EXP_stream.extension_start_indx(end+1, 1) = ...
                        Trials_encoder_events{1, i}.low_peaks.index(j+1);
                    Trials_Info{1, i}.Events.EXP_stream.extension_start_time(end+1, 1) = ...
                        Trials_encoder_events{1, i}.low_peaks.time(j+1);
                end
            elseif Case == 3
                Trials_Info{1, i}.Events.EXP_stream.extension_start_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.index(j);
                Trials_Info{1, i}.Events.EXP_stream.extension_start_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.low_peaks.time(j);
            elseif Case == 4
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    Trials_Info{1, i}.Events.EXP_stream.extension_start_indx(end+1, 1) = ...
                        Trials_encoder_events{1, i}.low_peaks.index(j);
                    Trials_Info{1, i}.Events.EXP_stream.extension_start_time(end+1, 1) = ...
                        Trials_encoder_events{1, i}.low_peaks.time(j);
                end
            end
            % extension end
            if Case == 1 || Case == 3
                Trials_Info{1, i}.Events.EXP_stream.extension_end_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.index(j+1);
                Trials_Info{1, i}.Events.EXP_stream.extension_end_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.time(j+1);
            else
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    Trials_Info{1, i}.Events.EXP_stream.extension_end_indx(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.index(j+1);
                    Trials_Info{1, i}.Events.EXP_stream.extension_end_time(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.time(j+1);
                end
            end
            % flextoflex start
            if Case == 1 || Case == 3
                Trials_Info{1, i}.Events.EXP_stream.flextoflex_start_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.index(j);
                Trials_Info{1, i}.Events.EXP_stream.flextoflex_start_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.time(j);
            else
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    Trials_Info{1, i}.Events.EXP_stream.flextoflex_start_indx(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.index(j);
                    Trials_Info{1, i}.Events.EXP_stream.flextoflex_start_time(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.time(j);
                end
            end
            % flextoflex end
            if Case == 1 || Case == 3
                Trials_Info{1, i}.Events.EXP_stream.flextoflex_end_indx(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.index(j+1);
                Trials_Info{1, i}.Events.EXP_stream.flextoflex_end_time(end+1, 1) = ...
                    Trials_encoder_events{1, i}.high_peaks.time(j+1);
            else
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    Trials_Info{1, i}.Events.EXP_stream.flextoflex_end_indx(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.index(j+1);
                    Trials_Info{1, i}.Events.EXP_stream.flextoflex_end_time(end+1, 1) = ...
                        Trials_encoder_events{1, i}.high_peaks.time(j+1);
                end
            end
    
            %% Events on EEG stream
            %% Flexion Start
            FlxS_in_All_EEG_time = All_EEG_time(FlxS_indx_in_All_EEG_time);
    
            % evevnts on EEG Raw stream
            Trials_Info{1, i}.Events.EEG_stream.Raw.flexion_start_indx(end+1, 1) = FlxS_indx_in_All_EEG_time;
            Trials_Info{1, i}.Events.EEG_stream.Raw.flexion_start_time(end+1, 1) = FlxS_in_All_EEG_time;
            
    
            t1 = (FlxS_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
            [~, t1_index_temp] = min(abs(EEG.times - t1));
    
            % events on EEG Preprocessed stream
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flexion_start_indx(end+1, 1) = t1_index_temp;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flexion_start_time(end+1, 1) = EEG.times(t1_index_temp);
    
            
    
            % events on EMG stream
            [~, t1_index_temp] = min(abs(All_EMG_time - FlxS_in_All_EEG_time));
            Trials_Info{1, i}.Events.EMG_stream.flexion_start_indx(end+1, 1) = t1_index_temp;
            Trials_Info{1, i}.Events.EMG_stream.flexion_start_time(end+1, 1) = All_EMG_time(t1_index_temp);
            
            
    
            
            %% Flexion End
            FlxE_in_All_EEG_time = All_EEG_time(FlxE_indx_in_All_EEG_time);
    
            % evevnts on EEG Raw stream
            Trials_Info{1, i}.Events.EEG_stream.Raw.flexion_end_indx(end+1, 1) = FlxE_indx_in_All_EEG_time;
            Trials_Info{1, i}.Events.EEG_stream.Raw.flexion_end_time(end+1, 1) = FlxE_in_All_EEG_time;
            
            
            t2 = (FlxE_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
            [~, t2_index_temp] = min(abs(EEG.times - t2));
    
            % events on EEG Preprocessed stream
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flexion_end_indx(end+1, 1) = t2_index_temp;
            Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flexion_end_time(end+1, 1) = EEG.times(t2_index_temp);
    
    
    
            % events on EMG stream
            [~, t2_index_temp] = min(abs(All_EMG_time - FlxE_in_All_EEG_time));
            Trials_Info{1, i}.Events.EMG_stream.flexion_end_indx(end+1, 1) = t2_index_temp;
            Trials_Info{1, i}.Events.EMG_stream.flexion_end_time(end+1, 1) = All_EMG_time(t2_index_temp);
            
           
    
    
            %% Extension Start
            if Case == 2 || Case == 4
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    ExtS_in_All_EEG_time = All_EEG_time(ExtS_indx_in_All_EEG_time);
    
                    % evevnts on EEG Raw stream
                    Trials_Info{1, i}.Events.EEG_stream.Raw.extension_start_indx(end+1, 1) = ExtS_indx_in_All_EEG_time;
                    Trials_Info{1, i}.Events.EEG_stream.Raw.extension_start_time(end+1, 1) = ExtS_in_All_EEG_time;
                
                    
                    t3 = (ExtS_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                    [~, t3_index_temp] = min(abs(EEG.times - t3));
                
                    % events on EEG Preprocessed stream
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_start_indx(end+1, 1) = t3_index_temp;
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_start_time(end+1, 1) = EEG.times(t3_index_temp);
                    
                        
    
                    % events on EMG stream
                    [~, t3_index_temp] = min(abs(All_EMG_time - ExtS_in_All_EEG_time));
                    Trials_Info{1, i}.Events.EMG_stream.extension_start_indx(end+1, 1) = t3_index_temp;
                    Trials_Info{1, i}.Events.EMG_stream.extension_start_time(end+1, 1) = All_EMG_time(t3_index_temp);
                    
                    
    
                end
            else
                ExtS_in_All_EEG_time = All_EEG_time(ExtS_indx_in_All_EEG_time);
    
                % evevnts on EEG Raw stream
                Trials_Info{1, i}.Events.EEG_stream.Raw.extension_start_indx(end+1, 1) = ExtS_indx_in_All_EEG_time;
                Trials_Info{1, i}.Events.EEG_stream.Raw.extension_start_time(end+1, 1) = ExtS_in_All_EEG_time;
            
                
                t3 = (ExtS_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                [~, t3_index_temp] = min(abs(EEG.times - t3));
            
                % events on EEG Preprocessed stream
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_start_indx(end+1, 1) = t3_index_temp;
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_start_time(end+1, 1) = EEG.times(t3_index_temp);
                
    
    
                % events on EMG stream
                [~, t3_index_temp] = min(abs(All_EMG_time - ExtS_in_All_EEG_time));
                Trials_Info{1, i}.Events.EMG_stream.extension_start_indx(end+1, 1) = t3_index_temp;
                Trials_Info{1, i}.Events.EMG_stream.extension_start_time(end+1, 1) = All_EMG_time(t3_index_temp);
                
               
    
            end
    
    
    
            %% Extension End
            if Case == 2 || Case == 4
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    ExtE_in_All_EEG_time = All_EEG_time(ExtE_indx_in_All_EEG_time);
    
                    % evevnts on EEG Raw stream
                    Trials_Info{1, i}.Events.EEG_stream.Raw.extension_end_indx(end+1, 1) = ExtE_indx_in_All_EEG_time;
                    Trials_Info{1, i}.Events.EEG_stream.Raw.extension_end_time(end+1, 1) = ExtE_in_All_EEG_time;
                
                    
                    t4 = (ExtE_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                    [~, t4_index_temp] = min(abs(EEG.times - t4));
                
                    % events on EEG Preprocessed stream
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_end_indx(end+1, 1) = t4_index_temp;
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_end_time(end+1, 1) = EEG.times(t4_index_temp);
                    
    
    
                    % events on EMG stream
                    [~, t4_index_temp] = min(abs(All_EMG_time - ExtE_in_All_EEG_time));
                    Trials_Info{1, i}.Events.EMG_stream.extension_end_indx(end+1, 1) = t4_index_temp;
                    Trials_Info{1, i}.Events.EMG_stream.extension_end_time(end+1, 1) = All_EMG_time(t4_index_temp);
                    
                   
    
                end
            else
                ExtE_in_All_EEG_time = All_EEG_time(ExtE_indx_in_All_EEG_time);
    
                % evevnts on EEG Raw stream
                Trials_Info{1, i}.Events.EEG_stream.Raw.extension_end_indx(end+1, 1) = ExtE_indx_in_All_EEG_time;
                Trials_Info{1, i}.Events.EEG_stream.Raw.extension_end_time(end+1, 1) = ExtE_in_All_EEG_time;
            
                
                t4 = (ExtE_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                [~, t4_index_temp] = min(abs(EEG.times - t4));
            
                % events on EEG Preprocessed stream
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_end_indx(end+1, 1) = t4_index_temp;
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_end_time(end+1, 1) = EEG.times(t4_index_temp);
                
    
    
                % events on EMG stream
                [~, t4_index_temp] = min(abs(All_EMG_time - ExtE_in_All_EEG_time));
                Trials_Info{1, i}.Events.EMG_stream.extension_end_indx(end+1, 1) = t4_index_temp;
                Trials_Info{1, i}.Events.EMG_stream.extension_end_time(end+1, 1) = All_EMG_time(t4_index_temp);
                
                
    
            end
    
    
    
            %% flextoflex Start
            if Case == 2 || Case == 4
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    flextoflexS_in_All_EEG_time = All_EEG_time(flextoflexS_indx_in_All_EEG_time);
    
                    % evevnts on EEG Raw stream
                    Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_start_indx(end+1, 1) = flextoflexS_indx_in_All_EEG_time;
                    Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_start_time(end+1, 1) = flextoflexS_in_All_EEG_time;
                
                    
                    t5 = (flextoflexS_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                    [~, t5_index_temp] = min(abs(EEG.times - t5));
                
                    % events on EEG Preprocessed stream
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_start_indx(end+1, 1) = t5_index_temp;
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_start_time(end+1, 1) = EEG.times(t5_index_temp);
                    
    
    
                    % events on EMG stream
                    [~, t5_index_temp] = min(abs(All_EMG_time - flextoflexS_in_All_EEG_time));
                    Trials_Info{1, i}.Events.EMG_stream.flextoflex_start_indx(end+1, 1) = t5_index_temp;
                    Trials_Info{1, i}.Events.EMG_stream.flextoflex_start_time(end+1, 1) = All_EMG_time(t5_index_temp);
                    
                    
    
                end
            else
                flextoflexS_in_All_EEG_time = All_EEG_time(flextoflexS_indx_in_All_EEG_time);
    
                % evevnts on EEG Raw stream
                Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_start_indx(end+1, 1) = flextoflexS_indx_in_All_EEG_time;
                Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_start_time(end+1, 1) = flextoflexS_in_All_EEG_time;
            
                
                t5 = (flextoflexS_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                [~, t5_index_temp] = min(abs(EEG.times - t5));
            
                % events on EEG Preprocessed stream
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_start_indx(end+1, 1) = t5_index_temp;
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_start_time(end+1, 1) = EEG.times(t5_index_temp);
                
    
    
                % events on EMG stream
                [~, t5_index_temp] = min(abs(All_EMG_time - flextoflexS_in_All_EEG_time));
                Trials_Info{1, i}.Events.EMG_stream.flextoflex_start_indx(end+1, 1) = t5_index_temp;
                Trials_Info{1, i}.Events.EMG_stream.flextoflex_start_time(end+1, 1) = All_EMG_time(t5_index_temp);
                
                
    
            end
    
    
    
            %% flextoflex End
            if Case == 2 || Case == 4
                if j < numel(Trials_encoder_events{1, i}.high_peaks.time)
                    flextoflexE_in_All_EEG_time = All_EEG_time(flextoflexE_indx_in_All_EEG_time);
    
                    % evevnts on EEG Raw stream
                    Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_end_indx(end+1, 1) = flextoflexE_indx_in_All_EEG_time;
                    Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_end_time(end+1, 1) = flextoflexE_in_All_EEG_time;
                
                    
                    t6 = (flextoflexE_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                    [~, t6_index_temp] = min(abs(EEG.times - t6));
                
                    % events on EEG Preprocessed stream
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_end_indx(end+1, 1) = t6_index_temp;
                    Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_end_time(end+1, 1) = EEG.times(t6_index_temp);
                    
    
    
                    % events on EMG stream
                    [~, t6_index_temp] = min(abs(All_EMG_time - flextoflexE_in_All_EEG_time));
                    Trials_Info{1, i}.Events.EMG_stream.flextoflex_end_indx(end+1, 1) = t6_index_temp;
                    Trials_Info{1, i}.Events.EMG_stream.flextoflex_end_time(end+1, 1) = All_EMG_time(t6_index_temp);
                    
                    
    
                end
            else
                flextoflexE_in_All_EEG_time = All_EEG_time(flextoflexE_indx_in_All_EEG_time);
    
                % evevnts on EEG Raw stream
                Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_end_indx(end+1, 1) = flextoflexE_indx_in_All_EEG_time;
                Trials_Info{1, i}.Events.EEG_stream.Raw.flextoflex_end_time(end+1, 1) = flextoflexE_in_All_EEG_time;
            
                
                t6 = (flextoflexE_in_All_EEG_time - All_EEG_time(start_beep_indx_All_EEG(i)))*1000 + start_beep_time_on_EEGtimes;                
                [~, t6_index_temp] = min(abs(EEG.times - t6));
            
                % events on EEG Preprocessed stream
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_end_indx(end+1, 1) = t6_index_temp;
                Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_end_time(end+1, 1) = EEG.times(t6_index_temp);
                
    
    
                % events on EMG stream
                [~, t6_index_temp] = min(abs(All_EMG_time - flextoflexE_in_All_EEG_time));
                Trials_Info{1, i}.Events.EMG_stream.flextoflex_end_indx(end+1, 1) = t6_index_temp;
                Trials_Info{1, i}.Events.EMG_stream.flextoflex_end_time(end+1, 1) = All_EMG_time(t6_index_temp);
                
               
    
            end
    
    
    
        end
    end
    
    
    % %% check timing on streams
    % err1 = Trials_Info{1, 150}.Events.EXP_stream.flexion_start_time  - ...
    %     Trials_Info{1, 150}.Events.EMG_stream.flexion_start_time;
    % err2 = Trials_Info{1, 100}.Events.EXP_stream.flextoflex_start_time  - ...
    %     Trials_Info{1, 100}.Events.EEG_stream.Raw.flextoflex_start_time;
    % err3 = Trials_Info{1, 30}.Events.EMG_stream.extension_end_time  - ...
    %     Trials_Info{1, 30}.Events.EEG_stream.Raw.extension_end_time;
    
    
    %% Save Trials_Info
    save_path = [data_path, '6_Trials_Info_and_Epoched_data\', ...
        'sub-', num2str(subject_id)];
    save(fullfile(save_path, 'Trials_Info.mat'), 'Trials_Info', '-v7.3');

end