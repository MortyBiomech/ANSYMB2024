function [EMG, removeindices] = import_same_event_on_EMG_stream(EMG, EEG, output, ...
    subject_id, processing_path, data_path)    

    % calling from run_CMC_allsubjects.m
    % [EMG_ev, ~]   = import_same_event_on_EMG_stream( ...
    %     EMG_tmpl, EEG_orig, output_timing, ...
    %     subject_id, processing_path, data_path);

    Trials_Encoder_Events_path = [data_path, ...
        '6_0_Trials_Info_and_Events\sub-', num2str(subject_id), '\'];

    
    All_Experiment = output.All_Exp;
    All_Experiment_time = output.All_Exp_time;
    All_EEG_time = output.All_EEG_time;
    All_EMG_time = output.All_EMG_time;

    

    %% Define and Add events
    % Compute latency values
    
    if subject_id > 9
        % start_beep event (first single beep)
        start_beep = find(diff(All_Experiment(6, :)) == 1);
        start_beep = reshape(start_beep, 2, []);
        start_beep_time_Expdata = All_Experiment_time(start_beep(1,:));
        start_beep_indx_EEG = ...
            knnsearch(All_EEG_time', start_beep_time_Expdata');
        start_beep_time_AllEEGtime = All_EEG_time(start_beep_indx_EEG);
        X_full = sum(All_EMG_time.' <= start_beep_time_AllEEGtime, 1);
        X_frac = start_beep_time_AllEEGtime - All_EMG_time(X_full); % second 
        start_beep_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        start_beep_latency_lastEMGpoint = EMG.times(X_full);
        start_beep_latency_EEG = EEG.times(start_beep_indx_EEG);
    else
        % start_beep event (first single beep)
        start_beep = find(diff(All_Experiment(6, :)) == 1);
        start_beep_time_Expdata = All_Experiment_time(start_beep(1,:));
        start_beep_indx_EEG = ...
            knnsearch(All_EEG_time', start_beep_time_Expdata');
        start_beep_time_AllEEGtime = All_EEG_time(start_beep_indx_EEG);
        X_full = sum(All_EMG_time.' <= start_beep_time_AllEEGtime, 1);
        X_frac = start_beep_time_AllEEGtime - All_EMG_time(X_full); % second 
        start_beep_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        start_beep_latency_lastEMGpoint = EMG.times(X_full);
        start_beep_latency_EEG = EEG.times(start_beep_indx_EEG);
    end


  
    if subject_id > 9
        % pressure_change event (2s after first single beep)
        pressure_change_time_Expdata = All_Experiment_time(start_beep(1,:)) + 2;
        pressure_change_indx_EEG = ...
            knnsearch(All_EEG_time', pressure_change_time_Expdata');
        pressure_change_time_AllEEGtime = All_EEG_time(pressure_change_indx_EEG);
        X_full = sum(All_EMG_time.' <= pressure_change_time_AllEEGtime, 1);
        X_frac = pressure_change_time_AllEEGtime - All_EMG_time(X_full); % second 
        pressure_change_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        pressure_change_latency_lastEMGpoint = EMG.times(X_full);
        pressure_change_latency_EEG = EEG.times(pressure_change_indx_EEG);
    else
        % pressure_change event (2s before single beep)
        pressure_change_time_Expdata = All_Experiment_time(start_beep(1,:)) - 2;
        pressure_change_indx_EEG = ...
            knnsearch(All_EEG_time', pressure_change_time_Expdata');
        pressure_change_time_AllEEGtime = All_EEG_time(pressure_change_indx_EEG);
        X_full = sum(All_EMG_time.' <= pressure_change_time_AllEEGtime, 1);
        X_frac = pressure_change_time_AllEEGtime - All_EMG_time(X_full); % second 
        pressure_change_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        pressure_change_latency_lastEMGpoint = EMG.times(X_full);
        pressure_change_latency_EEG = EEG.times(pressure_change_indx_EEG);
    end


    if subject_id > 9
        % start_move event (second single beep, 2s after pressure change)
        start_move = find(diff(All_Experiment(6, :)) == 1);
        start_move = reshape(start_move, 2, []);
        start_move_time_Expdata = All_Experiment_time(start_move(2,:));
        start_move_indx_EEG = ...
            knnsearch(All_EEG_time', start_move_time_Expdata');
        start_move_time_AllEEGtime = All_EEG_time(start_move_indx_EEG);
        X_full = sum(All_EMG_time.' <= start_move_time_AllEEGtime, 1);
        X_frac = start_move_time_AllEEGtime - All_EMG_time(X_full); % second 
        start_move_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        start_move_latency_lastEMGpoint = EMG.times(X_full);
        start_move_latency_EEG = EEG.times(start_move_indx_EEG);
    else
        % start_move event (second single beep, 2s after pressure change)
        start_move = find(diff(All_Experiment(6, :)) == 1);
        start_move_time_Expdata = All_Experiment_time(start_move(1,:));
        start_move_indx_EEG = ...
            knnsearch(All_EEG_time', start_move_time_Expdata');
        start_move_time_AllEEGtime = All_EEG_time(start_move_indx_EEG);
        X_full = sum(All_EMG_time.' <= start_move_time_AllEEGtime, 1);
        X_frac = start_move_time_AllEEGtime - All_EMG_time(X_full); % second 
        start_move_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        start_move_latency_lastEMGpoint = EMG.times(X_full);
        start_move_latency_EEG = EEG.times(start_move_indx_EEG);
    end
    

    if subject_id > 9
        % finish_beep event (20s after start_move event, double beep to stop movement)
        finish_beep = find(diff(All_Experiment(6, :)) == -2);
        finish_beep_time_Expdata = All_Experiment_time(finish_beep);
        finish_beep_indx_EEG = ...
            knnsearch(All_EEG_time', finish_beep_time_Expdata');
        finish_beep_time_AllEEGtime = All_EEG_time(finish_beep_indx_EEG);
        X_full = sum(All_EMG_time.' <= finish_beep_time_AllEEGtime, 1);
        X_frac = finish_beep_time_AllEEGtime - All_EMG_time(X_full); % second 
        finish_beep_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        finish_beep_latency_lastEMGpoint = EMG.times(X_full);
        finish_beep_latency_EEG = EEG.times(finish_beep_indx_EEG);
    else
        % finish_beep event (20s after start_move event, double beep to stop movement)
        finish_beep = find(diff(All_Experiment(6, :)) == -1);
        finish_beep_time_Expdata = All_Experiment_time(finish_beep);
        finish_beep_indx_EEG = ...
            knnsearch(All_EEG_time', finish_beep_time_Expdata');
        finish_beep_time_AllEEGtime = All_EEG_time(finish_beep_indx_EEG);
        X_full = sum(All_EMG_time.' <= finish_beep_time_AllEEGtime, 1);
        X_frac = finish_beep_time_AllEEGtime - All_EMG_time(X_full); % second 
        finish_beep_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        finish_beep_latency_lastEMGpoint = EMG.times(X_full);
        finish_beep_latency_EEG = EEG.times(finish_beep_indx_EEG);
    end
    

    if subject_id > 9
        % score_press event (experimenter presses the scores immidiately after subjects evaluate the task)
        score_press = find(diff(All_Experiment(7, :)) > 0);
        score_press_time_Expdata = All_Experiment_time(score_press);
        score_press_indx_EEG = ...
            knnsearch(All_EEG_time', score_press_time_Expdata');
        score_press_time_AllEEGtime = All_EEG_time(score_press_indx_EEG);
        X_full = sum(All_EMG_time.' <= score_press_time_AllEEGtime, 1);
        X_frac = score_press_time_AllEEGtime - All_EMG_time(X_full); % second 
        score_press_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        score_press_latency_lastEMGpoint = EMG.times(X_full);
        score_press_latency_EEG = EEG.times(score_press_indx_EEG);
    else
        score_press = find(diff(All_Experiment(6, :)) == -1) + 1;
        score_press_time_Expdata = All_Experiment_time(score_press);
        score_press_indx_EEG = ...
            knnsearch(All_EEG_time', score_press_time_Expdata');
        score_press_time_AllEEGtime = All_EEG_time(score_press_indx_EEG);
        X_full = sum(All_EMG_time.' <= score_press_time_AllEEGtime, 1);
        X_frac = score_press_time_AllEEGtime - All_EMG_time(X_full); % second 
        score_press_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
        score_press_latency_lastEMGpoint = EMG.times(X_full);
        score_press_latency_EEG = EEG.times(score_press_indx_EEG);
    end

    
    % load Trials_Encoder_Events
    fileName = [Trials_Encoder_Events_path, 'sub-', num2str(subject_id), ...
        '_Trials_encoder_events.mat'];
    Trials_encoder_events = load(fileName);
    Trials = Trials_encoder_events.Trials_encoder_events;



    % make events description
    desc1 = cell(1, numel(pressure_change_latency_EMG)); % {P(i-1), P(i), Trial}
    for i = 1:length(desc1)
        if i~=1
            desc1{1, i} = {Trials{1, i-1}.Pressure, Trials{1, i}.Pressure, i};
        else
            indx_temp = knnsearch(All_Experiment_time', pressure_change_time_Expdata(1) - 1);
            desc1{1, 1} = {All_Experiment(3, indx_temp), Trials{1, i}.Pressure, i};
        end
    end



    % Add start beep events
    type = repmat({'SB_Start_Beep'}, 1, size(start_beep,2));
    latency = start_beep_latency_EMG;
    latency_EEG = start_beep_latency_EEG;
    latency_lastEMGpoint = start_beep_latency_lastEMGpoint;
    desc = desc1;

    % Add pressure change events
    type = cat(2, type, repmat({'PC_Pressure_Change'}, 1, numel(pressure_change_latency_EMG)));
    latency = cat(2, latency, pressure_change_latency_EMG);
    latency_EEG = cat(2, latency_EEG, pressure_change_latency_EEG);
    latency_lastEMGpoint = ...
        cat(2, latency_lastEMGpoint, pressure_change_latency_lastEMGpoint);
    desc = cat(2, desc, desc1);

    % Add start move events
    type = cat(2, type, repmat({'SM_Start_Move'}, 1, size(start_beep,2)));
    latency = cat(2, latency, start_move_latency_EMG);
    latency_EEG = cat(2, latency_EEG, start_move_latency_EEG);
    latency_lastEMGpoint = ...
        cat(2, latency_lastEMGpoint, start_move_latency_lastEMGpoint);
    desc = cat(2, desc, desc1);
    
    % Add finish beep events
    type = cat(2, type, repmat({'FB_Finish_Beep'}, 1, numel(finish_beep)));
    latency = cat(2, latency, finish_beep_latency_EMG);
    latency_EEG = cat(2, latency_EEG, finish_beep_latency_EEG);
    latency_lastEMGpoint = ...
        cat(2, latency_lastEMGpoint, finish_beep_latency_lastEMGpoint);
    desc = cat(2, desc, desc1);

    % New description including the scores of current and previous trials
    % {P(i-1), P(i), S(i-1), S(i), Trial}
    desc2 = cell(1, numel(Trials)); 
    for i = 1:length(Trials)
        if i~=1
            desc2{1, i} = ...
                {Trials{1, i-1}.Pressure, Trials{1, i}.Pressure, ...
                Trials{1, i-1}.Score, Trials{1, i}.Score, i};
        else
            indx_temp = knnsearch(All_Experiment_time', ...
                pressure_change_time_Expdata(1) - 1);
            desc2{1, 1} = ...
                {All_Experiment(3, indx_temp), Trials{1, i}.Pressure, ...
                All_Experiment(4, indx_temp), Trials{1, i}.Score, i};
        end
    end

    % Add score press events
    type = cat(2, type, repmat({'SP_Score_Press'}, 1, numel(score_press_latency_EMG)));
    latency = cat(2, latency, score_press_latency_EMG);
    latency_EEG = cat(2, latency_EEG, score_press_latency_EEG);
    latency_lastEMGpoint = ...
        cat(2, latency_lastEMGpoint, score_press_latency_lastEMGpoint);
    desc = cat(2, desc, desc2);

    if subject_id > 9
        % Add Trial Start events
        type = cat(2, type, repmat({'TS_Trial_Start'}, 1, numel(start_beep_latency_EMG)));
        latency = cat(2, latency, start_beep_latency_EMG - 2); % 2ms before (one sample in EEG time domain)
        latency_EEG = cat(2, latency_EEG, start_beep_latency_EEG - 2);
        latency_lastEMGpoint = ...
            cat(2, latency_lastEMGpoint, start_beep_latency_lastEMGpoint - 2);
        desc = cat(2, desc, desc1);
    else
        % Add Trial Start events
        type = cat(2, type, repmat({'TS_Trial_Start'}, 1, numel(pressure_change_latency_EMG)));
        latency = cat(2, latency, pressure_change_latency_EMG - 2); % 2ms before (one sample in EEG time domain)
        latency_EEG = cat(2, latency_EEG, pressure_change_latency_EEG - 2);
        latency_lastEMGpoint = ...
            cat(2, latency_lastEMGpoint, pressure_change_latency_lastEMGpoint - 2);
        desc = cat(2, desc, desc1);
    end

    % Add Trial End events
    type = cat(2, type, repmat({'TE_Trial_End'}, 1, numel(score_press_latency_EMG)));
    latency = cat(2, latency, score_press_latency_EMG + 2); % 2ms after (one sample in EEG time domain)
    latency_EEG = cat(2, latency_EEG, score_press_latency_EEG + 2);
    latency_lastEMGpoint = ...
        cat(2, latency_lastEMGpoint, score_press_latency_lastEMGpoint + 2);
    desc = cat(2, desc, desc1);



    %% Write the TS_SB_PC_SM_FB_SP_TE_event.txt file
    % TS: Trial Start
    % SB: Start Beep
    % after 2s 
    % PC: Pressure Change
    % after 2s
    % SM: Start Move
    % after 20s
    % FB: Finish Beep
    % SP: Score Press
    % TE: Trial End

    folder = [processing_path, 'Events', filesep, 'EMG', ...
        filesep, 'sub-', num2str(subject_id)];

    % Ensure the folder exists, if not, create it
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    % File name
    filename = fullfile(folder, 'events_basic.txt');
    
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Cannot open file for writing: %s', filename);
    end
    
    % Write the header
    fprintf(fileID, 'type\tlatency\tdesc\tlatencyEEG\tlatencylastEMGpoint\n');
    
    % Write the data
    for i = 1:7*length(Trials) % each Trial has seven main events
        % Convert the nested cell array in desc to a string with underline separator
        desc_str = strjoin(cellfun(@num2str, desc{i}, 'UniformOutput', false), '_');
        fprintf(fileID, '%s\t%f\t%s\t%f\t%f\n', ...
            type{i}, latency(i), desc_str, ...
            latency_EEG(i)*0.001, latency_lastEMGpoint(i)*0.001);
    end
    
    % Close the file
    fclose(fileID);
    
    % Notify the user
    fprintf('File saved successfully: \n%s\n', filename);


    %% Add Events to the EEG file 
    [EMG, ~] = pop_importevent(EMG, 'event', ...
              filename, 'fields', ...
              {'type', 'latency', 'desc', 'init_time_EEG', 'init_time_lastEMGpoint'}, ...
              'append', 'no', 'align', NaN, 'skipline', 1, 'timeunit', 1E-3);

    % EEG = EMG;
    % eeglab redraw

    %% mark non-experimental segments for removal
    % it is stongly recommended to remove these segments because they may 
    % contain strong artifacts that confuse channel detection and AMICA
    removeindices = zeros(size(start_beep,2)+1 ,2);
    % remove from start to first event
    removeindices(1, :) = [0 EMG.event(1).latency-4]; 

    % add more removeIndices here for pauses or itnerruptions of the 
    % experiment if they have markers or you know their indices in the data

    start_beep_EEGLABevent = find(strcmp({EMG.event.type}, 'SB_Start_Beep'));
    for i = 1:size(start_beep_EEGLABevent,2)-1
        removeindices(i+1, :) = [EMG.event(7*i).latency+4 EMG.event(7*i+1).latency-4]; % this should be changed based on the subject
    end
    removeindices(end, :) = [EMG.event(end).latency+4 EMG.pnts]; % remove from last event to the end
    


    %% Add Flexion Extension Start/End events to the text file 
    % to import them in EEG data (EEGLAB-like structure)
    % Initialize variables
    
    n_trials = length(Trials);
    event_types = {'FlxS', 'FlxE', 'ExtS', 'ExtE'};
    event_fields = {'Flexion_Start', 'Flexion_End', 'Extension_Start', 'Extension_End'};
    
    % Pre-calculate total number of events to preallocate arrays
    total_events = 0;
    for i = 1:n_trials
        for j = 1:numel(event_fields)
            total_events = total_events + numel(Trials{1, i}.(event_fields{j}));
        end
    end
    
    % Preallocate arrays for efficiency
    all_event_times_Expdata = zeros(total_events, 1);
    all_event_labels = cell(total_events, 1);
    all_descs = cell(total_events, 1); 
    trial_indices = zeros(total_events, 1);
    
    idx = 1;
    for i = 1:n_trials
        for j = 1:numel(event_fields)
            event_field = event_fields{j};
            event_label = event_types{j};
            
            event_times = All_Experiment_time(Trials{1, i}.(event_field));
            n_events = numel(event_times);
            
            if n_events > 0
                range = idx:(idx + n_events - 1);
                all_event_times_Expdata(range) = event_times';
                all_event_labels(range) = repmat({event_label}, n_events, 1);
                all_descs(range) = repmat({desc1{1, i}}, n_events, 1);
                trial_indices(range) = i;
                idx = idx + n_events;
            end
        end
    end
    
    % Trim arrays to actual number of events
    all_event_times_Expdata = all_event_times_Expdata(1:idx-1);
    all_event_labels = all_event_labels(1:idx-1);
    all_descs = all_descs(1:idx-1);
    trial_indices = trial_indices(1:idx-1);
    
    % Map event times to EEG indices using interp1 for efficiency
    all_event_indx_EEG = interp1(All_EEG_time, 1:length(All_EEG_time), all_event_times_Expdata, 'nearest', 'extrap');
    
    all_event_time_AllEEGtime = All_EEG_time(all_event_indx_EEG);
    disp('calculating the EMG_time samples for Flx/Ext events (X_full) ...')
    X_full = arrayfun(@(t) find(All_EMG_time <= t, 1, 'last'), ...
        all_event_time_AllEEGtime);
    X_frac = all_event_time_AllEEGtime - All_EMG_time(X_full); % second 
    all_event_latency_EMG = EMG.times(X_full) + 1000*X_frac; % time unit: milisecond
    
    all_event_latency_lastEMGpoint = EMG.times(X_full); 

    all_event_latency_EEG = EEG.times(all_event_indx_EEG);

    % Let's find the gaps
    gap_indx_between_lastEMGpoint_AllEEGtime = X_frac >= 0.0005; 
    delta_t_lastEMPpoint_onesamplebefore = ...
        All_EMG_time(X_full) - All_EMG_time(X_full - 1);
    gap_indx_between_lastEMGpoint_prvsLastEMGpoint = ...
        delta_t_lastEMPpoint_onesamplebefore >= 0.00051;
    ignore_these_events = or(gap_indx_between_lastEMGpoint_AllEEGtime, ...
        gap_indx_between_lastEMGpoint_prvsLastEMGpoint);



    % Organize data back into per-trial cell arrays
    type_extended = cell(1, n_trials);
    latency_extended = cell(1, n_trials);
    desc_extended = cell(1, n_trials);
    latency_extended_lastEMGpoint = cell(1, n_trials);
    latency_extended_EEG = cell(1, n_trials);
    
    for i = 1:n_trials
        idx_trial = (trial_indices == i);
        type_extended{i} = all_event_labels(idx_trial)';
        latency_extended{i} = all_event_latency_EMG(idx_trial)';
        desc_extended{i} = all_descs(idx_trial)';
        latency_extended_lastEMGpoint{i} = ...
            all_event_latency_lastEMGpoint(idx_trial)';
        latency_extended_EEG{i} = all_event_latency_EEG(idx_trial)';
    end



    %% Appending to the previous event file and creating a new one
    % Define the file paths
    original_filename = [processing_path, 'Events', filesep, 'EMG', ...
        filesep, 'sub-', num2str(subject_id), filesep, 'events_basic.txt']; 
    new_filename = [processing_path, 'Events', filesep, 'EMG', ...
        filesep, 'sub-', num2str(subject_id), filesep, 'events_with_FlxExt.txt']; 
    
    
    % Read the existing content of the original file
    existing_lines = {};
    fid = fopen(original_filename, 'r');
    if fid ~= -1
        while ~feof(fid)
            line = fgetl(fid);
            existing_lines{end + 1} = line; %#ok<SAGROW> % Add each line to the cell array
        end
        fclose(fid);
    else
        error('Failed to open the original file for reading.');
    end


    % Ensure data arrays are column vectors
    event_labels = all_event_labels;          % Cell array of strings, size (9730 x 1)
    event_latency = all_event_latency_EMG';   % Numeric array, transpose to (9730 x 1)
    event_descs = all_descs;                  % Cell array of cell arrays, size (9730 x 1)
    event_latency_lastEMGpoint = 0.001*all_event_latency_lastEMGpoint';
    event_latency_EEG = 0.001*all_event_latency_EEG';


    % ignore events which were identified in the gaps (data loss included)
    event_labels = event_labels(~ignore_these_events);
    event_latency = event_latency(~ignore_these_events);
    event_descs = event_descs(~ignore_these_events);
    event_latency_lastEMGpoint = event_latency_lastEMGpoint(~ignore_these_events);
    event_latency_EEG = event_latency_EEG(~ignore_these_events);


    % Process event_descs to create a cell array of strings
    event_descs_str = cell(size(event_descs));  % Initialize cell array for descriptions
    for i = 1:length(event_descs)
        % Extract the cell array of numbers for the current event
        desc_numbers = event_descs{i};  % This is a 1x3 cell array of numbers
    
        % Convert each number to a string
        desc_strings = cellfun(@num2str, desc_numbers, 'UniformOutput', false);
    
        % Combine the strings with underscores
        desc_str = strjoin(desc_strings, '_');
    
        % Store the resulting string
        event_descs_str{i} = desc_str;
    end


    % Append the new lines to the existing lines
    for i = 1:length(event_latency)
        new_line = sprintf('%s\t%f\t%s\t%f\t%f', ...
            event_labels{i}, event_latency(i), ...
            event_descs_str{i}, event_latency_EEG(i), ...
            event_latency_lastEMGpoint(i));
        existing_lines{end + 1} = new_line; %#ok<SAGROW> % Add each new line to the cell array
    end
    
    % Write the combined content to a new file
    fid = fopen(new_filename, 'w');
    if fid == -1
        error('Failed to open the new file for writing.');
    end
    
    for i = 1:length(existing_lines)
        fprintf(fid, '%s\n', existing_lines{i});
    end
    
    % Close the new file
    fclose(fid);



    %% Add Events with Flexion & Extension to the EEG file 
    [EMG, ~] = pop_importevent(EMG, 'event', ...
              new_filename, 'fields', ...
              {'type', 'latency', 'desc', 'init_time_EEG', 'init_time_lastEMGpoint'}, ...
              'append', 'no', 'align', NaN, 'skipline', 1, 'timeunit', 1E-3);

    % EMG = eeg_checkset(EMG);
    % EEG = EMG;
    % eeglab redraw


    

end