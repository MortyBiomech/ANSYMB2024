clc
clear

%% Change these paths with respect to your system
study_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path = 'C:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
addpath(genpath('C:\Morteza\MyProjects\ANSYMB2024'))
addpath(genpath('C:\Morteza\LSL\xdf-Matlab-master'))

%% Add important paths
addpath('C:\Morteza\Toolboxes\Fieldtrip\fieldtrip-20231127')
addpath('C:\Morteza\Toolboxes\Fieldtrip\fieldtrip-20231127\fileio')

%% All signals from all sessions concatenated (it takes time!)
subject_id = 18;
rawdata_path = [study_path, '0_source_data\'];
output = runs_concatenated(subject_id, rawdata_path);

All_EEG_time = output.All_EEG_time;
All_Experiment = output.All_Exp;
All_Experiment_time = output.All_Exp_time;


%% Initialize EEGLAB 
if ~exist('ALLCOM','var')
	eeglab;
end

% initialize fieldtrip without adding alternative files to path 
% assuming FT is on your path already or is added via EEGlab plugin manager
global ft_default
ft_default.toolbox.signal = 'matlab';  % can be 'compat' or 'matlab'
ft_default.toolbox.stats  = 'matlab';
ft_default.toolbox.image  = 'matlab';
ft_defaults % this sets up the FieldTrip path


%% [OPTIONAL] check the .xdf data to explore the structure
ftPath = fileparts(which('ft_defaults')); 
addpath(fullfile(ftPath, 'external','xdf')); 
xdfPath = [rawdata_path, 'sub-', num2str(subject_id), filesep, ...
    'ses-S001\eeg\sub-', num2str(subject_id), ...
    '_ses-S001_task-Default_run-001_eeg.xdf']; % enter full path to your .xdf file 

% load .xdf data to check what is in there
streams         = load_xdf(xdfPath);
streamnames     = cellfun(@(x) x.info.name, streams, 'UniformOutput', 0)' % will display names of streams contained in .xdf

% display names of all channels in the .xdf data
for Si = 1:numel(streamnames)
    if isfield( streams{Si}.info.desc, 'channels')
        channelnames    = cellfun(@(x) x.label, streams{Si}.info.desc.channels.channel, 'UniformOutput', 0)'
    end
end


%% [OPTIONAL] enter metadata about the data set, data modalities, and participants

% information about the eeg recording system 
% will be saved in BIDS-folder/sub-XX/[ses-XX]/eeg/*_eeg.json and *_coordsystem.json
% see "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=MAY%20be%20specified.-,Sidecar%20JSON%20(*_eeg.json),-Generic%20fields%20MUST"
% and "https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/03-electroencephalography.html#:~:text=after%20the%20recording.-,Coordinate%20System%20JSON,-(*_coordsystem.json"
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
eegInfo     = [];
% eegInfo.coordsystem.EEGCoordinateSystem     = 'enter the name of your coordinate system'; % only needed when you share eloc
% eegInfo.coordsystem.EEGCoordinateUnits      = 'enter the unit of your coordinate system'; % only needed when you share eloc
% eegInfo.coordsystem.EEGCoordinateSystemDescription = 'enter description of your coordinate system'; % only needed when you share eloc
eegInfo.eeg.SamplingFrequency               = 500; % nominal sampling frequency                            


%% Import data

% iterate over sessions to import xdf files one by one
% full path to your study folder  
studyFolder   = study_path(1:end-1);     

% replace with names of your sessions (if there are no multiple sessions, 
% remove confg.ses and the session loop in the following)
sessionNames  = {'S001', 'S002', 'S003', 'S004'};             

% loop over sessions 
for session = 1:length(sessionNames)
    
    % reset for each loop
    config = [];  
    % required, replace with the folder where you want to store your bids data
    config.bids_target_folder = [study_path, '1_BIDS_data'];  
    
    % required, replace with your xdf file full path
    config.filename  = ...
        fullfile([rawdata_path,'sub-', num2str(subject_id),'\ses-', ...
        sessionNames{session}, '\eeg\sub-', num2str(subject_id), ...
        '_ses-', sessionNames{session}, '_task-Default_run-001_eeg.xdf']);    
    
    % optional, if you have electrode location file, replace with the full 
    % path to the file.
    % config.eeg.chanloc = fullfile([processing_path ,'chanlocs.ced']);
    
    % optional, replace with your task name
    config.task                   = 'KneeSwingingwithExo';  
    config.subject                = subject_id;  % required
    config.session                = sessionNames{session}; % optional
    config.overwrite              = 'on'; % optional
    
    % required, replace with the unique keyword in your eeg stream in 
    % the .xdf file
    config.eeg.stream_name        = 'LiveAmpSN-102108-1125'; 
    
    %------------------------------------------------------------------
    
    % config.motion.streams{1}.xdfname            = 'YourStreamNameInXDF'; % replace with name of the stream corresponding to the first tracking system
    % config.motion.streams{1}.bidsname           = tracking_systems{1};  % a comprehensible name to represent the tracking system 
    % config.motion.streams{1}.tracked_points     = 'headRigid';          % name of the point that is being tracked in the tracking system, the keyword has to be containted in the channel name (see "bemobil_bids_motionconvert")
    % config.motion.streams{1}.tracked_points_anat= 'head';               % example of how the tracked point can be renamed to body part name for metadata
    % 
    % % names of position and quaternion channels in each stream
    % config.motion.streams{1}.positions.channel_names    = {'headRigid_Rigid_headRigid_X';  'headRigid_Rigid_headRigid_Y' ; 'headRigid_Rigid_headRigid_Z' };
    % config.motion.streams{1}.quaternions.channel_names  = {'headRigid_Rigid_headRigid_quat_W';'headRigid_Rigid_headRigid_quat_Z';...
    %                                                        'headRigid_Rigid_headRigid_quat_X';'headRigid_Rigid_headRigid_quat_Y'};
    % 
    % config.motion.streams{2}.xdfname            = 'YourStreamNameInXDF2';
    % config.motion.streams{2}.bidsname           = tracking_systems{2};
    % config.motion.streams{2}.tracked_points     = {'Rigid1', 'Rigid2', 'Rigid3', 'Rigid4'}; % example when there are multiple points tracked by the system
    % config.motion.streams{2}.positions.channel_names = {'Rigid1_X', 'Rigid2_X', 'Rigid3_X', 'Rigid4_X';... % each column is one tracked point and rows are different coordinates
    %                                                     'Rigid1_Y', 'Rigid2_Y', 'Rigid3_Y', 'Rigid4_Y';...
    %                                                     'Rigid1_Z', 'Rigid2_Z', 'Rigid3_Z', 'Rigid4_Z'};
    % config.motion.streams{2}.quaternions.channel_names = {'Rigid1_A', 'Rigid2_A', 'Rigid3_A', 'Rigid4_A';...
    %                                                     'Rigid1_B', 'Rigid2_B', 'Rigid3_B', 'Rigid4_B';...
    %                                                     'Rigid1_C', 'Rigid2_C', 'Rigid3_C', 'Rigid4_C'; ...
    %                                                     'Rigid1_D', 'Rigid2_D', 'Rigid3_D', 'Rigid4_D'};
    % 

    % bemobil_xdf2bids(config, ...
    %     'general_metadata', generalInfo,...
    %     'participant_metadata', subjectInfo,...
    %     'motion_metadata', motionInfo, ...
    %     'eeg_metadata', eegInfo);

    bemobil_xdf2bids(config, ...
        'eeg_metadata', eegInfo);
end

fclose all;


%% configuration for bemobil bids2set
%----------------------------------------------------------------------
config.set_folder               = fullfile(studyFolder,'2_raw-EEGLAB');
config.session_names            = sessionNames;

% config.other_data_types = {'motion'};  % specify which other data type than eeg is there (only 'motion' and 'physio' supported atm)
bemobil_bids2set(config);


%% Configuration of the BeMoBIL pipeline
bemobil_config = BeMoBIL_Configuration(study_path);

% enter all subjects to process here (you can split it up in more MATLAB instances if you have more CPU power and RAM)
subjects = 18; 

% set to 1 if all files should be computed, independently of whether they are present on disk or not
force_recompute = 0; 


%% processing loop

for subject = subjects
    
    %% prepare filepaths and check if already done
    
	disp(['Subject #' num2str(subject)]);
    
	STUDY = []; CURRENTSTUDY = 0; ALLEEG = [];  CURRENTSET=[]; EEG=[]; EEG_interp_avref = []; EEG_single_subject_final = [];
	
	input_filepath = [bemobil_config.study_folder bemobil_config.raw_EEGLAB_data_folder bemobil_config.filename_prefix num2str(subject)];
	output_filepath = [bemobil_config.study_folder bemobil_config.single_subject_analysis_folder bemobil_config.filename_prefix num2str(subject)];
	
	try
		% load completely processed file
		EEG_single_subject_final = ...
            pop_loadset('filename', [ bemobil_config.filename_prefix num2str(subject)...
			'_' bemobil_config.single_subject_cleaned_ICA_filename], 'filepath', output_filepath);
    catch
        disp('...failed. Computing now.')
    end

	
	if ~force_recompute && exist('EEG_single_subject_final','var') ...
            && ~isempty(EEG_single_subject_final)
		clear EEG_single_subject_final
		disp('Subject is completely preprocessed already.')
		continue
    end

	%% load data in EEGLAB .set structure
    % make sure the data is stored in double precision, large datafiles are
    % supported, no memory mapped objects are used but data is processed 
    % locally, and two files are used for storing sets (.set and .fdt)
	try 
        pop_editoptions('option_saveversion6', 0, 'option_single', 0, ...
            'option_memmapdata', 0, 'option_savetwofiles', 1, ...
            'option_storedisk', 0);
    catch
        warning('Could NOT edit EEGLAB memory options!!'); 
    end
    
    % load files that were created from xdf to BIDS to EEGLAB
    EEG = pop_loadset('filename', ...
        [ bemobil_config.filename_prefix num2str(subject) '_' ...
        bemobil_config.merged_filename],'filepath',input_filepath);
    

    %% reselect EEG channels and remove ACC channels
    EEG = pop_select(EEG, 'rmchannel', [65, 66, 67]);


    %% Load channels location file
    EEG = pop_chanedit(EEG,'load', [processing_path, 'chanlocs.ced']); 


    %% Define and Add events

    % %% Add Events to EEG data
    % [EEG, removeindices] = Main_add_events(EEG, output, ...
    %     subject, subject_id, processing_path, input_filepath, bemobil_config);

    
    % Compute latency values
    
    % start_beep event (first single beep)
    start_beep = find(diff(All_Experiment(6, :)) == 1);
    start_beep = reshape(start_beep, 2, []);
    start_beep_time_Expdata = All_Experiment_time(start_beep(1,:));
    start_beep_indx_EEG = ...
        knnsearch(All_EEG_time', start_beep_time_Expdata');
    start_beep_latency_EEG = EEG.times(start_beep_indx_EEG); % time unit: milisecond
    

    % pressure_change event (2s after first single beep)
    pressure_change_time_Expdata = All_Experiment_time(start_beep(1,:)) + 2;
    pressure_change_indx_EEG = ...
        knnsearch(All_EEG_time', pressure_change_time_Expdata');
    pressure_change_latency_EEG = EEG.times(pressure_change_indx_EEG); % time unit: milisecond


    % start_move event (second single beep, 2s after pressure change)
    start_move = find(diff(All_Experiment(6, :)) == 1);
    start_move = reshape(start_move, 2, []);
    start_move_time_Expdata = All_Experiment_time(start_move(2,:));
    start_move_indx_EEG = ...
        knnsearch(All_EEG_time', start_move_time_Expdata');
    start_move_latency_EEG = EEG.times(start_move_indx_EEG); % time unit: milisecond
    

    % finish_beep event (20s after start_move event, double beep to stop movement)
    finish_beep = find(diff(All_Experiment(6, :)) == -2);
    finish_beep_time_Expdata = All_Experiment_time(finish_beep);
    finish_beep_indx_EEG = ...
        knnsearch(All_EEG_time', finish_beep_time_Expdata');
    finish_beep_latency_EEG = EEG.times(finish_beep_indx_EEG); % time unit: milisecond
    

    % score_press event (experimenter presses the scores immidiately after subjects evaluate the task)
    score_press = find(diff(All_Experiment(7, :)) > 0);
    score_press_time_Expdata = All_Experiment_time(score_press);
    score_press_indx_EEG = ...
        knnsearch(All_EEG_time', score_press_time_Expdata');
    score_press_latency_EEG = EEG.times(score_press_indx_EEG); % time unit: milisecond
    
    % score_press(43) = [];
    % score_press_time_Expdata(43) = [];
    % score_press_indx_EEG(43) = [];
    % score_press_latency_EEG(43) = [];
    
    %% Define and Add Events
    if subject_id == 18
        no_PAM_trials = [1:3, 40:42, 43:45, 76:78, 79:81, 112:114, 115:117, 148:150];
        familiarization_trials = 4:9;
    end
    % Import Trials information
    Trials = cell(1, size(start_beep,2));
    for i = 1:numel(Trials)
        if ismember(i, no_PAM_trials)
            Trials{1, i}.Description = 'No_PAM';
        elseif ismember(i, familiarization_trials)
            Trials{1, i}.Description = 'Familiarizattion';
        else
            Trials{1, i}.Description = 'Experiment';
        end
        Trials{1, i}.Pressure = All_Experiment(3, start_beep(2, i));
        if i ~= numel(Trials)
            Trials{1, i}.Score = All_Experiment(4, start_beep(1, i+1));
        else
            Trials{1, i}.Score = All_Experiment(4, end);
        end
    end
    
    
    % make events description
    desc1 = cell(1, numel(pressure_change_latency_EEG)); % {P(i-1), P(i), Trial}
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
    latency = start_beep_latency_EEG;
    desc = desc1;

    % Add pressure change events
    type = cat(2, type, repmat({'PC_Pressure_Change'}, 1, numel(pressure_change_latency_EEG)));
    latency = cat(2, latency, pressure_change_latency_EEG);
    desc = cat(2, desc, desc1);

    % Add start move events
    type = cat(2, type, repmat({'SM_Start_Move'}, 1, size(start_beep,2)));
    latency = cat(2, latency, start_move_latency_EEG);
    desc = cat(2, desc, desc1);
    
    % Add finish beep events
    type = cat(2, type, repmat({'FB_Finish_Beep'}, 1, numel(finish_beep)));
    latency = cat(2, latency, finish_beep_latency_EEG);
    desc = cat(2, desc, desc1);

    % New description including the scores of current and previous trials
    desc2 = cell(1, numel(score_press_latency_EEG)); % {P(i-1), P(i), Trial}
    for i = 1:length(desc2)
        if i~=1
            desc2{1, i} = {Trials{1, i-1}.Pressure, Trials{1, i}.Pressure, Trials{1, i-1}.Score, Trials{1, i}.Score, i};
        else
            indx_temp = knnsearch(All_Experiment_time', pressure_change_time_Expdata(1) - 1);
            desc2{1, 1} = {All_Experiment(3, indx_temp), Trials{1, i}.Pressure, All_Experiment(4, indx_temp), Trials{1, i}.Score, i};
        end
    end

    % Add score press events
    type = cat(2, type, repmat({'SP_Score_Press'}, 1, numel(score_press_latency_EEG)));
    latency = cat(2, latency, score_press_latency_EEG);
    desc = cat(2, desc, desc2);

    % Add Trial Start events
    type = cat(2, type, repmat({'TS_Trial_Start'}, 1, numel(start_beep_latency_EEG)));
    latency = cat(2, latency, start_beep_latency_EEG - 2); % 2ms before (one sample)
    desc = cat(2, desc, desc1);

    % Add Trial End events
    type = cat(2, type, repmat({'TE_Trial_End'}, 1, numel(score_press_latency_EEG)));
    latency = cat(2, latency, score_press_latency_EEG + 2); % 2ms after (one sample)
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

    folder = [processing_path, 'Events', ...
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
    fprintf(fileID, 'type\tlatency\tdesc\n');
    
    % Write the data
    for i = 1:numel(type)
        % Convert the nested cell array in desc to a string with underline separator
        desc_str = strjoin(cellfun(@num2str, desc{i}, 'UniformOutput', false), '_');
        fprintf(fileID, '%s\t%d\t%s\n', type{i}, latency(i), desc_str);
    end
    
    % Close the file
    fclose(fileID);
    
    % Notify the user
    fprintf('File saved successfully: \n%s\n', filename);
    
    
    %% Add Events to the EEG file 
    [EEG, eventnumbers] = pop_importevent(EEG, 'event', ...
              filename, 'fields', {'type', 'latency','desc' }, ...
              'append', 'no', 'align', NaN, 'skipline', 1, 'timeunit', 1E-3);


    %% individual EEG processing to remove non-exp segments
    % it is stongly recommended to remove these segments because they may contain strong artifacts that confuse channel
    % detection and AMICA
    removeindices = zeros(size(start_beep,2)+1 ,2);
    % remove from start to first event
    removeindices(1, :) = [0 EEG.event(1).latency]; 

    % add more removeIndices here for pauses or itnerruptions of the 
    % experiment if they have markers or you know their indices in the data
    for i = 1:size(start_beep,2)-1
        removeindices(i+1, :) = [EEG.event(7*i).latency EEG.event(7*i+1).latency];
    end
    removeindices(end, :) = [EEG.event(end).latency EEG.pnts]; % remove from last event to the end
    

    %%
    % filter for plot
    EEG_plot = pop_eegfiltnew(EEG, 'locutoff',0.5, 'hicutoff', 40,'plotfreqz',0);
    
    % plot
    fig1 = figure; set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
    plot(normalize(EEG_plot.data') + [1:10:10*EEG_plot.nbchan], 'color', [78 165 216]/255)
    yticks([])
    
    xlim([0 EEG.pnts])
    ylim([-10 10*EEG_plot.nbchan+10])
    
    hold on
    
    % plot lines for valid times
    
    for i = 1:size(removeindices,1)
        plot([removeindices(i,1) removeindices(i,1)],ylim,'k', 'LineStyle','-')
        plot([removeindices(i,2) removeindices(i,2)],ylim,'k', 'LineStyle','--')
    end
    title(['Subject ', num2str(subject), ', Non-Exp Segments: From Solid Line to Next Dashed Line on the Right'])
    
    % save plot
    print(gcf,fullfile(input_filepath,[bemobil_config.filename_prefix num2str(subject) '_raw-full_EEG.png']),'-dpng')
    close

    %% Add other flexion/extension start/end events
    % Important Note: before rejecting the non-exp segments you must add
    % the flexion/Extension start events

    % Use peak selection app to find the flexion and extension start events
    output_peaks = find_peaks_and_select_events(output);

    start_beep = start_beep(1, :);
    start_move = start_move(2, :);
    XLimits = output_peaks.XLimits;

    trial_pks_high_peaks = output_peaks.trial_pks_high_peaks;
    trial_locs_high_peaks = output_peaks.trial_locs_high_peaks;

    trial_pks_low_peaks = output_peaks.trial_pks_low_peaks;
    trial_locs_low_peaks = output_peaks.trial_locs_low_peaks;

    N_Trials = output_peaks.N_Trials;

    Trials_encoder_events = output_peaks.Trials_encoder_events;

    for i = 1:numel(Trials_encoder_events)
        Trials_encoder_events{1, i}.Description = Trials{1, i}.Description;
        Trials_encoder_events{1, i}.Pressure = Trials{1, i}.Pressure;
        Trials_encoder_events{1, i}.Score = Trials{1, i}.Score;
    end


    %% run the app to select unwanted peaks 
    app = find_flexion_extension_events;


    %% load the Trials_encoder_events fully defined
    % Trials_encoder_events = load([study_path, '6_0_Trials_Info_and_Events\sub-', ...
    %     num2str(subject_id), filesep, 'sub-', num2str(subject_id), ...
    %     '_Trials_encoder_events.mat']);


    %% Add flexion and extension indexes to Trials_encoder_events structure (based on Exp stream indexes)
    for i = 1:length(Trials_encoder_events)
        if length(Trials_encoder_events{1, i}.high_peaks.index) > 1
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
    
            
            if flag1 == 1 && flag2 == 1 % case 1 - should not happen!
                
                Trials_encoder_events{1, i}.Case = 1;
    
                Trials_encoder_events{1, i}.low_peaks.index(1) = [];
                Trials_encoder_events{1, i}.low_peaks.time(1) = [];
                Trials_encoder_events{1, i}.low_peaks.value(1) = [];
    
                Trials_encoder_events{1, i}.Flexion_Start = Trials_encoder_events{1, i}.high_peaks.index(1:end-1);
                Trials_encoder_events{1, i}.Flexion_End   = Trials_encoder_events{1, i}.low_peaks.index;
    
                Trials_encoder_events{1, i}.Extension_Start = Trials_encoder_events{1, i}.low_peaks.index;
                Trials_encoder_events{1, i}.Extension_End   = Trials_encoder_events{1, i}.high_peaks.index(2:end);
        
            elseif flag1 == 1 && flag2 == 0 % case 2 - should not happen!
                
                Trials_encoder_events{1, i}.Case = 2;
    
                Trials_encoder_events{1, i}.low_peaks.index(1) = [];
                Trials_encoder_events{1, i}.low_peaks.time(1) = [];
                Trials_encoder_events{1, i}.low_peaks.value(1) = [];
    
                Trials_encoder_events{1, i}.Flexion_Start = Trials_encoder_events{1, i}.high_peaks.index;
                Trials_encoder_events{1, i}.Flexion_End   = Trials_encoder_events{1, i}.low_peaks.index;
    
                Trials_encoder_events{1, i}.Extension_Start = Trials_encoder_events{1, i}.low_peaks.index(1:end-1);
                Trials_encoder_events{1, i}.Extension_End   = Trials_encoder_events{1, i}.high_peaks.index(2:end);
        
            elseif flag1 == 0 && flag2 == 1 % case 3
                
                Trials_encoder_events{1, i}.Case = 3;
    
                Trials_encoder_events{1, i}.Flexion_Start = Trials_encoder_events{1, i}.high_peaks.index(1:end-1);
                Trials_encoder_events{1, i}.Flexion_End   = Trials_encoder_events{1, i}.low_peaks.index;
    
                Trials_encoder_events{1, i}.Extension_Start = Trials_encoder_events{1, i}.low_peaks.index;
                Trials_encoder_events{1, i}.Extension_End   = Trials_encoder_events{1, i}.high_peaks.index(2:end);
    
        
            elseif flag1 == 0 && flag2 == 0 % case 4
                
                Trials_encoder_events{1, i}.Case = 4;
    
                Trials_encoder_events{1, i}.Flexion_Start = Trials_encoder_events{1, i}.high_peaks.index;
                Trials_encoder_events{1, i}.Flexion_End   = Trials_encoder_events{1, i}.low_peaks.index;
    
                Trials_encoder_events{1, i}.Extension_Start = Trials_encoder_events{1, i}.low_peaks.index(1:end-1);
                Trials_encoder_events{1, i}.Extension_End   = Trials_encoder_events{1, i}.high_peaks.index(2:end);
    
            end
        end
    end


    save([study_path, '6_0_Trials_Info_and_Events\sub-', ...
        num2str(subject_id), filesep, 'sub-', num2str(subject_id), ...
        '_Trials_encoder_events.mat'], 'Trials_encoder_events', '-v7.3')


    %% Add events to the text file to import them in EEG data (EEGLAB)
    % Initialize variables
    
    n_trials = length(Trials_encoder_events);
    event_types = {'FlxS', 'FlxE', 'ExtS', 'ExtE'};
    event_fields = {'Flexion_Start', 'Flexion_End', 'Extension_Start', 'Extension_End'};
    
    % Pre-calculate total number of events to preallocate arrays
    total_events = 0;
    for i = 1:n_trials
        for j = 1:numel(event_fields)
            total_events = total_events + numel(Trials_encoder_events{1, i}.(event_fields{j}));
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
            
            event_times = All_Experiment_time(Trials_encoder_events{1, i}.(event_field));
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
    all_event_latency_EEG = EEG.times(all_event_indx_EEG);
    
    % Organize data back into per-trial cell arrays
    type_extended = cell(1, n_trials);
    latency_extended = cell(1, n_trials);
    desc_extended = cell(1, n_trials);
    
    for i = 1:n_trials
        idx_trial = (trial_indices == i);
        type_extended{i} = all_event_labels(idx_trial)';
        latency_extended{i} = all_event_latency_EEG(idx_trial)';
        desc_extended{i} = all_descs(idx_trial)';
    end
    

    %% Appending to the text file
    % Define the file paths
    original_filename = [processing_path, 'Events\sub-', ...
        num2str(subject_id), filesep, 'events_basic.txt'];
    new_filename = [processing_path, 'Events\sub-', ...
        num2str(subject_id), filesep, 'events_with_FlxExt.txt'];
    
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



    % Ensure your data arrays are column vectors
    event_labels = all_event_labels;          % Cell array of strings, size (9730 x 1)
    event_latency = all_event_latency_EEG';   % Numeric array, transpose to (9730 x 1)
    event_descs = all_descs;                  % Cell array of cell arrays, size (9730 x 1)
    
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
        new_line = sprintf('%s\t%f\t%s', event_labels{i}, event_latency(i), event_descs_str{i});
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
    [EEG, ~] = pop_importevent(EEG, 'event', ...
              new_filename, 'fields', {'type', 'latency','desc' }, ...
              'append', 'no', 'align', NaN, 'skipline', 1, 'timeunit', 1E-3);



    %% Reject non-experimental periods
    EEG = eeg_eegrej(EEG, removeindices);

    
    %% processing wrappers for basic processing and AMICA
    
    % do basic preprocessing, line noise removal, and channel interpolation
	[ALLEEG, EEG_preprocessed, CURRENTSET] = bemobil_process_all_EEG_preprocessing(subject, bemobil_config, ALLEEG, EEG, CURRENTSET, force_recompute);

    % start the processing pipeline for AMICA
	bemobil_process_all_AMICA(ALLEEG, EEG_preprocessed, CURRENTSET, subject, bemobil_config, force_recompute);

end

bemobil_copy_plots_in_one(bemobil_config)

subjects
subject

disp('PROCESSING DONE! YOU CAN CLOSE THE WINDOW NOW!')


% rmpath(genpath('C:\Morteza\MyProjects\ANSYMB2024'))
