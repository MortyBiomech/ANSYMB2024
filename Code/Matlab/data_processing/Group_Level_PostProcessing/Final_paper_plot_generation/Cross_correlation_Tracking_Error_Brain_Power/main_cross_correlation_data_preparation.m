clc
clear


%% Add necessary paths
main_project_path = 'D:\Morteza\MyProjects\ANSYMB2024\';

addpath(genpath([main_project_path, 'Code']));
addpath(genpath([main_project_path, 'data\7_STUDY\Epoched_data']));

data_path         = [main_project_path, 'data\'];
Code_path         = [main_project_path, 'Code\Matlab\data_processing\'];
all_STUDY_PATH    = [data_path, '7_STUDY\Epoched_data\', ...
                        'multiple_clustering\'];

icatimef_path     = [data_path, '5_single-subject-EEG-analysis\', ...
                        'timewarp_test\Epoched_data'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
ersp_data_path    = [data_path, '7_STUDY\Epoched_data\Final_figures', ...
                        '\ERSP\Three Pressure Conditions\', ...
                        'p 0.01 ersp results\'];

eeglabPath = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';

Subjects_ICs_in_cluster_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\', ...
    'extracting Subjects and ICs in the brain clusters\']
Detailed_Analysis_on_TF_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\'];

current_path = ['D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab', ...
    '\data_processing\Group_Level_PostProcessing\', ...
    'Final_paper_plot_generation\', ...
    'Cross_correlation_Tracking_Error_Brain_Power'];


%% In the next section, file was created and saved. 
% Please load the Subjects_corresponding_TrialsEpochs_to_icatimef.mat


%% Load data to extract [trial-epoch] pairs on non-EEG streams
% because every stream was processed separately and each time specific
% criteria were considered to remove back epochs, now that we want to do
% the cross-correlation it is not easy to track the removed epochs and this
% makes it not possible to have EEG_power-Tracking_Error associated pair.
% Solution: Use the init_time in the icatimef.trialinfo structure and find
% that init_time in the epoched data structure (Epochs_FlextoFlex_based).
% If the first value in the time vector of the epochs inside 
% Epochs_FlextoFlex_based{1, i}.EEG_stream.Preprocessed
% is equal to init_time then we mark that epoch and later redo the tracking
% error calculation (across time) for the selected epochs which would be
% the same number as icatimef epochs and this way we are sure that the
% cross-correlation is done on the same pairs.

% Run eeglab
cd(eeglabPath)
if ~exist("ALLEEG", "var")
    eeglab
end

% Path and filename of your STUDY
% Choose between:
% Left_Parieto_Occipital; Right_Parieto_Occipital
studyPath = [all_STUDY_PATH, 'Left_Parieto_Occipital'];   
studyFile = 'Left_Parieto_Occipital.study';              

% Load the STUDY (and its dataset list)
[STUDY, ALLEEG] = pop_loadstudy( ...
    'filename', studyFile, ...
    'filepath', studyPath);


%% Load icatimef data structure and Epochs_FlextoFlex_based to extract the
% actual [Trial, epochs] pair which are present in the icatimef structure
% then we can have a correct cross-correlation anaylsis
subject_list = 5:18;
Subjects_Trials_Epochs = cell(length(subject_list), 1);
Subjects_icatimef_trialinfo = cell(length(subject_list), 1);
% check the EXP data and if there is no epochs in a specific trial remove
% the epochs which has the same trial id from icatimef structure
trials_to_remove_from_icatimef = cell(size(Subjects_Trials_Epochs));
for sub = 1:length(subject_list)

    cd(icatimef_path)
    disp(['Loading icatimef structure for subject ', ...
        num2str(subject_list(sub)), ' ...']);
    icatimef = load(['S', num2str(subject_list(sub)), '.icatimef'], ...
        "-mat");
    Subjects_icatimef_trialinfo{sub, 1} = icatimef.trialinfo;
    init_index_icatimef = {icatimef.trialinfo.init_index}.';
    init_index_icatimef = cell2mat(init_index_icatimef);



    %%% extracting the urevent.init_index coloumn from ALLEEG structure
    c = {ALLEEG(sub).urevent.init_index}.';
    init_index_urevent = nan(size(c));
    emptyMask1  = cellfun(@isempty, c);        
    init_index_urevent(~emptyMask1) = [c{~emptyMask1}];  

    %%% extracting the urevent.latency coloumn from ALLEEG structure
    c = {ALLEEG(sub).urevent.latency}.';
    latency_urevent = nan(size(c));
    emptyMask2  = cellfun(@isempty, c);        
    latency_urevent(~emptyMask2) = [c{~emptyMask2}];  

    init_index_urevent(emptyMask1) = [];
    latency_urevent(emptyMask1) = [];



    %%% Load Epochs_FlextoFlex_based
    cd([epoched_data_path, 'sub-', num2str(subject_list(sub))]);
    disp(['Loading Epochs_FlextoFlex_based structure for subject ', ...
        num2str(subject_list(sub)), ' ...']);
    load("Epochs_FlextoFlex_based.mat")
    epochedData = Epochs_FlextoFlex_based;
    clear Epochs_FlextoFlex_based

    time_FlxS_EEG = [];
    trial_vector = [];
    epochs_vector = [];
    for trial = 1:length(epochedData)
        if isempty(epochedData{1, trial}.EEG_stream.Preprocessed.Times)
            trials_to_remove_from_icatimef{sub, 1} = cat(1, ...
                trials_to_remove_from_icatimef{sub, 1}, trial);
            continue
        end
        trial_times = cellfun(@(x) x(1), ...
            epochedData{1, trial}.EEG_stream.Preprocessed.Times);
        time_FlxS_EEG = cat(1, time_FlxS_EEG, trial_times');

        trial_vector = cat(1, trial_vector, ...
            repmat(trial, size(trial_times, 2), 1));
        epochs_vector = cat(2, epochs_vector, 1:length(trial_times));
    end
    epochs_vector = epochs_vector';



    %%% Finding the init_index_icatimef in init_index_urevent
    idxAll = arrayfun(@(x) find(init_index_urevent == x), ...
        init_index_icatimef, 'UniformOutput', false);
    lenPerCell = cellfun(@numel, idxAll);   
    dupMask    = lenPerCell > 1;       
    anyDup = any(dupMask);   
    if anyDup 
        disp(['duplicate init_index values in subject', ...
            num2str(subject_list(sub))]); 
    end
    idxAll = cell2mat(idxAll);

    %%% Corresponding latencies in urevent
    times_icatimef = latency_urevent(idxAll)*2; % multiplication to 2 (ms)


    %%% Selecting the trials and epochs which match the ones from icatimef
    [min_values, idx_min] =  arrayfun(@(x) min(abs(time_FlxS_EEG - x)), ...
        times_icatimef);
    trials_epochs_mins = [trial_vector(idx_min), epochs_vector(idx_min), min_values];


    Subjects_Trials_Epochs{sub} = trials_epochs_mins;

end


% save the corresponding trials-epochs pairs
cd(current_path)
save('Subjects_corresponding_TrialsEpochs_to_icatimef.mat', ...
    "Subjects_Trials_Epochs");
save('trials_to_remove_from_icatimef.mat', ...
    "trials_to_remove_from_icatimef");
save('Subjects_icatimef_trialinfo.mat', ...
    "Subjects_icatimef_trialinfo");





%% Load the tracking errors and scores 
subject_list = 5:18;

inner_struct = struct('time', [], 'tracking_error', [], 'pressure', [], ...
    'score', [], 'trial_epoch', [], 'events', []);
Subject_Tracking_Error = repmat({inner_struct}, length(subject_list), 1); 

for sub = 1:length(subject_list)
    
    % Load the epoched experimental data
    folderName = [epoched_data_path, 'sub-', ...
        num2str(subject_list(sub))];
    cd(folderName)
    disp(['Loading data from subject ', num2str(subject_list(sub)), ' ...'])
    load Epochs_FlextoFlex_based.mat
    load Trials_Info.mat
    % extract the EXP_stream
    EXP_data = cellfun(@(x) x.EXP_stream, Epochs_FlextoFlex_based, ...
        'UniformOutput', false);
    % free up the memory
    clear Epochs_FlextoFlex_based
    cd(current_path)

    % Store the tracking error per epoch 
    for trial = 1:length(EXP_data)
    
        % check if we are at the experimental trial
        % if ~strcmp(Trials_Info{trial}.General.Description, 'Experiment')
        %     continue
        % end

        % non-empty trial
        if isempty(EXP_data{trial}.Times)
            continue
        end


        % Epoch time
        Subject_Tracking_Error{sub}.time = vertcat( ...
            Subject_Tracking_Error{sub}.time, ...
            EXP_data{trial}.Times');


        % Tracking error
        trial_tracking_error = ...
            cellfun(@(x1, x2) abs(x1 - x2), ...
                      EXP_data{trial}.Encoder_angle, ... % x1
                      EXP_data{trial}.Ref_angle, ...     % x2
                      'UniformOutput', false);
        Subject_Tracking_Error{sub}.tracking_error = vertcat( ...
            Subject_Tracking_Error{sub}.tracking_error, ...
            trial_tracking_error');


        % Pressure level
        P = Trials_Info{trial}.General.Pressure;
        Subject_Tracking_Error{sub}.pressure = cat(1, ...
            Subject_Tracking_Error{sub}.pressure, ...
            repmat(P, length(trial_tracking_error), 1));

        
        % Trial score
        score = Trials_Info{trial}.General.Score;
        Subject_Tracking_Error{sub}.score = cat(1, ...
            Subject_Tracking_Error{sub}.score, ...
            repmat(score, length(trial_tracking_error), 1));


        % Trial-Epoch number
        Subject_Tracking_Error{sub}.trial_epoch = cat(1, ...
            Subject_Tracking_Error{sub}.trial_epoch, ...
            [repmat(trial, length(trial_tracking_error), 1), ...
                (1:length(trial_tracking_error)).']);

        
        % Event indexes
        events = [...
        Trials_Info{trial}.Events.EXP_stream.flextoflex_start_indx', ...
        Trials_Info{trial}.Events.EXP_stream.extension_start_indx', ...
        Trials_Info{trial}.Events.EXP_stream.flextoflex_end_indx'] - ...
        repmat(Trials_Info{trial}.Events.EXP_stream.flextoflex_start_indx'-1, ...
               1, 3);
        Subject_Tracking_Error{sub}.events = cat(1, ...
            Subject_Tracking_Error{sub}.events, ...
            events);
        
    end
   
end


%% Time-warping the tracking error data
warpto_subjects = zeros(length(subject_list), 3);
for sub = 1:length(subject_list)


    epochs_error = Subject_Tracking_Error{sub}.tracking_error;
    epochs_time  = Subject_Tracking_Error{sub}.time;
    epochs_event = Subject_Tracking_Error{sub}.events;
    epochs_event = epochs_event(:, 2);
    epochs_event = num2cell(epochs_event, 2);


    epochs_error = cellfun(@(x1, x2) ...
        interp1(x1, x2, linspace(x1(1), x1(end), 3*size(x1, 2)), "linear"), ...
        epochs_time, epochs_error, 'UniformOutput', false);
    Subject_Tracking_Error{sub}.tracking_error = epochs_error;
    
    epochs_event_time = cellfun(@(x1, x2) x1(x2), ...
        epochs_time, epochs_event, 'UniformOutput', false);

    epochs_time = cellfun(@(x1, x2) ...
        interp1(x1, x2, linspace(x1(1), x1(end), 3*size(x1, 2)), "linear"), ...
        epochs_time, epochs_time, 'UniformOutput', false);
    Subject_Tracking_Error{sub}.time = epochs_time;

    [~, epochs_event_2_indx] = cellfun(@(x1, x2) min(abs(x1 - x2)), ...
        epochs_time, epochs_event_time, 'UniformOutput', false);

    epochs_event = cellfun(@(x1, x2) [1, x1, size(x2, 2)], ...
        epochs_event_2_indx, epochs_time, 'UniformOutput', false);
    epochs_event = cell2mat(epochs_event);
    Subject_Tracking_Error{sub}.events = epochs_event;




    % mark outliers and remove their data
    outlier_indx1 = isoutlier(Subject_Tracking_Error{sub}.events(:, 2));
    outlier_indx2 = isoutlier(Subject_Tracking_Error{sub}.events(:, 3));
    outlier_indx  = or(outlier_indx1, outlier_indx2);

    % Subject_Tracking_Error{sub}.time(outlier_indx) = [];
    % Subject_Tracking_Error{sub}.tracking_error(outlier_indx) = [];
    % Subject_Tracking_Error{sub}.pressure(outlier_indx) = [];
    % Subject_Tracking_Error{sub}.trial_epoch(outlier_indx) = [];
    % Subject_Tracking_Error{sub}.score(outlier_indx) = [];
    % Subject_Tracking_Error{sub}.events(outlier_indx, :) = [];

    warpto_subjects(sub, :) = ...
        [1, median(Subject_Tracking_Error{sub}.events(~outlier_indx, 2)), ...
            median(Subject_Tracking_Error{sub}.events(~outlier_indx, 3))];
end


roundNear = 50; % round numbers to the closest multiple of this value
warpingvalues = round(median(warpto_subjects)/roundNear)*roundNear;
warpingvalues = warpingvalues + [1, 0, 0];



%% warping
Flx_L = warpingvalues(2);
Ext_L = warpingvalues(3)-warpingvalues(2);
Subject_Tracking_Error_warped = cell(length(subject_list), 1);
for sub = 1:length(subject_list)
    
    epochs_N = length(Subject_Tracking_Error{sub}.tracking_error);
    epochs_error = Subject_Tracking_Error{sub}.tracking_error;
    epochs_time  = Subject_Tracking_Error{sub}.time;
    epochs_event = Subject_Tracking_Error{sub}.events;
    for i = 1:epochs_N
        
        error = epochs_error{i};
        t     = epochs_time{i};
        event = epochs_event(i, :);
       
        t_new_Flx = linspace(t(1), t(event(2)), Flx_L);
        t_new_Ext = linspace(t(event(2)+1), t(end), Ext_L);
        t_new     = [t_new_Flx, t_new_Ext];

        error_warped = interp1(t, error, t_new, "linear");
        Subject_Tracking_Error_warped{sub} = cat(1, ...
            Subject_Tracking_Error_warped{sub}, error_warped);
        
    end

end



%% Select the epochs which were identified in icatimef structure
Subject_Tracking_Error_warped_final = Subject_Tracking_Error_warped;
Subject_Tracking_Error_final = Subject_Tracking_Error;
Subject_Epochs_toRemove_in_icatimef = cell(length(subject_list), 1);

for sub = 1:length(subject_list)
    
    big = Subject_Tracking_Error{sub, 1}.trial_epoch;
    small = Subjects_Trials_Epochs{sub, 1}(:, [1,2]);

    % remove unwanted trials from small
    % (those which was extracted from icatimef)
    idx_TrialsToRemove = [];
    if ~isempty(trials_to_remove_from_icatimef{sub, 1})
        RT = trials_to_remove_from_icatimef{sub, 1};
        for i = 1:length(RT)
            alltrial_icatimef = ...
                {Subjects_icatimef_trialinfo{sub, 1}.trial}.';
            alltrial_icatimef = ...
                cellfun(@(x) str2num(x), alltrial_icatimef);
            idx_TrialsToRemove = cat(1, idx_TrialsToRemove, ...
                find(alltrial_icatimef == RT(i)));
        end
    end
    small(idx_TrialsToRemove, :) = [];
    

    idxAll = arrayfun(@(i) find(ismember(big, small(i, :), 'rows')), ...
                  (1:size(small,1))', 'UniformOutput', false);
    idxAll = cell2mat(idxAll);
    

    if ~iscell(idxAll)
        keep = nan(size(big, 1), 1);
        keep(idxAll) = 1; 
        Subject_Tracking_Error_warped_final{sub, 1}(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.time(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.tracking_error(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.pressure(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.score(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.trial_epoch(isnan(keep), :) = [];
        Subject_Tracking_Error_final{sub, 1}.events(isnan(keep), :) = [];
    end


    % find not-wanted trials in incatimef.trialinfo structure
    idx_TrialsToRemove = [];
    if ~isempty(trials_to_remove_from_icatimef{sub, 1})
        RT = trials_to_remove_from_icatimef{sub, 1};
        for i = 1:length(RT)
            alltrial_icatimef = ...
                {Subjects_icatimef_trialinfo{sub, 1}.trial}.';
            alltrial_icatimef = ...
                cellfun(@(x) str2num(x), alltrial_icatimef);
            idx_TrialsToRemove = cat(1, idx_TrialsToRemove, ...
                find(alltrial_icatimef == RT(i)));
        end
    end
    Subject_Epochs_toRemove_in_icatimef{sub, 1} = idx_TrialsToRemove;

end

save('Subjects_Epochs_toRemove_in_icatimef.mat', ...
    "Subject_Epochs_toRemove_in_icatimef");
save('Subjects_Tracking_Error_final.mat', ...
    "Subject_Tracking_Error_final");
save('Subjects_Tracking_Error_warped_final.mat', ...
    "Subject_Tracking_Error_warped_final");






