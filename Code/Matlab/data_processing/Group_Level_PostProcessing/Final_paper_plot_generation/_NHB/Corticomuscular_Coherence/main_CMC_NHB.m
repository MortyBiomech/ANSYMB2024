clc
clear


%% Add paths
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\LSL\xdf-Matlab-master'));

EEGLAB_path = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';
data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
rawdata_path = [data_path, '0_source_data\'];
rawEEGLAB_path = [data_path, '2_raw-EEGLAB\']; 
epoched_EEG_path = [data_path, '5_single-subject-EEG-analysis\', ...
    'timewarp_test\Epoched_data'];
cleanedICA_EEG_path = [data_path, '5_single-subject-EEG-analysis\'];


%% load EEGLAB first
this_path = pwd; 
cd(EEGLAB_path)
if ~exist('ALLCOM','var')
	eeglab;
end
cd(this_path)

subject_id = 18;


%% First step: create an EEGLAB-like EMG dataset
% Original continuous EEG dataset (before cuts)
cd([rawEEGLAB_path, 'sub-', num2str(subject_id)])
fileName = ['sub-', num2str(subject_id), '_merged_EEG.set'];
EEG_orig = pop_loadset('filename', fileName);
cd(this_path)
eeg_emptyset = structfun(@(x) [], EEG_orig, 'UniformOutput', false);



% % Cleaned with ICA EEG dataset (after non-exp periods removal; no epochs)
% cd([cleanedICA_EEG_path, 'sub-', num2str(subject_id)]);
% fileName = ['sub-', num2str(subject_id), '_cleaned_with_ICA.set'];
% EEG_cut_old = pop_loadset('filename', fileName);
% cd(this_path)




% Epoched EEG dataset (after non-experimental periods were removed)
cd(epoched_EEG_path)
fileName = ['sub-', num2str(subject_id), '_cleaned_with_ICA_epoched.set'];
EEG_epoched_main = pop_loadset('filename', fileName);
cd(this_path)



% Original EMG data (.xdf files)
% Size: [nEMGch x nSamples]
output = runs_concatenated(subject_id, rawdata_path);
emg_ids = [2, 3, 4, 5]; % should be defined for each subject
EMG_raw = output.All_EMG(emg_ids, :);
emgLabels = {'VL', 'RF', 'GM', 'BF'}; % EMG channel labels
emg_srate = 2000; % EMG sampling rate

[b,a] = butter(4, [20 250]/(emg_srate/2), 'bandpass'); % band-pass filter
EMG_bandpass = filtfilt(b, a, double(EMG_raw') )';   % channels x samples



%% ===  if replicating Roeder & Boonstra et al. 2020  =====================
%  ========================================================================
% Full-wave rectification via Hilbert transform
% The analytic signal: z(t) = x(t) + j*H{x(t)}
% Its magnitude |z(t)| is the instantaneous amplitude (envelope)
[b, a] = butter(4, 20 / (emg_srate/2), 'high');  % High-pass filter
EMG_hpf = filtfilt(b, a, double(EMG_raw') )';
EMG_analytic = hilbert(EMG_hpf);       % [samples x channels], complex
EMG_rect     = abs(EMG_analytic);      % Hilbert amplitude (real, positive)

% Decompose analytic signal: z(t) = A(t) * exp(j*phi(t))
% Normalize amplitude to 1:  z_norm(t) = exp(j*phi(t)) = z(t) / A(t)
% Take real part:             cos(phi(t)) — unit-amplitude, phase only
EMG_demod = real(EMG_analytic ./ abs(EMG_analytic));



% %% test: adding events to EEG_orig
% [EEG_orig_events, removeindices_eeg_orig] = ...
%     import_same_event_on_EEG_stream(EEG_orig, output, ...
%     subject_id, processing_path, data_path);
% 
% EEG_orig_cut = eeg_eegrej(EEG_orig_events, removeindices_eeg_orig);
% 
% EEG_orig_cut = pop_selectevent(EEG_orig_cut, 'type', {'FlxS', 'ExtS', 'ExtE'}, 'deleteevents', 'on');
% EEG_orig_cut = eeg_checkset(EEG_orig_cut, 'eventconsistency');
% 
% for e = 1:length(EEG_orig_cut.event)
% 
%     nums = str2double(strsplit(EEG_orig_cut.event(e).desc,'_')).';
%     if strcmp(EEG_orig_cut.event(e).type, 'boundary')
%         EEG_orig_cut.event(e).trial = 'none';
%         EEG_orig_cut.event(e).cond = 'none';
%     else
%         EEG_orig_cut.event(e).trial = nums(end);
%         EEG_orig_cut.event(e).cond = nums(2);
%     end
% 
% end
% 
% EEG_orig_epoched = pop_epoch(EEG_orig_cut, {'FlxS'}, [-0.5 3.5], 'newname', 'FlxStartEvents epochs', 'epochinfo', 'yes');
% EEG_orig_epoched = eeg_checkset( EEG_orig_epoched );





%% ---------- CREATE EMG EEGLAB-LIKE DATASET ----------
EMG               = eeg_emptyset;
% EMG.data          = double(EMG_bandpass);  % [channels x samples]
EMG.data          = double(EMG_demod);  % [channels x samples]
% EMG.data          = double(EMG_rect);  % [channels x samples]
EMG.nbchan        = size(EMG.data,1); 
EMG.pnts          = size(EMG.data,2);
EMG.trials        = 1;
EMG.srate         = emg_srate;
EMG.xmin          = 0;
EMG.xmax          = (EMG.pnts - 1) / EMG.srate;
EMG.setname       = 'EMG_continuous';
EMG.filename      = '';
EMG.filepath      = '';

for ch = 1:EMG.nbchan
    EMG.chanlocs(ch).labels = emgLabels{ch};
end

EMG = eeg_checkset(EMG);



%% Define and Add events on EMG structure
% I used the same events that was defined based on EEG stream timing
% steps:
% - find the nearest indexes in EEG timing to the EXP stream event timing
% - based on the associated EEG events define the type, latency and desc
%   for EMG events
% - Something I learned in this process was that EEGLAB can store
%   fractional latencies and we do not necessarily need to define the events
%   on the time stamps.
[EMG_events, removeindices_EMG] = ...
    import_same_event_on_EMG_stream(EMG, EEG_orig, output, ...
    subject_id, processing_path, data_path);


% %% Reject non-exp periods from EMG structure
% EMG_nonExpRej = eeg_eegrej(EMG_events, removeindices_EMG);
% 
% 
% %% remove unwanted events
% EMG_nonExpRej = pop_selectevent(EMG_nonExpRej, 'type', {'FlxS', 'ExtS', 'ExtE'}, 'deleteevents', 'on');
% EMG_nonExpRej = eeg_checkset(EMG_nonExpRej, 'eventconsistency');



%% Add trial and cond columns in the event structure
% for e = 1:length(EMG_nonExpRej.event)
% 
%     nums = str2double(strsplit(EMG_nonExpRej.event(e).desc,'_')).';
%     if strcmp(EMG_nonExpRej.event(e).type, 'boundary')
%         EMG_nonExpRej.event(e).trial = 'none';
%         EMG_nonExpRej.event(e).cond = 'none';
%     else
%         EMG_nonExpRej.event(e).trial = nums(end);
%         EMG_nonExpRej.event(e).cond = nums(2);
%     end
% 
% end

for e = 1:length(EMG_events.event)

    nums = str2double(strsplit(EMG_events.event(e).desc,'_')).';
    if strcmp(EMG_events.event(e).type, 'boundary')
        EMG_events.event(e).trial = 'none';
        EMG_events.event(e).cond = 'none';
    else
        EMG_events.event(e).trial = nums(end);
        EMG_events.event(e).cond = nums(2);
    end

end



%% Epoching
EMG_events = pop_selectevent(EMG_events, 'type', {'FlxS', 'ExtS', 'ExtE'}, 'deleteevents', 'on');
EMG_epoched = pop_epoch(EMG_events, {'FlxS'}, [-0.5 3.5], 'newname', 'FlxStartEvents epochs', 'epochinfo', 'yes');

% EMG_epoched = pop_epoch(EMG_nonExpRej, {'FlxS'}, [-0.5 3.5], 'newname', 'FlxStartEvents epochs', 'epochinfo', 'yes');
EMG_epoched = eeg_checkset( EMG_epoched );
EMG_epoched.setname = strcat('EMG sub-', num2str(subject_id), ' Epoched');
EMG_epoched = eeg_checkset(EMG_epoched);

% because of weired issue during pop_saveset that changes the dimension of
% EMG_epoched.data from 3d to 2d (combining all epochs)
EMG_epoched_copy = EMG_epoched; 


%% save EMG_epoched dataset
% filePath = [this_path, filesep, 'EMG_epoched', filesep, ...
%     'sub-', num2str(subject_id)];
% 
% if ~exist(filePath, "dir")
%     mkdir(filePath)
% end
% 
% EMG_epoched.filename = ['sub-', num2str(subject_id), '_EMG_Epoched.set'];
% EMG_epoched.filepath = filePath;
% EMG_epoched = eeg_checkset(EMG_epoched);
% 
% EMG_epoched = pop_saveset(EMG_epoched, ...
%     'filename', EMG_epoched.filename, ...
%     'filepath', EMG_epoched.filepath);


% % after savingn the dataset in 3d data shape, let's bring back the 2d
% % version using the copied backup
% EMG_epoched = EMG_epoched_copy;


% % %% Time warping and save on EMG_epoched_new
% % EMG_epoched_new = EMG_epoched;
% % events = {'FlxS', 'ExtS', 'ExtE'}; %{'RHS', 'LTO', 'LHS', 'RTO', 'RHS'}; %My previously defined events, I specify which to warp to
% % timewarp = make_timewarp(EMG_epoched_new, events, 'baselineLatency', 0, ...
% %     'maxSTDForAbsolute', 3,...
% %     'maxSTDForRelative', 3);
% % timewarp.warpto = median(timewarp.latencies); % Will be used in newtimef, group analysis uses median of these warpto values in mod_std_precomp_v10_1_5_5a.m
% % EMG_epoched_new.timewarp = timewarp;
% % median_latency = median(timewarp.latencies(:,3)); %Warping to the median latency of my 5 events
% % EMG_epoched_new.timewarp.medianlatency = median_latency;
% % 
% % % Getting rid of bad epochs
% % goodepochs = sort([timewarp.epochs]);
% % EMG_epoched_new = eeg_checkset(EMG_epoched_new);
% % notneeded = [];
% % badepochs = setdiff(1:length(EMG_epoched_new.epoch),goodepochs);
% % EMG_epoched_new.etc.badepochs = badepochs;
% % EMG_epoched_new = pop_select(EMG_epoched_new, 'notrial', badepochs );




%% Find matched epochs between EMG dataset and EEG_epoched_main dataset

% find the first FlxS event init_time_EEG on the EMG_epoched dataset
eventtypes_EMG = {EMG_epoched.epoch.eventtype};
eventlatency_EMG = cellfun(@(x) cell2mat(x), ...
    {EMG_epoched.epoch.eventlatency}, 'UniformOutput', false);


eventtype_FlxS_logical = cellfun(@(x) strcmpi(x, 'FlxS'), ...
    eventtypes_EMG, 'UniformOutput', false);
eventlatency_FlxS_logical = cellfun(@(x) x == 0, ...
    eventlatency_EMG, 'UniformOutput', false);
FlxS_indx_EMG = cellfun(@(x, y) find(and(x, y)), ...
    eventtype_FlxS_logical, eventlatency_FlxS_logical, 'UniformOutput', false);
more_than_one_FlxS_indx = cellfun(@(x) numel(x) > 1, FlxS_indx_EMG, ...
    'UniformOutput', false);
more_than_one_FlxS_indx = cell2mat(more_than_one_FlxS_indx);
epochs_with_morethanone_FlxS = find(more_than_one_FlxS_indx);
disp(['sub-', num2str(subject_id), ' EMG_epoched dataset: epochs ', ...
    num2str(epochs_with_morethanone_FlxS), ' should be ingnored'])


init_time_EEG_onEMG_cell = cellfun(@(x) cell2mat(x), ...
    {EMG_epoched.epoch.eventinit_time_EEG}, 'UniformOutput', false);
init_time_EEG_onEMG = cellfun(@(x, y) x(y(1)), ...
    init_time_EEG_onEMG_cell, FlxS_indx_EMG, 'UniformOutput', false);
init_time_EEG_onEMG = cell2mat(init_time_EEG_onEMG);

nNaN_init_time_EEG_onEMG = sum(isnan(init_time_EEG_onEMG));


[uv, ~, ic] = unique(init_time_EEG_onEMG);
counts = accumarray(ic, 1);
repeatedVals = uv(counts > 1);
idxByValue = arrayfun(@(x) find(init_time_EEG_onEMG == x), ...
    repeatedVals, 'UniformOutput', false);

N_EEG_epochs = length(EEG_epoched_main.epoch);
% Note: 
% matched_vector has nEpochsEEG rows and each row shows the corresponding
% EMG epoch number to the EEG epoch number. For example is the 1000th row
% of matched_vector is 1002, this means 1002th epoch from EMG should be
% selected for CMC calculation with 1000th epoch from EEG. 
matched_vector = NaN(N_EEG_epochs, 1);
condition_vector = zeros(N_EEG_epochs, 1);
for i = 1:N_EEG_epochs
    
    % find the first FlxS event
    eventtypes = EEG_epoched_main.epoch(i).eventtype;
    FlxS_type_indx = strcmpi(eventtypes, 'FlxS');
    eventlatencies = EEG_epoched_main.epoch(i).eventlatency;
    FlxS_latency_indx = cell2mat(eventlatencies) == 0;
    FlxS_indx = find(FlxS_type_indx & FlxS_latency_indx);

    % FlxS event init_time from the EEG_epoched_main
    init_time_EEG_main = EEG_epoched_main.epoch(i).eventinit_time(FlxS_indx);
    init_time_EEG_main = init_time_EEG_main{1};

    matched_indx = find(abs(init_time_EEG_onEMG - init_time_EEG_main) < 1e-10);

    if numel(matched_indx) > 1
        disp(['More then one epoch in EMG dataset was found', ...
            ' for the EEG epoch ', num2str(i), ': ', num2str(matched_indx)]);
    elseif any(ismember(matched_indx, epochs_with_morethanone_FlxS))
        matched_vector(i) = NaN;
    elseif isempty(matched_indx)
        matched_vector(i) = NaN;
    else
        matched_vector(i) = matched_indx;
        condition_vector(i) = ...
            str2num(EEG_epoched_main.epoch(i).eventcond{FlxS_indx});
    end

end 
nNaN = sum(isnan(matched_vector));





%% Interpolate the EMG epochs based on the EEG timing
EEG_time = EEG_epoched_main.times;
EMG_time_main = EMG_epoched.times;
EMG_time_new_all = zeros(length(matched_vector), length(EMG_time_main));

EMG_epoched_data_intrp = zeros(4, 2000, length(matched_vector));

for i = 1:length(matched_vector)
    % here "i" is the EEG epoch number
    % and "matched_vector(i)" is the equivalent EMG epoch number
    if isnan(matched_vector(i)) 
        EMG_epoched_data_intrp(:, :, i) = NaN;
        continue
    end

    EEG_epoch_info = EEG_epoched_main.epoch(i);
    EMG_epoch_info = EMG_epoched.epoch(matched_vector(i));
    eventtypes = EMG_epoch_info.eventtype;
    FlxS_type_indx = strcmpi(eventtypes, 'FlxS');
    eventlatencies = EMG_epoch_info.eventlatency;
    FlxS_latency_indx = cell2mat(eventlatencies) == 0;
    FlxS_indx = find(FlxS_type_indx & FlxS_latency_indx);

    if numel(FlxS_indx) > 1 || isempty(FlxS_indx)
        disp(['Check EEG epoch ', num2str(i), ' and EMG epoch ', ...
            num2str(matched_vector(i))]);
        EMG_epoched_data_intrp(:, :, i) = NaN;
        continue
    end

    init_time_lastEMGpoint = ...
        EMG_epoch_info.eventinit_time_lastEMGpoint{FlxS_indx};
    init_time_FlxS_EMG = EMG_epoch_info.eventinit_time{FlxS_indx};
    delta_t = 1000*(init_time_FlxS_EMG - init_time_lastEMGpoint); % mili-second

    % interpolate to EEG timing
    EMG_time = EMG_time_main - delta_t;
    EMG_time_new_all(i, :) = EMG_time;
    EMG_epoched_data = EMG_epoched.data(:, :, matched_vector(i));
    EMG_epoched_data_intrp(:, :, i) = ...
        interp1(EMG_time', EMG_epoched_data', EEG_time', "linear")';


end

nan_trials = find(any(isnan(EMG_epoched_data_intrp(2, :, :)), 2));
epochs_to_keep = ~any(1:size(EMG_epoched_data_intrp, 3) == nan_trials);

condition_vector_final = condition_vector(epochs_to_keep);
EMG_epochs_CMC = EMG_epoched_data_intrp(:, :, epochs_to_keep);

% load the Subjects_ICs_in_clusters and select the Right and Left Primary
% Motor IC-id for the subject and use that for selecting the final epochs
load("Subjects_ICs_in_clusters.mat")
% Left_Prim_Motor
cluster = 'Left_Prim_Motor';
cluster_indx = find(strcmpi(SUBJECTS_ICS(:, 1), cluster));
subject_indx_in_cluster = ...
    find(SUBJECTS_ICS{cluster_indx, 2}.Subjects == subject_id-4);
ic_indx = SUBJECTS_ICS{cluster_indx, 2}.ICs(subject_indx_in_cluster);
LeftPrimMotor_EEG_epochs_CMC = ...
    squeeze(EEG_epoched_main.icaact(ic_indx, :, epochs_to_keep));
% Right_Prim_Motor
cluster = 'Right_Prim_Motor';
cluster_indx = find(strcmpi(SUBJECTS_ICS(:, 1), cluster));
subject_indx_in_cluster = ...
    find(SUBJECTS_ICS{cluster_indx, 2}.Subjects == subject_id-4);
ic_indx = SUBJECTS_ICS{cluster_indx, 2}.ICs(subject_indx_in_cluster);
RightPrimMotor_EEG_epochs_CMC = ...
    squeeze(EEG_epoched_main.icaact(ic_indx, :, epochs_to_keep));
% Left_Parieto_Occipital
cluster = 'Left_Parieto_Occipital';
cluster_indx = find(strcmpi(SUBJECTS_ICS(:, 1), cluster));
subject_indx_in_cluster = ...
    find(SUBJECTS_ICS{cluster_indx, 2}.Subjects == subject_id-4);
ic_indx = SUBJECTS_ICS{cluster_indx, 2}.ICs(subject_indx_in_cluster);
LeftParietoOccipital_EEG_epochs_CMC = ...
    squeeze(EEG_epoched_main.icaact(ic_indx, :, epochs_to_keep));
% Right_Parieto_Occipital
cluster = 'Right_Parieto_Occipital';
cluster_indx = find(strcmpi(SUBJECTS_ICS(:, 1), cluster));
subject_indx_in_cluster = ...
    find(SUBJECTS_ICS{cluster_indx, 2}.Subjects == subject_id-4);
ic_indx = SUBJECTS_ICS{cluster_indx, 2}.ICs(subject_indx_in_cluster);
RightParietoOccipital_EEG_epochs_CMC = ...
    squeeze(EEG_epoched_main.icaact(ic_indx, :, epochs_to_keep));




% %% Set parameters
% params.srate             = EEG_epoched_main.srate;   % 500
% params.cycles            = [3 0.8];
% params.freqs             = [3 60];
% params.nfreqs            = 100;
% params.subject_id_offset = 4;     % subject 18 stored as index 14
% params.motor_clusters    = {'Left_Prim_Motor', 'Right_Prim_Motor', ...
%                             'Left_PreMot_SuppMot', 'Right_PreMot_SuppMot',   ...
%                             'Left_Parieto_Occipital', 'Right_Parieto_Occipital'};
% params.emg_labels        = {'VL', 'RF', 'GM', 'BF'};
% 
% savePath = [this_path, '\CMC_results'];
% 
% warpingvalues = EEG_epoched_main.timewarp.warpto;
% 
% 
% %% Run
% cmc = compute_CMC_subject(18, EEG_epoched_main, EMG_epochs_CMC, ...
%           epochs_to_keep, condition_vector_final, SUBJECTS_ICS, ...
%           warpingvalues, params, savePath);








%% ===  if replicating the Fauvet et al. 2019  ============================
%  ========================================================================
% Up-sampling to the maximum segment of interest length
srate = 500;
% check the timefreq function output with/without up-sampling the epoch
tlimits = [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000;
% find the maximum length Segment of Interest (SOI)
timewarp_latencies = EEG_epoched_main.timewarp.latencies(epochs_to_keep, :);
% find the indeces of events
timewarp_indices = round(eeg_lat2point(reshape(timewarp_latencies, 1, []), ...
    1, 500, tlimits, 1E-3));
timewarp_indices = reshape(timewarp_indices, ...
    size(timewarp_latencies, 1), size(timewarp_latencies, 2));

SOI_lengths = timewarp_latencies(:, 3);
[~, max_SOI_indx] = max(SOI_lengths)
fs_eff_orig = SOI_lengths(max_SOI_indx)./SOI_lengths;

t_start = EEG_epoched_main.xmin;

EMG_epoch_upsampled = cell(length(fs_eff_orig), 1);
EMG_epoch_upsampled_time = cell(length(fs_eff_orig), 1);
timewarp_indices_upsampled = zeros(length(fs_eff_orig), 3);
for i = 1:length(fs_eff_orig)
    [p, q] = rat(fs_eff_orig(i), 1e-8);   % 1e-6 = tolerance
    EMG_epoch_upsampled{i} = resample(EMG_epochs_CMC(:, :, i)', p, q);
    EMG_epoch_upsampled{i} = EMG_epoch_upsampled{i}';

    timewarp_indices_upsampled(i, :) = round(timewarp_indices(i, :) * p/q);
    
    fs_new = 500*p/q;
    % method 1
    EMG_epoch_upsampled_time{i} = ...
        (0:length(EMG_epoch_upsampled{i})-1) / fs_new + t_start;
    
end








%% first plot six epochs including the maximum length SOI
% v = 1:length(fs_eff_orig);
% exclude = max_SOI_indx;
% pool = v(v ~= exclude);
% epochs_to_plot = randsample(pool, 5); 

EMG_time_orig = EEG_epoched_main.times;
epochs_to_plot = [358, 442, 775, 934, 1026]; 
muscle = 3;

figure()
ylimits = [];
for j = 1:5
    subplot(3, 2, j); hold on
    % plot(EMG_epochs_CMC(muscle, :, epochs_to_plot(j)), 'k')
    % % xlim([1, length(EMG_epochs_CMC(muscle, :, epochs_to_plot(j)))])
    % xlim([1, 2500])

    plot(EMG_time_orig, EMG_epochs_CMC(muscle, :, epochs_to_plot(j)), 'k')
    xlim([EMG_time_orig(1) EMG_time_orig(end)])
    
    title(['Epoch ', num2str(epochs_to_plot(j))], 'FontWeight', 'normal')

    n = timewarp_indices(epochs_to_plot(j), :);
    % xline(n([1, 3]), 'b', 'LineWidth', 2);
    % xline(n(2), '--b', 'LineWidth', 1.5);
    xline(EMG_time_orig(n([1, 3])), 'b', 'LineWidth', 2);
    xline(EMG_time_orig(n(2)), '--b', 'LineWidth', 1.5);

    set(gca, 'FontSize', 12)
    
    c = get(gca, 'YLim');
    % set(gca, 'YLim', [-max(abs(c)) max(abs(c))]);
    % set(gca, 'YLim', [0 max(abs(c))]); % in case of using EMG_rectified(t) --> Hilbert transform on EMG data
    set(gca, 'YLim', [-2 2]); % in case of using phase(t) --> Hilbert transform on EMG data
    ylimits = cat(1, ylimits, c);
end
subplot(3, 2, 6); hold on
% plot(EMG_epochs_CMC(muscle, :, max_SOI_indx), 'k')
% % xlim([1, length(EMG_epochs_CMC(muscle, :, max_SOI_indx))])
% xlim([1, 2500])
    
plot(EMG_time_orig, EMG_epochs_CMC(muscle, :, max_SOI_indx), 'k')
xlim([EMG_time_orig(1) EMG_time_orig(end)])

title(['Epoch ', num2str(max_SOI_indx), ' (Max Segment of Interest)'], ...
    'FontWeight', 'normal')
set(gca, 'FontSize', 12)
c = get(gca, 'YLim');
% set(gca, 'YLim', [0 max(abs(c))]); % in case of using EMG_rectified(t) --> Hilbert transform on EMG data
set(gca, 'YLim', [-2 2]); % in case of using phase(t) --> Hilbert transform on EMG data

% xline(timewarp_indices(max_SOI_indx, [1, 3]), 'r', 'LineWidth', 2);
% xline(timewarp_indices(max_SOI_indx, 2), '--r', 'LineWidth', 1.5);
xline(EMG_time_orig(timewarp_indices(max_SOI_indx, [1, 3])), 'r', 'LineWidth', 2);
xline(EMG_time_orig(timewarp_indices(max_SOI_indx, 2)), '--r', 'LineWidth', 1.5);
sgtitle('Before Up-sampling', 'FontSize', 18)

% xlabel('sample', 'FontSize', 14)
xlabel('time (ms)', 'FontSize', 14)



%%
figure()
for j = 1:5
    subplot(3, 2, j); hold on
    % plot(EMG_epoch_upsampled{epochs_to_plot(j)}(muscle, :), 'k')
    % % xlim([1, length(EMG_epoch_upsampled{epochs_to_plot(j)}(muscle, :))])
    % xlim([1, 2500])
    plot(1000*EMG_epoch_upsampled_time{epochs_to_plot(j)}, ...
        EMG_epoch_upsampled{epochs_to_plot(j)}(muscle, :), 'k')
    xlim(1000*[EMG_epoch_upsampled_time{epochs_to_plot(j)}(1) ...
        EMG_epoch_upsampled_time{epochs_to_plot(j)}(end)])
    
    [p, q] = rat(fs_eff_orig(epochs_to_plot(j)), 1e-8);   % 1e-6 = tolerance
    title(['Epoch ', num2str(epochs_to_plot(j)), ...
        ', Fs ~ ', num2str(round(500*p/q))], ...
        'FontWeight', 'normal');

    n = timewarp_indices_upsampled(epochs_to_plot(j), :);
    % xline(n([1, 3]), 'b', 'LineWidth', 2);
    % xline(n(2), '--b', 'LineWidth', 1.5);
    nn = 1000*EMG_epoch_upsampled_time{epochs_to_plot(j)}(n);
    xline(nn([1, 3]), 'b', 'LineWidth', 2);
    xline(nn(2), '--b', 'LineWidth', 1.5);

    set(gca, 'FontSize', 12)
    % set(gca, 'YLim', [-max(abs(ylimits(j, :))) max(abs(ylimits(j, :)))])
    % set(gca, 'YLim', [0 max(abs(ylimits(j, :)))]); % in case of using EMG_rectified(t) --> Hilbert transform on EMG data
    set(gca, 'YLim', [-2 2]); % in case of using phase(t) --> Hilbert transform on EMG data
end
subplot(3, 2, 6); hold on
% plot(EMG_epoch_upsampled{max_SOI_indx}(muscle, :), 'k')
% % xlim([1, length(EMG_epoch_upsampled{max_SOI_indx}(muscle, :))])
% xlim([1, 2500])
plot(1000*EMG_epoch_upsampled_time{max_SOI_indx}, ...
    EMG_epoch_upsampled{max_SOI_indx}(muscle, :), 'k')
xlim(1000*[EMG_epoch_upsampled_time{max_SOI_indx}(1) ...
    EMG_epoch_upsampled_time{max_SOI_indx}(end)])

title(['Epoch: ', num2str(max_SOI_indx), ...
    ' Fs = ', num2str(round(500)), ...
    ' (Max Segment of Interest)'], ...
    'FontWeight', 'normal')
set(gca, 'FontSize', 12)

% xline(timewarp_indices_upsampled(max_SOI_indx, [1, 3]), 'r', 'LineWidth', 2);
% xline(timewarp_indices_upsampled(max_SOI_indx, 2), '--r', 'LineWidth', 1.5);
n = timewarp_indices_upsampled(max_SOI_indx, :);
nn = 1000*EMG_epoch_upsampled_time{max_SOI_indx}(n);
xline(nn([1, 3]), 'r', 'LineWidth', 2);
xline(nn(2), '--r', 'LineWidth', 1.5);


% set(gca, 'YLim', [0 max(abs(ylimits(j, :)))]); % in case of using EMG_rectified(t) --> Hilbert transform on EMG data
set(gca, 'YLim', [-2 2]); % in case of using phase(t) --> Hilbert transform on EMG data

sgtitle('After Up-sampling', 'FontSize', 18)

% xlabel('sample', 'FontSize', 14)
xlabel('time (ms)', 'FontSize', 14)



%% test: timefreq function of the EMG epochs including max_SOI

% on the max SOI epoch
[tmpall_maxSOI, freqs, timesout_maxSOI, ~] = ...
    timefreq( ...
    squeeze(EMG_epochs_CMC(muscle, :, max_SOI_indx))', ...
    srate, ...
    'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
    'freqs', [3, 60], ...
    'nfreqs', 100, 'ntimesout', 400, ...
    'freqscale', 'log', ...
    'cycles', [2, 0.8]);

figure()
subplot(3, 2, 6)
contourf(timesout_maxSOI, freqs, abs(tmpall_maxSOI), 200, 'linecolor', 'none');
set(gca, ...
    'yscale',   'log', ...
    'ydir',     'norm', ...
    'ylim',     [3, 60], ...
    'ytick',    [4 8 14 30 60], ...
    'box',      'on', ...
    'FontSize', 14);
colormap("hot")
colorbar


event_times_TF = timesout_maxSOI(knnsearch(timesout_maxSOI', timewarp_latencies(max_SOI_indx, :)'));

xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')

ylabel('Frequency (Hz)')
xlabel('Time (ms)')
title(['Epoch ', num2str(max_SOI_indx), ' - max Segment of Interest TF plot'], ...
    'FontWeight', 'bold')


tmpall_e = cell(5, 1);
timesout_e = cell(5, 1);
for i = 1:5

    % on the randomely selected epochs before up-sampling
    [tmpall_epoch, freqs, timesout_epoch, ~] = ...
        timefreq( ...
        squeeze(EMG_epochs_CMC(muscle, :, epochs_to_plot(i)))', ...
        srate, ...
        'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
        'freqs', [3, 60], ...
        'nfreqs', 100, 'ntimesout', 400, ...
        'freqscale', 'log', ...
        'cycles', [2, 0.8]);

    tmpall_e{i} = tmpall_epoch;
    timesout_e{i} = timesout_epoch;

    subplot(3, 2, i)
    contourf(timesout_epoch, freqs, abs(tmpall_epoch), 200, 'linecolor', 'none');
    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     [3, 60], ...
        'ytick',    [4 8 14 30 60], ...
        'box',      'on', ...
        'FontSize', 14);
    colormap("hot")
    colorbar
    
    event_times_TF = timesout_epoch(knnsearch(timesout_epoch', ...
        timewarp_latencies(epochs_to_plot(i), :)'));
    
    xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
    xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')
    
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    title(['Epoch ', num2str(epochs_to_plot(i))], ...
        'FontWeight', 'normal')

end

sgtitle('Before Up-sampling', 'FontSize', 18)




% %% test: timefreq function on the up-sampled EMG epochs with 500 Hz srate
% 
% % on the max SOI epoch
% [tmpall_maxSOI, freqs, timesout_maxSOI, ~] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(muscle, :, max_SOI_indx))', ...
%     srate, ...
%     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, 'ntimesout', 400, ...
%     'freqscale', 'log', ...
%     'cycles', [2, 0.8]);
% 
% figure()
% subplot(3, 2, 6)
% contourf(timesout_maxSOI, freqs, abs(tmpall_maxSOI), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on', ...
%     'FontSize', 14);
% colormap("hot")
% colorbar
% 
% 
% event_times_TF = timesout_maxSOI(knnsearch(timesout_maxSOI', timewarp_latencies(max_SOI_indx, :)'));
% 
% xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
% xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')
% 
% ylabel('Frequency (Hz)')
% xlabel('Time (ms)')
% title(['Epoch ', num2str(max_SOI_indx), ' - max Segment of Interest TF plot'], ...
%     'FontWeight', 'bold')
% 
% 
% tmpall_e_upsample = cell(5, 1);
% timesout_e_upsample = cell(5, 1);
% % [EMG_epoch_upsampled_time{epochs_to_plot(i)}(1) ...
% %  EMG_epoch_upsampled_time{epochs_to_plot(i)}(end)]*1000
% for i = 1:5
% 
%     % on the randomely selected epochs before up-sampling
%     [tmpall_epoch, freqs, timesout_epoch, ~] = ...
%         timefreq( ...
%         squeeze(EMG_epoch_upsampled{epochs_to_plot(i)}(muscle, :))', ...
%         srate, ...
%         'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
%         'freqs', [3, 60], ...
%         'nfreqs', 100, 'ntimesout', 400, ...
%         'freqscale', 'log', ...
%         'cycles', [2, 0.8]);
% 
%     tmpall_e_upsample{i} = tmpall_epoch;
%     timesout_e_upsample{i} = timesout_epoch;
% 
%     subplot(3, 2, i)
%     contourf(timesout_epoch, freqs, abs(tmpall_epoch), 200, 'linecolor', 'none');
%     set(gca, ...
%         'yscale',   'log', ...
%         'ydir',     'norm', ...
%         'ylim',     [3, 60], ...
%         'ytick',    [4 8 14 30 60], ...
%         'box',      'on', ...
%         'FontSize', 14);
%     colormap("hot")
%     colorbar
% 
%     event_times_TF = timesout_epoch(knnsearch(timesout_epoch', ...
%         timewarp_latencies(epochs_to_plot(i), :)'));
% 
%     xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
%     xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')
% 
%     ylabel('Frequency (Hz)')
%     xlabel('Time (ms)')
%     title(['Epoch ', num2str(epochs_to_plot(i))], ...
%         'FontWeight', 'normal')
% 
% end
% 
% sgtitle('After Up-sampling with same srate in timefreq function (500 Hz)', 'FontSize', 18)




%% test: timefreq function on the up-sampled EMG epochs with fs_effective

% on the max SOI epoch
[tmpall_maxSOI, freqs, timesout_maxSOI, ~] = ...
    timefreq( ...
    squeeze(EMG_epochs_CMC(muscle, :, max_SOI_indx))', ...
    srate, ...
    'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
    'freqs', [3, 60], ...
    'nfreqs', 100, 'ntimesout', 400, ...
    'freqscale', 'log', ...
    'cycles', [2, 0.8]);

figure()
subplot(3, 2, 6)
contourf(timesout_maxSOI, freqs, abs(tmpall_maxSOI), 200, 'linecolor', 'none');
set(gca, ...
    'yscale',   'log', ...
    'ydir',     'norm', ...
    'ylim',     [3, 60], ...
    'ytick',    [4 8 14 30 60], ...
    'box',      'on', ...
    'FontSize', 14);
colormap("hot")
colorbar


event_times_TF = timesout_maxSOI(knnsearch(timesout_maxSOI', timewarp_latencies(max_SOI_indx, :)'));

xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')

ylabel('Frequency (Hz)')
xlabel('Time (ms)')
title(['Epoch ', num2str(max_SOI_indx), ' - max Segment of Interest TF plot'], ...
    'FontWeight', 'bold')


tmpall_e_upsample_eff = cell(5, 1);
timesout_e_upsample_eff = cell(5, 1);
% [EMG_epoch_upsampled_time{epochs_to_plot(i)}(1) ...
%  EMG_epoch_upsampled_time{epochs_to_plot(i)}(end)]*1000
for i = 1:5

    [p, q] = rat(fs_eff_orig(epochs_to_plot(i)), 1e-8);   % 1e-6 = tolerance
    fs_new = 500*p/q;

    % on the randomely selected epochs before up-sampling
    [tmpall_epoch, freqs, timesout_epoch, ~] = ...
        timefreq( ...
        squeeze(EMG_epoch_upsampled{epochs_to_plot(i)}(muscle, :))', ...
        fs_new, ...
        'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
        'freqs', [3, 60], ...
        'nfreqs', 100, 'ntimesout', 400, ...
        'freqscale', 'log', ...
        'cycles', [2, 0.8]);

    tmpall_e_upsample_eff{i} = tmpall_epoch;
    timesout_e_upsample_eff{i} = timesout_epoch;

    subplot(3, 2, i)
    contourf(timesout_epoch, freqs, abs(tmpall_epoch), 200, 'linecolor', 'none');
    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     [3, 60], ...
        'ytick',    [4 8 14 30 60], ...
        'box',      'on', ...
        'FontSize', 14);
    colormap("hot")
    colorbar
    
    event_times_TF = timesout_epoch(knnsearch(timesout_epoch', ...
        timewarp_latencies(epochs_to_plot(i), :)'));
    
    xline(event_times_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
    xline(event_times_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')
    
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    title(['Epoch ', num2str(epochs_to_plot(i)), ...
        ' - srate effective: ', num2str(round(fs_new))], ...
        'FontWeight', 'normal')

end

sgtitle('After Up-sampling with effective sampling rate in timefreq function', 'FontSize', 18)





%% test: timefreq function with cycle (%) instead of time

event1_maxSOI_indx = timewarp_indices(max_SOI_indx, 1);
event3_maxSOI_indx = timewarp_indices(max_SOI_indx, 3);
cycle_percent_maxSOI = linspace(0, 100, ...
    event3_maxSOI_indx - event1_maxSOI_indx + 1);

delta_percent = cycle_percent_maxSOI(2) - cycle_percent_maxSOI(1);
before_event1 = -delta_percent*(event1_maxSOI_indx-1):delta_percent:-delta_percent;
after_event3  = 100+delta_percent : delta_percent : 100+(2000-event3_maxSOI_indx)*delta_percent;

cycle_full_percent_maxSOI = [before_event1, cycle_percent_maxSOI, after_event3];

% plot the max_SOI TF considering the percent min max as tlimit
[tmpall_maxSOI_percent, freqs, percentsout_maxSOI, ~] = ...
    timefreq( ...
    squeeze(EMG_epochs_CMC(muscle, :, max_SOI_indx))', ...
    srate, ...
    'tlimits', [cycle_full_percent_maxSOI(1) cycle_full_percent_maxSOI(end)], ...
    'freqs', [3, 60], ...
    'nfreqs', 100, 'ntimesout', 400, ...
    'freqscale', 'log', ...
    'cycles', [2, 0.8]);

figure()
subplot(3, 2, 6)
contourf(percentsout_maxSOI, freqs, abs(tmpall_maxSOI_percent), 200, 'linecolor', 'none');
set(gca, ...
    'yscale',   'log', ...
    'ydir',     'norm', ...
    'ylim',     [3, 60], ...
    'ytick',    [4 8 14 30 60], ...
    'box',      'on', ...
    'FontSize', 14);
colormap("hot")
colorbar


percent_events = timewarp_indices(max_SOI_indx, :)*delta_percent + ...
    cycle_full_percent_maxSOI(1) - delta_percent;
event_percents_TF = percentsout_maxSOI(knnsearch(percentsout_maxSOI', percent_events'));

xline(event_percents_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
xline(event_percents_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')

ylabel('Frequency (Hz)')
xlabel('Cycle (%)')
title(['Epoch ', num2str(max_SOI_indx), ' - max Segment of Interest TF plot'], ...
    'FontWeight', 'bold')



tmpall_e_upsample_eff_percent = cell(5, 1);
timesout_e_upsample_eff_percent = cell(5, 1);
% [EMG_epoch_upsampled_time{epochs_to_plot(i)}(1) ...
%  EMG_epoch_upsampled_time{epochs_to_plot(i)}(end)]*1000
for i = 1:5

    % create the cycle_percent for the epoch
    event1_indx = timewarp_indices_upsampled(epochs_to_plot(i), 1);
    event3_indx = timewarp_indices_upsampled(epochs_to_plot(i), 3);
    cycle_percent_epoch = linspace(0, 100, ...
        event3_indx - event1_indx + 1);
    
    delta_percent = cycle_percent_epoch(2) - cycle_percent_epoch(1);
    before_event1 = -delta_percent*(event1_indx-1):delta_percent:-delta_percent;
    L = size(EMG_epoch_upsampled{epochs_to_plot(i)}, 2);
    after_event3  = 100+delta_percent:delta_percent:100+(L-event3_indx)*delta_percent;
    
    cycle_full_percent_epoch = [before_event1, cycle_percent_epoch, after_event3];


    [p, q] = rat(fs_eff_orig(epochs_to_plot(i)), 1e-8);   % 1e-6 = tolerance
    fs_new = 500*p/q;

    % on the randomely selected epochs before up-sampling
    [tmpall_epoch, freqs, percentsout_epoch, ~] = ...
        timefreq( ...
        squeeze(EMG_epoch_upsampled{epochs_to_plot(i)}(muscle, :))', ...
        fs_new, ...
        'tlimits', [cycle_full_percent_epoch(1) cycle_full_percent_epoch(end)], ...
        'timesout', percentsout_maxSOI, ...
        'freqs', [3, 60], ...
        'nfreqs', 100, 'ntimesout', 400, ...
        'freqscale', 'log', ...
        'cycles', [2, 0.8]);

    tmpall_e_upsample_eff_percent{i} = tmpall_epoch;
    timesout_e_upsample_eff_percent{i} = percentsout_epoch;

    subplot(3, 2, i)
    contourf(percentsout_epoch, freqs, abs(tmpall_epoch), 200, 'linecolor', 'none');
    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     [3, 60], ...
        'ytick',    [4 8 14 30 60], ...
        'box',      'on', ...
        'FontSize', 14);
    colormap("hot")
    colorbar
    
    percent_events = timewarp_indices_upsampled(epochs_to_plot(i), :)*delta_percent + ...
        cycle_full_percent_epoch(1) - delta_percent;
    event_percents_TF = percentsout_maxSOI(knnsearch(percentsout_maxSOI', percent_events'));
   
    xline(event_percents_TF([1, 3]), 'LineWidth', 2, 'Color', 'g')
    xline(event_percents_TF(2), 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--')
    
    ylabel('Frequency (Hz)')
    xlabel('Cycle (%)')
    title(['Epoch ', num2str(epochs_to_plot(i)), ...
        ' - srate effective: ', num2str(round(fs_new))], ...
        'FontWeight', 'normal')

end

sgtitle('After Up-sampling with effective sampling rate in timefreq function', 'FontSize', 18)



%% Claude function to test the percent of cycle technique 27.04.2026

params.srate     = 500;
params.freqs     = [3 60]; % [3 60]
params.nfreqs    = 60; % 100
params.cycles    = [2 0.8];
params.freqscale = 'log';
params.ntimesout = 200; % 400

% timewarp_indices: nEpochs x 3, FlxS / ExtS / ExtE indices
[TF_X, TF_Y, percent_grid, freqs, fs_eff, mid_event_pct] = ...
    compute_TF_timenorm( ...
        LeftPrimMotor_EEG_epochs_CMC, ...     % [2000 x nEpochs]
        EMG_epochs_CMC, ... % [4 x 2000 x nEpochs]
        timewarp_indices, ...                 % [nEpochs x 3]
        params);



opts.title = sprintf('CMC sub-%d (Left M1 IC)', subject_id);
[coh, fh] = compute_and_plot_CMC(TF_X, TF_Y, percent_grid, freqs, opts);



opts.title = sprintf('CMC sub-%d (Left M1 IC)', subject_id);
[coh, sig_mask, lambda, fig_handle] = compute_and_plot_CMC_with_shuffle( ...
    TF_X, TF_Y, percent_grid, freqs, opts);


opts.title = sprintf('CMC sub-%d (Left M1 IC)', subject_id);
[coh, sig_mask, lambda, cluster_info, fig_handle] = ...
    compute_and_plot_CMC_cluster_perm(TF_X, TF_Y, percent_grid, freqs, opts);


% vastus_TF = TF_Y(:, :, :, 1);
% recfem_TF = TF_Y(:, :, :, 2);
% gastroc_TF = TF_Y(:, :, :, 3);
% bicepfem_TF = TF_Y(:, :, :, 4);
% 
% 
% a = abs(sum( TF_X.*conj(vastus_TF) , 3, "omitmissing"));
% b = sqrt( sum( abs(TF_X) , 3, "omitmissing") .* ...
%           sum( abs(vastus_TF) , 3, "omitmissing"));
% cmc_leftPrim_vastus = a ./ b;

figure()
contourf(percent_grid, freqs, abs(TF_X(:, :, 358)), 200, 'linecolor', 'none');
set(gca, ...
    'yscale',   'log', ...
    'ydir',     'norm', ...
    'ylim',     [3, 60], ...
    'ytick',    [4 8 14 30 60], ...
    'box',      'on');
colormap('turbo')




% TF_EEG: [100 x 400 x nEpochs]
% TF_EMG: [100 x 400 x nEpochs x 4]
% percent_grid: 1 x 400, in % of cycle








%%
% %%
% subplot(2, 1, 2)
% plot(EMG_epoch_i_upsampled)
% xlim([1, length(EMG_epoch_i_upsampled)])
% 
% 
% 
% %% test
% 
% [tmpall, freqs, timesout, itcvals] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
%     500, ...
%     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, ...
%     'freqscale', 'log', ...
%     'cycles', [3, 0.8]);
% 
% [tmpall2, freqs, timesout2, itcvals2] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
%     500, ...
%     'tlimits', [-1 EEG_epoched_main.xmax]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, ...
%     'freqscale', 'log', ...
%     'cycles', [3, 0.8]);
% 
% [tmpall3, freqs, timesout3, itcvals3] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
%     500, ...
%     'tlimits', [-1 EEG_epoched_main.xmax]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, ...
%     'freqscale', 'log', ...
%     'cycles', [3, 0.8], ...
%     'ntimesout', 300);
% 
% [tmpall4, freqs, timesout4, itcvals4] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
%     500, ...
%     'tlimits', [-0.2 2.5]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, ...
%     'freqscale', 'log', ...
%     'cycles', [3, 0.8], ...
%     'ntimesout', 300);
% 
% [tmpall5, freqs, timesout5, itcvals5] = ...
%     timefreq( ...
%     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
%     500, ...
%     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
%     'freqs', [3, 60], ...
%     'nfreqs', 100, ...
%     'freqscale', 'log', ...
%     'cycles', [3, 0.8], ...
%     'ntimesout', 1000);
% 
% figure()
% contourf(timesout, freqs, abs(tmpall(:, :, 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% title('tlimits: [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ntimes: 200')
% colormap("turbo")
% colorbar
% 
% figure()
% contourf(timesout2, freqs, abs(tmpall2(:, :, 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% title('tlimits: [-1 EEG_epoched_main.xmax]*1000, ntimes: 200')
% colormap("turbo")
% colorbar
% 
% figure()
% contourf(timesout3, freqs, abs(tmpall3(:, :, 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% title('tlimits: [-1 EEG_epoched_main.xmax]*1000, ntimes: 300')
% colormap("turbo")
% colorbar
% 
% figure()
% contourf(timesout4, freqs, abs(tmpall4(:, :, 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% title('tlimits: [-0.2 2.5]*1000, ntimes: 300')
% colormap("turbo")
% colorbar
% 
% figure()
% contourf(timesout5, freqs, abs(tmpall5(:, :, 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% title('tlimits: [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ntimes: 1000')
% colormap("turbo")
% colorbar
% 
% 
% %% test
% % how EMG rectification affects the results
% 
% tlimits = [EEG_epoched_main.xmin, EEG_epoched_main.xmax]*1000;
% 
% x = reshape(LeftPrimMotor_EEG_epochs_CMC, 1, []);
% 
% muscle = 1;
% y = reshape(EMG_epochs_CMC(muscle, :, :), 1, []);
% y_rect = abs(y);
% 
% [coherres, mbase, timesout, freqs, ~, Rangle, ...
%  crossspec_warped, crossspec_noWarp, alltfX, alltfX_pow_warped, ...
%  alltfY, alltfY_pow_warped] =                 ...
%         my_newcrossf(x, y, 2000, tlimits, 500, [3, 0.8], ...
%         'freqs', [3 60], 'nfreqs', 100,                   ...
%         'timewarp',   EEG_epoched_main.timewarp.latencies(epochs_to_keep,:), ...
%         'timewarpms', EEG_epoched_main.timewarp.warpto,     ...
%         'plotamp', 'off', 'plotphase', 'off', 'type', 'coher');
% 
% 
% figure()
% coh_num = sum(crossspec_warped, 3);
% coh_denum = sqrt(sum(alltfX_pow_warped, 3) .* sum(alltfY_pow_warped, 3));
% coh = abs(coh_num ./ coh_denum);
% contourf(timesout, freqs, coh, 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - No EMG Rectification (manual calc)', 'FontSize', 18)
% 
% 
% 
% [coherres_rect, mbase, timesout, freqs, ~, Rangle, ...
%  crossspec_warped_rect, crossspec_noWarp_rect, alltfX, alltfX_pow_warped, ...
%  alltfY_rect, alltfY_pow_warped_rect] =                 ...
%         my_newcrossf(x, y_rect, 2000, tlimits, 500, [3, 0.8], ...
%         'freqs', [3 60], 'nfreqs', 100,                   ...
%         'timewarp',   EEG_epoched_main.timewarp.latencies(epochs_to_keep,:), ...
%         'timewarpms', EEG_epoched_main.timewarp.warpto,     ...
%         'plotamp', 'off', 'plotphase', 'off', 'type', 'coher');
% 
% 
% 
% figure()
% coh_num = sum(crossspec_warped_rect, 3);
% coh_denum = sqrt(sum(alltfX_pow_warped, 3) .* sum(alltfY_pow_warped_rect, 3));
% coh = abs(coh_num ./ coh_denum);
% contourf(timesout, freqs, coherres_rect, 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - With EMG Rectification (manual calc)', 'FontSize', 18)
% 
% 
% 
% 
% %% test
% % check the alltfY_pow_warped and abs(alltfY)
% 
% figure()
% contourf(timesout, freqs, sqrt(sum(alltfX_pow_warped(:, :, :), 3)), 200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Power', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - EEG Power', 'FontSize', 18)
% 
% 
% % figure()
% % contourf(timesout, freqs, sqrt(sum(abs(alltfX(:, :, :)).^2, 3)), 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 60], ...
% %     'ytick',    [4 8 14 30 60], ...
% %     'box',      'on');
% % 
% % colormap(hot);
% % c = colorbar;
% % ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% % xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % title('No Timewarp - EEG Power', 'FontSize', 18)
% 
% 
% 
% 
% %% test
% % plot the denumenator of coherres of my_newcrossf
% 
% figure()
% % norm_denum = sqrt(mean(alltfX_pow_warped./max(max(alltfX_pow_warped, [], 1), [], 2), 3) .* ...
% %          mean(alltfY_pow_warped./max(max(alltfY_pow_warped, [], 1), [], 2), 3));
% denum = sqrt(sum(alltfX_pow_warped, 3) .* ...
%          sum(alltfY_pow_warped_rect, 3));
% contourf(timesout, freqs, ...
%     denum, ...
%     200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Power', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - Denumenator [sqrt(EEG*EMGRect)]', 'FontSize', 18)
% 
% 
% 
% 
% %% test
% % plot the cross spectrum warped (numerator)
% 
% figure()
% contourf(timesout, freqs, ...
%     abs( sum(crossspec_warped, 3) ), ...
%     200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Power', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - Numerator [crossspec]', 'FontSize', 18)
% 
% 
% 
% %% test
% % plot the coherence 
% 
% figure()
% contourf(timesout, freqs, ...
%     abs(sum(crossspec_warped, 3)) ./ (sqrt(sum(alltfX_pow_warped, 3) .* ...
%                                            sum(alltfY_pow_warped, 3))), ...
%     200, 'linecolor', 'none');
% set(gca, ...
%     'yscale',   'log', ...
%     'ydir',     'norm', ...
%     'xlim',     [0, 3000], ...
%     'ylim',     [3, 60], ...
%     'ytick',    [4 8 14 30 60], ...
%     'box',      'on');
% 
% colormap(hot);
% c = colorbar;
% ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% title('Timewarped - Coherence (manual calc)', 'FontSize', 18)
% 
% 
% 
% 
% 
% 
% 
% %% test timefreq of interpolated EMG vs. non-imterpolated EMG
% % [tmpall, freqs, timesout, itcvals] = ...
% %     timefreq( ...
% %     squeeze(EMG_epoched.data(1, :, matched_vector(1:10))), ...
% %     2000, ...
% %     'tlimits', [EMG_epoched.xmin EMG_epoched.xmax]*1000, ...
% %     'freqs', [3, 200], ...
% %     'nfreqs', 150, ...
% %     'freqscale', 'log', ...
% %     'cycles', [3, 0.8]);
% % 
% % figure()
% % contourf(timesout, freqs, abs(tmpall(:, :, 3)), 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 200], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% % 
% % 
% % [tmpall_intrp, freqs_intrp, timesout_intrp, itcvals] = ...
% %     timefreq( ...
% %     squeeze(EMG_epochs_CMC(1, :, 1:10)), ...
% %     500, ...
% %     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
% %     'freqs', [3, 200], ...
% %     'nfreqs', 150, ...
% %     'freqscale', 'log', ...
% %     'cycles', [3, 0.8]);
% % 
% % figure()
% % contourf(timesout_intrp, freqs, abs(tmpall_intrp(:, :, 3)) - abs(tmpall(:, :, 3)), 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 200], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% 
% % g = finputcheck(varargin, ...
% %     { 'ntimesout'     'integer'  []                     []; ...
% %     'timesout'      'real'     []                       []; ...
% %     'winsize'       'integer'  [0 Inf]                  []; ...
% %     'tlimits'       'real'     []                       []; ...
% %     'detrend'       'string'   {'on','off'}             'off'; ...
% %     'causal'        'string'   {'on','off'}             'off'; ...
% %     'verbose'       'string'   {'on','off'}             'on'; ...
% %     'freqs'         'real'     [0 Inf]                  []; ...
% %     'nfreqs'        'integer'  [0 Inf]                  []; ...
% %     'freqscale'     'string'   { 'linear','log','' }    'linear'; ...
% %     'ffttaper'      'string'   { 'hanning','hamming','blackmanharris','none' }  'hanning';
% %     'wavelet'       'real'     [0 Inf]                  0; ...
% %     'cycles'        {'real','integer'}    [0 Inf]       0; ...
% %     'padratio'      'integer'  [1 Inf]                  2; ...
% %     'itctype'       'string'   {'phasecoher','phasecoher2','coher'}  'phasecoher'; ...
% %     'subitc'        'string'   {'on','off'}             'off'; ...
% %     'timestretch'   'cell'     []                       {}; ...
% %     'wletmethod'    'string'   {'dftfilt2','dftfilt3'}    'dftfilt3'; ...
% %     });
% 
% 
% 
% 
% 
% 
% %% Time-resolved CMC with NEWCROSSF function
% % first try the function
% % x = reshape(LeftPrimMotor_EEG_epochs_CMC, 1, []);
% % figure()
% % plot(LeftPrimMotor_EEG_epochs_CMC(:, 2), 'LineWidth', 2, 'Color', 'b')
% % hold on
% % plot(x(2001:4000), 'r--', 'LineWidth', 2)
% % muscle = 1;
% % y = reshape(EMG_epochs_CMC(muscle, :, :), 1, []);
% % figure()
% % plot(EMG_epochs_CMC(muscle, :, 2), 'LineWidth', 2, 'Color', 'b')
% % hold on
% % plot(y(2001:4000), 'r--', 'LineWidth', 2)
% % frames = 2000;
% % tlimits = [EEG_epoched_main.xmin, EEG_epoched_main.xmax]*1000;
% % srate = 500;
% % cycles = [3, 0.8];
% 
% % [coh, mcoh, timesout, freqsout, ~, cohangles, allcoher, alltfX, alltfY] = ...
% %     newcrossf(x, y, 2000, [-500 3500], 500, [3 0.8], ...
% %     'freqscale', 'log', 'nfreqs', 250, 'freqs', [3 100], 'type', 'coher', ...
% %     'plotamp', 'off', 'plotphase', 'off');
% 
% 
% % timewarp alltfX and alltfY then compute coherence
% % [coh, mcoh, timesout, freqsout, ~, cohangles, allcoher, alltfX, alltfY] = ...
% %     my_newcrossf_v1(x, y, frames, tlimits, srate, cycles, ...
% %     'freqscale', 'log', 'nfreqs', 150, 'freqs', [3 200], ...
% %     'padratio', 2, 'type', 'coher', ...
% %     'timewarp',   EEG_epoched_main.timewarp.latencies(epochs_to_keep, :), ...
% %     'timewarpms', EEG_epoched_main.timewarp.warpto, ...
% %     'plotamp', 'off', 'plotphase', 'off');
% 
% 
% % [STUDY, ALLEEG] = mod_std_precomp_v_forEEGlabv2021(STUDY, ALLEEG, 'components', ...
% %     'ersp','on','itc','off', ...
% %     'erspparams', {'cycles', [3 0.8], 'freqs', [3 130], 'nfreqs', 250, ...
% %                    'padratio', 2, 'alpha', NaN, 'freqscale', 'log', ...
% %                    'savetrials', 'off', 'baseline', 'median latency baseline', ...
% %                    'basenorm', 'off', 'trialbase', 'off', ...
% %                    'timewarp', [0], 'timewarpms', warpingvalues}, 'recompute','on');
% 
% 
% % frames = 2000;
% % tlimits = [EEG_epoched_main.xmin, EEG_epoched_main.xmax]*1000;
% % srate = 500;
% % cycles = [3, 0.8];
% % 
% % [coherres, mbase, timesout, freqs, ~, Rangle, ...
% %  crossspec_warped, crossspec_noWarp, alltfX, alltfX_pow_warped, ...
% %  alltfY, alltfY_pow_warped] =                 ...
% %         my_newcrossf(x, y, frames, tlimits, srate, cycles, ...
% %         'freqs', [3 200], 'nfreqs', 150,                   ...
% %         'timewarp',   EEG_epoched_main.timewarp.latencies(epochs_to_keep,:), ...
% %         'timewarpms', EEG_epoched_main.timewarp.warpto,     ...
% %         'plotamp', 'off', 'plotphase', 'off', 'type', 'coher');
% % 
% % 
% % figure()
% % crosspec_abs_warped = abs(sum(crossspec_warped, 3));
% % crosspec_abs_nowarp = abs(sum(crossspec_noWarp, 3));
% % alltfY_abs = abs(sum(alltfY, 3));
% % alltfX_abs = sum(sum(alltfX, 3));
% % contourf(timesout, freqsout, sum(alltfX_pow_warped, 3), 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 50], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% 
% 
% 
% 
% % figure()
% % plot(freqs, sum(abs(crossspec_warped(:, :, 1)), 2))
% % hold on
% % plot(freqs, sum(abs( alltfX(:, :, 1).*conj(alltfY(:, :, 1)) ) , 2))
% % 
% % 
% % figure()
% % coh_num = alltfX(:, :, 1).*conj(alltfY(:, :, 1));
% % coh_denum = sqrt((abs(alltfX(:, :, 1)).^2) .* (abs(alltfY(:, :, 1)).^2));
% % coh = angle(coh_num ./ coh_denum);
% % contourf(timesout, freqsout, coh, 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 200], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% % 
% % colormap(hot);
% % c = colorbar;
% % ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% % xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % title('No Timewarp', 'FontSize', 18)
% % 
% % 
% % 
% % figure()
% % coh_num = sum(crossspec_warped(:, :, :), 3);
% % % coh_denum = sqrt(alltfX_pow_warped(:, :, 1) .* alltfY_pow_warped(:, :, 1));
% % coh_denum = sqrt(sum(abs(alltfX_pow_warped(:, :, :)), 3) .* sum(abs(alltfY_pow_warped(:, :, :)), 3));
% % coh = abs(coh_num ./ coh_denum);
% % contourf(timesout, freqsout, angle(itcvals), 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 200], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% % 
% % colormap("parula");
% % c = colorbar;
% % ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
% % xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% % title('With Timewarp', 'FontSize', 18)
% 
% 
% 
% 
% 
% 
% 
% %% Let's try something new:  (spoiler alert: it didn't work! I'm tired!)
% % defining sub-phases coherence plots
% % we divide between events into 3 sub phases
% % FlxS to FlxE: early-flexion, mid-flexion, late-flexion
% % ExtS to ExtE: early-extension, mid-extension, late-extension
% 
% % x = LeftPrimMotor_EEG_epochs_CMC;
% % [tmpall_EEG, freqsoutEEG, timesoutEEG, itcvalsEEG] = ...
% %     timefreq( x, 500, ...
% %     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
% %     'freqs', [3, 200], ...
% %     'nfreqs', 150, ...
% %     'freqscale', 'log', ...
% %     'cycles', [3, 0.8]);
% % 
% % y = squeeze(EMG_epochs_CMC(1, :, :));
% % [tmpall_EMG, freqsoutEMG, timesoutEMG, itcvalsEMG] = ...
% %     timefreq( y, 500, ...
% %     'tlimits', [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000, ...
% %     'freqs', [3, 200], ...
% %     'nfreqs', 150, ...
% %     'freqscale', 'log', ...
% %     'cycles', [3, 0.8]);
% % 
% % events_latency = EEG_epoched_main.timewarp.latencies(indx_to_keep, :); % ms
% % % early flexion 
% % t1_k_indx = ones(length(events_latency), 1);
% % t2_k = (events_latency(:, 2) - events_latency(:, 1))/3 + events_latency(:, 1);
% % t2_k_indx = sum(timesoutEMG - t2_k < 0, 2);
% % early_Flx_indx = [t1_k_indx, t2_k_indx];
% % % mid flexion
% % t3_k_indx = t2_k_indx + 1;
% % t4_k = 2*(events_latency(:, 2) - events_latency(:, 1))/3 + events_latency(:, 1);
% % t4_k_indx = sum(timesoutEMG - t4_k < 0, 2);
% % mid_Flx_indx = [t3_k_indx, t4_k_indx];
% % % late flexion
% % t5_k_indx = t4_k_indx + 1;
% % t6_k = events_latency(:, 2);
% % t6_k_indx = sum(timesoutEMG - t6_k < 0, 2);
% % late_Flx_indx = [t5_k_indx, t6_k_indx];
% % 
% % % early extension
% % t7_k_indx = t6_k_indx + 1;
% % t8_k = (events_latency(:, 3) - events_latency(:, 2))/3 + events_latency(:, 2);
% % t8_k_indx = sum(timesoutEMG - t8_k < 0, 2);
% % early_Ext_indx = [t7_k_indx, t8_k_indx];
% % % mid extension
% % t9_k_indx = t8_k_indx + 1;
% % t10_k = 2*(events_latency(:, 3) - events_latency(:, 2))/3 + events_latency(:, 2);
% % t10_k_indx = sum(timesoutEMG - t10_k < 0, 2);
% % mid_Ext_indx = [t9_k_indx, t10_k_indx];
% % % late extension
% % t11_k_indx = t10_k_indx + 1;
% % t12_k = events_latency(:, 3);
% % t12_k_indx = sum(timesoutEMG - t12_k < 0, 2);
% % late_Ext_indx = [t11_k_indx, t12_k_indx];
% % 
% % 
% % crosspec = tmpall_EEG.*conj(tmpall_EMG);
% % autospecEEG = abs(tmpall_EEG).^2;
% % autospecEMG = abs(tmpall_EMG).^2;
% % 
% % % early flexion
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(mean(crosspec(i, early_Flx_indx(i,1):early_Flx_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( mean(autospecEEG(i, early_Flx_indx(i, 1):early_Flx_indx(i, 2), :), 2) , 3) .* ...
% %          sum( mean(autospecEMG(i, early_Flx_indx(i, 1):early_Flx_indx(i, 2), :), 2) , 3));
% % end
% % % coher_early_flexion = abs(numerator ./ denumerator);
% % coher_early_flexion = angle(numerator ./ denumerator);
% % % mid flexion
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(sum(crosspec(i, mid_Flx_indx(i,1):mid_Flx_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( sum(autospecEEG(i, mid_Flx_indx(i, 1):mid_Flx_indx(i, 2), :), 2) , 3) .* ...
% %          sum( sum(autospecEMG(i, mid_Flx_indx(i, 1):mid_Flx_indx(i, 2), :), 2) , 3));
% % end
% % % coher_mid_flexion = abs(numerator ./ denumerator);
% % coher_mid_flexion = angle(numerator ./ denumerator);
% % % late flexion
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(sum(crosspec(i, late_Flx_indx(i,1):late_Flx_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( sum(autospecEEG(i, late_Flx_indx(i, 1):late_Flx_indx(i, 2), :), 2) , 3) .* ...
% %          sum( sum(autospecEMG(i, late_Flx_indx(i, 1):late_Flx_indx(i, 2), :), 2) , 3));
% % end
% % % coher_late_flexion = abs(numerator ./ denumerator);
% % coher_late_flexion = angle(numerator ./ denumerator);
% % % early extension
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(sum(crosspec(i, early_Ext_indx(i,1):early_Ext_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( sum(autospecEEG(i, early_Ext_indx(i, 1):early_Ext_indx(i, 2), :), 2) , 3) .* ...
% %          sum( sum(autospecEMG(i, early_Ext_indx(i, 1):early_Ext_indx(i, 2), :), 2) , 3));
% % end
% % % coher_early_extension = abs(numerator ./ denumerator);
% % coher_early_extension = angle(numerator ./ denumerator);
% % % mid extension
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(sum(crosspec(i, mid_Ext_indx(i,1):mid_Ext_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( sum(autospecEEG(i, mid_Ext_indx(i, 1):mid_Ext_indx(i, 2), :), 2) , 3) .* ...
% %          sum( sum(autospecEMG(i, mid_Ext_indx(i, 1):mid_Ext_indx(i, 2), :), 2) , 3));
% % end
% % % coher_mid_extension = abs(numerator ./ denumerator);
% % coher_mid_extension = angle(numerator ./ denumerator);
% % % late extension
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(sum(crosspec(i, late_Ext_indx(i,1):late_Ext_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( sum(autospecEEG(i, late_Ext_indx(i, 1):late_Ext_indx(i, 2), :), 2) , 3) .* ...
% %          sum( sum(autospecEMG(i, late_Ext_indx(i, 1):late_Ext_indx(i, 2), :), 2) , 3));
% % end
% % % coher_late_extension = abs(numerator ./ denumerator);
% % coher_late_extension = angle(numerator ./ denumerator);
% % 
% % 
% % figure();
% % hold on
% % % coher_early_flexion(coher_early_flexion < 0) = ...
% % %     coher_early_flexion(coher_early_flexion < 0) + 2*pi;
% % plot(freqsout, coher_early_flexion)
% % % coher_mid_flexion(coher_mid_flexion < 0) = ...
% % %     coher_mid_flexion(coher_mid_flexion < 0) + 2*pi;
% % plot(freqsout, coher_mid_flexion)
% % % coher_late_flexion(coher_late_flexion < 0) = ...
% % %     coher_late_flexion(coher_late_flexion < 0) + 2*pi;
% % plot(freqsout, coher_late_flexion)
% % % coher_early_extension(coher_early_extension < 0) = ...
% % %     coher_early_extension(coher_early_extension < 0) + 2*pi;
% % plot(freqsout, coher_early_extension)
% % % coher_mid_extension(coher_mid_extension < 0) = ...
% % %     coher_mid_extension(coher_mid_extension < 0) + 2*pi;
% % plot(freqsout, coher_mid_extension)
% % % coher_late_extension(coher_late_extension < 0) = ...
% % %     coher_late_extension(coher_late_extension < 0) + 2*pi;
% % plot(freqsout, coher_late_extension)
% % 
% % set(gca, 'XScale', 'log')
% % xlim([14, 50])
% % 
% % 
% % 
% % 
% % 
% % % late flexion early extension
% % numerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     numerator(i) = sum(mean(crosspec(i, late_Flx_indx(i,1):early_Ext_indx(i,2), :), 2), 3);
% % end
% % denumerator = zeros(size(crosspec,1),1);
% % for i = 1:size(crosspec,1)
% %     denumerator(i) = ...
% %     sqrt(sum( mean(autospecEEG(i, late_Flx_indx(i, 1):early_Ext_indx(i, 2), :), 2) , 3) .* ...
% %          sum( mean(autospecEMG(i, late_Flx_indx(i, 1):early_Ext_indx(i, 2), :), 2) , 3));
% % end
% % coher_late_flx_early_ext = angle(numerator ./ denumerator);
% % 
% % ylim([0 0.2])
% % 
% % 
% % 
% % coherres_whole = angle( sum(crosspec, 3) ./ sqrt(sum(autospecEEG, 3) .* sum(autospecEMG, 3)));
% % crosspec_amp_warped = abs(sum(crossspec_warped, 3));
% % crosspec_amp_nowarp = abs(sum(crosspec, 3));
% % figure();
% % contourf(timesout, freqsout, crosspec_amp_nowarp, 200, 'linecolor', 'none');
% % set(gca, ...
% %     'yscale',   'log', ...
% %     'ydir',     'norm', ...
% %     'xlim',     [0, 3000], ...
% %     'ylim',     [3, 200], ...
% %     'ytick',    [4 8 14 30 60 100 200], ...
% %     'box',      'on');
% % 
% % colormap("turbo");
% % c = colorbar;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %% check issues!
% % sum(isnan(EMG_epoched_data_intrp(:)))
% % 
% % nan_trials = find(any(isnan(EMG_epoched_data_intrp(2, :, :)), 2));
% % disp(nan_trials)
% 
% 
% 
% 
% % %% check further issues!!
% % figure()
% % e = 176;
% % plot(EMG_time, EMG_epoched.data(1, :, matched_vector(e)))
% % hold on
% % plot(EEG_time, EMG_epoched_data_intrp(1, :, e))
% % title([num2str(e), ' EMG filt'])
% % % plot(EEG_time, 0.001*EEG_epoched_main.icaact(1, :, 1100))
% 
% 
% 
% 
% % %% check the power spectrum of epoched EMGs
% % fs = 500;   % after interpolation to EEG grid
% % 
% % % Example:
% % nan_trials = find(any(isnan(EMG_epoched_data_intrp(2, :, :)), 2));
% % indx_to_keep = ~any(1:size(EMG_epoched_data_intrp, 3) == nan_trials);
% % 
% % EMG_epoched_data_intrp_toplot = EMG_epoched_data_intrp(:, :, indx_to_keep);
% % y = squeeze(EMG_epoched_data_intrp_toplot(1, :, :)); 
% % 
% % win = round(0.5*fs);          % 500 ms
% % noverlap = round(0.25*fs);
% % nfft = max(512, 2^nextpow2(win));
% % 
% % Pyy_all = [];
% % 
% % for tr = 1:size(y,2)
% %     sig = double(y(:,tr));
% %     [Pyy,F] = pwelch(sig, hamming(win), noverlap, nfft, fs);
% %     Pyy_all(:,tr) = Pyy;
% % end
% % 
% % meanPyy = mean(Pyy_all,2);
% % 
% % figure;
% % plot(F,10*log10(meanPyy),'LineWidth',1.5);
% % xlabel('Frequency (Hz)');
% % ylabel('EMG Power (dB)');
% % title('Average EMG power spectrum');
% % xlim([0 60]);
% % grid on;
% 
% 
% 
% 
% % %% check the interpolation procedure
% % epch = 117;
% % figure();
% % hold on
% % plot(EMG_time_new_all(epch, :), ...
% %     EMG_epoched.data(1, :, matched_vector(epch)), 'b');
% % plot(EEG_time, EMG_epoched_data_intrp(1, :, epch), 'r');
% 
% 
% 
% 
% 
% 
% % %% check more - Test 1: PSD before epoching and before interpolation
% % fs_emg = 2000;
% % sig = double(EMG_filt(1,:));   % one EMG channel
% % 
% % [P,F] = pwelch(sig, hamming(1000), 500, 2048, fs_emg);
% % 
% % figure;
% % plot(F,10*log10(P),'LineWidth',1.5);
% % xlim([0 60]);
% % xlabel('Frequency (Hz)');
% % ylabel('Power (dB)');
% % title('Continuous filtered EMG');
% % grid on;
% 
% 
% 
% 
% 
% 
% 
% 
% % %% label early and late stage
% %use function EEG = tag_splitbelt_subconditions(EEG, beginningString, numEpochs, sub_condition_name)
% %%use function EEG = tag_splitbelt_subconditions(EEG, StartString,begOrEnd,numEpochs,epoch_start_delay,sub_condition_name)
% % EEG = tag_exo_subconditions(EEG, 'noExo','end',-30,0,'noExo');
% % EEG = tag_exo_subconditions(EEG, 'unpow','end',-30,0,'unpow');  
% % EEG = tag_exo_subconditions(EEG, 'pow_1','beginning',30,0,'early adapt');
% % EEG = tag_exo_subconditions(EEG, 'pow_3','end',-30,0,'late adapt');
% % EEG = tag_exo_subconditions(EEG, 'deadapt','beginning',30,0,'early post-adapt');
% % EEG = tag_exo_subconditions(EEG, 'deadapt','end',-30,0,'late post-adapt');   
% % 
% % disp('Sub conditions have been labeled in EEG.event')
% % %fill in empty cells so STUDY is happy later
% % empty_ind = find(cellfun(@isempty,{EEG.event.subcond}));
% % [EEG.event(empty_ind).subcond] = deal('none');
% % empty_ind = find(cellfun(@isempty,{EEG.event.cond}));
% % [EEG.event(empty_ind).cond] = deal('none');
% 
% % %% save dataset
% % EEG.filename = string( strcat(filename,'_epoched') );
% % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, CURRENTSET);
% % eeglab redraw;
% % EEG = eeg_checkset(EEG);
% % 
% % disp(EEG.filename)
% % EEG = pop_saveset(EEG);



