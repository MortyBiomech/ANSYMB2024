%% run_CMC_conditions.m
%
% Systematically evaluates cortico-muscular coherence (CMC) across all
% combinations of EMG preprocessing options and TF normalisation methods,
% as laid out in the design diagram.
%
% Conditions:
%   Filter     : bp 20-450 | bp 20-200 | hp 20 | hp 10 | raw (no filter)
%   EMG proc   : abs rectification | Hilbert amplitude | Hilbert demod | none
%   TF method  : timenorm_upsample | locked_FlxS | timewarp_TF_posthoc
%   EEG cluster: L Prim Motor | R Prim Motor | L Parieto-Occipital | R Parieto-Occipital
%
% Strategy
% ─────────
%   import_same_event_on_EMG_stream() is called ONCE on a neutral (bandpass)
%   EMG signal to extract event timestamps and per-epoch timing metadata.
%   All subsequent conditions reuse this metadata, re-applying only the
%   filter + processing step to the continuous raw EMG before interpolating
%   each epoch onto the EEG time grid.
%
% Prerequisites on MATLAB path:
%   compute_TF_timenorm.m              (upsampling / cycle-% method)
%   compute_TF_locked.m                (new — provided alongside this file)
%   compute_TF_timewarp_posthoc.m      (new — provided alongside this file)
%   compute_and_plot_CMC_cluster_perm.m
%   runs_concatenated.m
%   import_same_event_on_EMG_stream.m
%   EEGLAB (eeglab2026.0.0 or later)

clc; clear;

% ── Paths ────────────────────────────────────────────────────────────────
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\LSL\xdf-Matlab-master'));

EEGLAB_path      = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';
data_path        = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path  = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
rawdata_path     = [data_path, '0_source_data\'];
rawEEGLAB_path   = [data_path, '2_raw-EEGLAB\'];
epoched_EEG_path = [data_path, '5_single-subject-EEG-analysis\', ...
                    'timewarp_test\Epoched_data'];
cmc_folder_path  = [processing_path, 'Group_Level_PostProcessing\Final_paper_plot_generation\_NHB\Corticomuscular_Coherence\']
results_path     = [cmc_folder_path, 'CMC_results\CMC_conditions\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

% ── EEGLAB ───────────────────────────────────────────────────────────────
this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

% ── Subject ───────────────────────────────────────────────────────────────
subject_id = 18;
emg_ids    = [2, 3, 4, 5];
emg_labels = {'VL', 'RF', 'GM', 'BF'};
emg_srate  = 2000;
eeg_srate  = 500;

epoch_limits_s  = [-0.5, 3.5];          % pop_epoch window in seconds
epoch_limits_ms = epoch_limits_s * 1000; % ms (for timefreq tlimits)

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 1 — Load EEG datasets
% ═════════════════════════════════════════════════════════════════════════
fprintf('Loading EEG datasets...\n');

cd([rawEEGLAB_path, 'sub-', num2str(subject_id)])
EEG_orig = pop_loadset('filename', ...
    ['sub-', num2str(subject_id), '_merged_EEG.set']);
cd(this_path)

cd(epoched_EEG_path)
EEG_epoched_main = pop_loadset('filename', ...
    ['sub-', num2str(subject_id), '_cleaned_with_ICA_epoched.set']);
cd(this_path)

EEG_times = EEG_epoched_main.times;   % [1 x nSamplesEEG] in ms
N_EEG     = EEG_epoched_main.nbchan;

% ── Timewarp metadata (shared across all conditions) ──────────────────────
tlimits_tw = [EEG_epoched_main.xmin, EEG_epoched_main.xmax] * 1000; % ms
tw_latencies_all = EEG_epoched_main.timewarp.latencies; % [nEpochs x nEvents]

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 2 — Load raw EMG and run the timing-metadata extraction ONCE
% ═════════════════════════════════════════════════════════════════════════
fprintf('Loading raw EMG and computing epoch timing metadata...\n');

output  = runs_concatenated(subject_id, rawdata_path);
EMG_raw = double(output.All_EMG(emg_ids, :));   % [4 x nSamples]
nChan   = size(EMG_raw, 1);

% --- Apply a neutral bandpass filter for the metadata pass only ---
[b_meta, a_meta] = butter(4, [20, 450] / (emg_srate/2), 'bandpass');
EMG_meta = filtfilt(b_meta, a_meta, EMG_raw')';

% --- Build EEGLAB-like EMG structure template ---
eeg_emptyset = structfun(@(x) [], EEG_orig, 'UniformOutput', false);
EMG_tmpl               = eeg_emptyset;
EMG_tmpl.data          = double(EMG_meta);
EMG_tmpl.nbchan        = nChan;
EMG_tmpl.pnts          = size(EMG_meta, 2);
EMG_tmpl.trials        = 1;
EMG_tmpl.srate         = emg_srate;
EMG_tmpl.xmin          = 0;
EMG_tmpl.xmax          = (EMG_tmpl.pnts - 1) / emg_srate;
EMG_tmpl.setname       = 'EMG_continuous';
EMG_tmpl.filename      = '';
EMG_tmpl.filepath      = '';
for ch = 1:nChan
    EMG_tmpl.chanlocs(ch).labels = emg_labels{ch};
end
EMG_tmpl = eeg_checkset(EMG_tmpl);

% --- Import events (run only once; events contain all timing metadata) ---
[EMG_events_struct, ~] = import_same_event_on_EMG_stream( ...
    EMG_tmpl, EEG_orig, output, subject_id, processing_path, data_path);

% Add trial / cond fields
for e = 1:length(EMG_events_struct.event)
    nums = str2double(strsplit(EMG_events_struct.event(e).desc, '_')).';
    if strcmp(EMG_events_struct.event(e).type, 'boundary')
        EMG_events_struct.event(e).trial = 'none';
        EMG_events_struct.event(e).cond  = 'none';
    else
        EMG_events_struct.event(e).trial = nums(end);
        EMG_events_struct.event(e).cond  = nums(2);
    end
end

% --- Epoch using the stored event structure ---
EMG_for_epoch = EMG_events_struct;
EMG_for_epoch = pop_selectevent(EMG_for_epoch, ...
    'type', {'FlxS','ExtS','ExtE'}, 'deleteevents', 'on');
EMG_epoched_meta = pop_epoch(EMG_for_epoch, {'FlxS'}, epoch_limits_s, ...
    'newname', 'FlxStart meta', 'epochinfo', 'yes');
EMG_epoched_meta = eeg_checkset(EMG_epoched_meta);

% --- Find matched EEG-EMG epoch pairs ---
[matched_vector, condition_vector] = find_matched_epochs( ...
    EEG_epoched_main, EMG_epoched_meta);

% --- Pre-extract per-epoch timing metadata (delta_t, epoch windows) ---
% delta_t corrects for sub-sample offset between EEG and EMG clocks.
% It is independent of EMG signal content — valid for all conditions.
EMG_time_meta    = EMG_epoched_meta.times;   % [1 x nSamplesEMGepoch] ms
N_EMG_ep_samples = EMG_epoched_meta.pnts;

delta_t_ms     = NaN(length(matched_vector), 1);
epoch_start_s  = NaN(size(EMG_epoched_meta.epoch));  % FlxS abs time in seconds

for i = 1:length(matched_vector)
    if isnan(matched_vector(i)), continue; end
    ep_info = EMG_epoched_meta.epoch(matched_vector(i));

    % Find the FlxS event at latency = 0
    FlxS_idx = find( ...
        strcmpi(ep_info.eventtype, 'FlxS') & ...
        cell2mat(ep_info.eventlatency) == 0, 1);

    if isempty(FlxS_idx), continue; end

    t_last  = ep_info.eventinit_time_lastEMGpoint{FlxS_idx};
    t_FlxS  = ep_info.eventinit_time{FlxS_idx};
    delta_t_ms(i)    = 1000 * (t_FlxS - t_last);   % ms, sub-sample correction
    epoch_start_s(matched_vector(i)) = t_FlxS + epoch_limits_s(1); % abs start
end

% Convert epoch start times to continuous EMG sample indices
epoch_start_smp = round(epoch_start_s * emg_srate) + 1;  % 1-based
epoch_end_smp   = epoch_start_smp + N_EMG_ep_samples - 1;

% ── Timewarp indices (sample indices at EEG srate, nEpochs × nEvents) ────
% These are used by compute_TF_timenorm and compute_TF_timewarp_posthoc.
timewarp_indices_all = round(eeg_lat2point( ...
    reshape(tw_latencies_all, 1, []), 1, eeg_srate, tlimits_tw, 1e-3));
timewarp_indices_all = reshape(timewarp_indices_all, size(tw_latencies_all));

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 3 — Load IC cluster activations
% ═════════════════════════════════════════════════════════════════════════
load('Subjects_ICs_in_clusters.mat')   % -> SUBJECTS_ICS

cluster_defs = { ...
    'Left_Prim_Motor',        'L_PrimMotor'; ...
    'Right_Prim_Motor',       'R_PrimMotor'; ...
    'Left_Parieto_Occipital', 'L_ParOcc';    ...
    'Right_Parieto_Occipital','R_ParOcc'     };

EEG_IC = struct();
for c = 1:size(cluster_defs, 1)
    cname  = cluster_defs{c, 1};
    cfield = cluster_defs{c, 2};
    cidx   = find(strcmpi(SUBJECTS_ICS(:, 1), cname));
    sidx   = find(SUBJECTS_ICS{cidx, 2}.Subjects == subject_id - 4);
    ic_id  = SUBJECTS_ICS{cidx, 2}.ICs(sidx);
    % Store as [nSamples x nEpochs]
    EEG_IC.(cfield) = squeeze(EEG_epoched_main.icaact(ic_id, :, :));
end

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 4 — Define conditions
% ═════════════════════════════════════════════════════════════════════════

% 1. EMG filter conditions
filter_conds = struct( ...
    'label',  {'bp20_450',   'bp20_200',   'hp20',    'hp10',    'raw'}, ...
    'type',   {'bandpass',   'bandpass',   'highpass','highpass','none'}, ...
    'cutoff', {[20 450],     [20 200],     20,         10,        []}    );

% 2. EMG signal processing (applied after filtering)
proc_conds = struct( ...
    'label', {'abs_rect',             'hilbert_amp',       'hilbert_demod',         'none'}, ...
    'desc',  {'Abs rectification',    'Hilbert amplitude', 'Hilbert demodulation',  'No processing'} );

% 3. TF time-normalisation / CMC method
%    'timenorm_upsample'   – Fauvet-style: upsample to max SOI, cycle-% axis
%    'locked_FlxS'         – no timewarp; epochs locked to FlxS
%    'timewarp_TF_posthoc' – compute TF at native srate, then warp TF matrix
%                            to common %-grid post-hoc (per-epoch interpolation)
%
%  TODO: two further methods from the diagram are not yet implemented here:
%    'timewarp_crossspectrum' – compute cross-/auto-spectra per trial, then
%                               timewarp the spectra before summing.
%    're_epoch'               – re-epoch around each of ExtS, FlxS, ExtE
%                               without timewarping and compute separate CMC.
tf_methods = {'timenorm_upsample', 'locked_FlxS', 'timewarp_TF_posthoc'};

% 4. EEG IC clusters
eeg_clusters = {'L_PrimMotor', 'R_PrimMotor', 'L_ParOcc', 'R_ParOcc'};

% ── Shared TF and CMC parameters ─────────────────────────────────────────
tf_params.srate     = eeg_srate;
tf_params.freqs     = [3 60];
tf_params.nfreqs    = 60;
tf_params.cycles    = [2 0.8];
tf_params.freqscale = 'log';
tf_params.ntimesout = 200;

cmc_opts.nShuffles   = 500;
cmc_opts.alpha_pix   = 0.05;
cmc_opts.alpha_clust = 0.05;
cmc_opts.verbose     = false;
cmc_opts.plot        = true;

%% ════════════════════════════════════════════════════════════════════════
%  PHASE 5 — Main condition loop
% ═════════════════════════════════════════════════════════════════════════
total_conds = numel(filter_conds) * numel(proc_conds) * ...
              numel(tf_methods)   * numel(eeg_clusters);
fprintf('\nTotal conditions: %d\n\n', total_conds);
run_count = 0;

for fi = 1:numel(filter_conds)
for pi = 1:numel(proc_conds)

    %% ── 5a. Apply EMG filter ────────────────────────────────────────────
    EMG_filt = apply_emg_filter(EMG_raw, emg_srate, filter_conds(fi));

    %% ── 5b. Apply EMG processing (rectification / demodulation) ─────────
    EMG_proc = apply_emg_processing(EMG_filt, proc_conds(pi).label);

    %% ── 5c. Extract and interpolate EMG epochs to EEG time grid ─────────
    %   Reuse pre-computed epoch windows and delta_t; no need to re-epoch
    %   via EEGLAB — just slice from the new continuous processed signal.
    N_EEG_epochs = length(EEG_epoched_main.epoch);
    EMG_intrp = NaN(nChan, length(EEG_times), N_EEG_epochs);

    for i = 1:N_EEG_epochs
        mi = matched_vector(i);
        if isnan(mi) || isnan(delta_t_ms(i)), continue; end

        s1 = epoch_start_smp(mi);
        s2 = epoch_end_smp(mi);
        if s1 < 1 || s2 > size(EMG_proc, 2), continue; end

        emg_ep   = EMG_proc(:, s1:s2);           % [4 x nSamplesEMG]
        emg_time = EMG_time_meta - delta_t_ms(i); % shift for sub-sample alignment

        EMG_intrp(:, :, i) = interp1(emg_time', emg_ep', EEG_times', 'linear')';
    end

    % Identify valid (non-NaN) EEG epochs after matching & interpolation
    nan_eps    = any(any(isnan(EMG_intrp), 1), 2);
    nan_eps    = squeeze(nan_eps);             % [N_EEG_epochs x 1] logical
    keep       = ~nan_eps;                     % logical index
    keep_idx   = find(keep);                  % numeric indices

    EMG_epochs_CMC = EMG_intrp(:, :, keep);   % [4 x nSampEEG x nValidEpochs]
    tw_idx         = timewarp_indices_all(keep, :); % [nValid x nEvents]

    fprintf('  %s + %s: %d / %d valid epochs\n', ...
        filter_conds(fi).label, proc_conds(pi).label, ...
        sum(keep), N_EEG_epochs);

    %% ── 5d. TF method × EEG cluster loop ────────────────────────────────
    for ti = 1:numel(tf_methods)
    for ci = 1:numel(eeg_clusters)

        run_count = run_count + 1;
        cname       = eeg_clusters{ci};
        tfm         = tf_methods{ti};
        cond_label  = sprintf('sub%02d__%s__%s__%s__%s', ...
            subject_id, filter_conds(fi).label, proc_conds(pi).label, tfm, cname);

        fprintf('  [%d/%d] %s\n', run_count, total_conds, cond_label);

        fig_path = fullfile(results_path, [cond_label, '.png']);
        if exist(fig_path, 'file')
            fprintf('    already saved — skipping.\n');
            continue;
        end

        % EEG IC activations for this cluster [nSamples x nEpochs_valid]
        X_epochs = EEG_IC.(cname)(:, keep);   % [nSamples x nValid]
        Y_epochs = EMG_epochs_CMC;             % [4 x nSamples x nValid]

        try
            %% ── Compute TF ───────────────────────────────────────────────
            switch tfm

                case 'timenorm_upsample'
                    % Upsample each epoch to the max-SOI srate, express time
                    % as cycle %, and call timefreq with the same % output
                    % grid for every epoch (existing compute_TF_timenorm).
                    [TF_X, TF_Y, time_axis, freqs] = compute_TF_timenorm( ...
                        X_epochs, Y_epochs, tw_idx, tf_params);

                case 'locked_FlxS'
                    % No time-normalisation. All epochs share the original
                    % epoch window; the x-axis is time in ms.
                    % NOTE: the CMC plot x-axis label will read "Cycle (%)"
                    % but the values are actually milliseconds.
                    [TF_X, TF_Y, time_axis, freqs] = compute_TF_locked( ...
                        X_epochs, Y_epochs, tf_params, epoch_limits_ms);

                case 'timewarp_TF_posthoc'
                    % Compute TF at native srate for every epoch, then warp
                    % each complex TF matrix to the common cycle-% grid by
                    % linear interpolation of real and imaginary parts.
                    [TF_X, TF_Y, time_axis, freqs] = ...
                        compute_TF_timewarp_posthoc( ...
                            X_epochs, Y_epochs, tw_idx, tf_params, epoch_limits_ms);
            end

            %% ── Compute CMC and plot (no stats — phase 1 only) ───────────────
            cmc_opts.title         = strrep(cond_label, '_', '\_');
            cmc_opts.muscle_labels = emg_labels;
            
            [coh, fig_h] = compute_and_plot_CMC( ...
                TF_X, TF_Y, time_axis, freqs, cmc_opts);
            
            %% ── Save figure and raw MSC matrix ──────────────────────────────
            if ~isempty(fig_h) && isvalid(fig_h)
                try
                    exportgraphics(fig_h, fig_path, 'Resolution', 150);
                catch
                    saveas(fig_h, fig_path, 'png');
                end
                close(fig_h);
            end
            
            % Save the small MSC matrix alongside the figure for later review
            mat_path = fullfile(results_path, [cond_label, '.mat']);
            save(mat_path, 'coh', 'time_axis', 'freqs', 'cond_label', '-v7.3');

        catch ME
            warning('Condition %s failed:\n  %s', cond_label, ME.message);
        end

    end % cluster loop
    end % TF method loop

end % proc loop
end % filter loop

fprintf('\nAll %d conditions done.\nResults in: %s\n', total_conds, results_path);


% ══════════════════════════════════════════════════════════════════════════
%  LOCAL HELPER FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

% function EMG_filt = apply_emg_filter(EMG_raw, fs, fcond)
% % APPLY_EMG_FILTER  Apply band-pass, high-pass, or no filter to raw EMG.
% %   Input / output: [nChan x nSamples]
%     switch fcond.type
%         case 'bandpass'
%             [b, a]   = butter(4, fcond.cutoff / (fs/2), 'bandpass');
%             EMG_filt = filtfilt(b, a, double(EMG_raw'))';
%         case 'highpass'
%             [b, a]   = butter(4, fcond.cutoff / (fs/2), 'high');
%             EMG_filt = filtfilt(b, a, double(EMG_raw'))';
%         case 'none'
%             EMG_filt = double(EMG_raw);
%         otherwise
%             error('apply_emg_filter: unknown filter type "%s".', fcond.type);
%     end
% end


% function EMG_proc = apply_emg_processing(EMG_filt, proc_label)
% % APPLY_EMG_PROCESSING  Rectify or demodulate filtered EMG.
% %   The Hilbert transform is applied along the time dimension.
% %   Input / output: [nChan x nSamples]
%     switch proc_label
%         case 'abs_rect'
%             % Full-wave rectification
%             EMG_proc = abs(EMG_filt);
% 
%         case 'hilbert_amp'
%             % Instantaneous amplitude: A(t) = |z(t)|  where z = analytic signal
%             z        = hilbert(double(EMG_filt'))';   % hilbert expects columns
%             EMG_proc = abs(z);
% 
%         case 'hilbert_demod'
%             % Phase-only (demodulated): cos(phi(t)) = Re{ z(t) / |z(t)| }
%             % Removes amplitude; unit amplitude, carries phase information only.
%             z        = hilbert(double(EMG_filt'))';
%             EMG_proc = real(z ./ abs(z));
% 
%         case 'none'
%             EMG_proc = double(EMG_filt);
% 
%         otherwise
%             error('apply_emg_processing: unknown proc label "%s".', proc_label);
%     end
% end


% function [matched_vector, condition_vector] = find_matched_epochs( ...
%         EEG_epoched_main, EMG_epoched)
% % FIND_MATCHED_EPOCHS  Replicate the epoch-matching logic from the main script.
% %   Returns matched_vector [nEEGEpochs x 1] (EMG epoch index per EEG epoch,
% %   NaN if no match) and condition_vector [nEEGEpochs x 1].
% 
%     % --- collect FlxS init_time for each EMG epoch ---
%     eventtypes_EMG   = {EMG_epoched.epoch.eventtype};
%     eventlatency_EMG = cellfun(@(x) cell2mat(x), ...
%         {EMG_epoched.epoch.eventlatency}, 'UniformOutput', false);
% 
%     FlxS_type_lgl = cellfun(@(x) strcmpi(x, 'FlxS'), eventtypes_EMG, ...
%         'UniformOutput', false);
%     FlxS_lat_lgl  = cellfun(@(x) x == 0, eventlatency_EMG, ...
%         'UniformOutput', false);
%     FlxS_idx_EMG  = cellfun(@(x, y) find(x & y), ...
%         FlxS_type_lgl, FlxS_lat_lgl, 'UniformOutput', false);
% 
%     more_one = cellfun(@(x) numel(x) > 1, FlxS_idx_EMG);
% 
%     init_time_cell = cellfun(@(x) cell2mat(x), ...
%         {EMG_epoched.epoch.eventinit_time_EEG}, 'UniformOutput', false);
%     init_time_EMG  = cellfun(@(x, y) x(y(1)), init_time_cell, FlxS_idx_EMG, ...
%         'UniformOutput', false);
%     init_time_EMG  = cell2mat(init_time_EMG);
% 
%     N_EEG = length(EEG_epoched_main.epoch);
%     matched_vector   = NaN(N_EEG, 1);
%     condition_vector = zeros(N_EEG, 1);
% 
%     for i = 1:N_EEG
%         evtypes = EEG_epoched_main.epoch(i).eventtype;
%         evlat   = cell2mat(EEG_epoched_main.epoch(i).eventlatency);
%         FlxS_i  = find(strcmpi(evtypes, 'FlxS') & evlat == 0, 1);
%         if isempty(FlxS_i), continue; end
% 
%         t_EEG = EEG_epoched_main.epoch(i).eventinit_time{FlxS_i};
%         hit   = find(abs(init_time_EMG - t_EEG) < 1e-10);
% 
%         if numel(hit) > 1 || any(more_one(hit)) || isempty(hit)
%             matched_vector(i) = NaN;
%         else
%             matched_vector(i)   = hit;
%             condition_vector(i) = str2double( ...
%                 EEG_epoched_main.epoch(i).eventcond{FlxS_i});
%         end
%     end
% end
