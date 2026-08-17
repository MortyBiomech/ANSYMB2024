%% run_stats_selected.m
%
% Phase 2: run cluster-based permutation statistics only on conditions
% that look promising after visual inspection of the phase-1 CMC figures.
%
% Fill in the 'selected_labels' cell array below with the condition labels
% you want to test (copied from the figure filenames, without the .png).
% Everything else (data loading, TF recomputation) is identical to the
% phase-1 script — only the CMC function changes.

% ── Conditions to test ────────────────────────────────────────────────────
selected_labels = { ...
    'sub18__bp20_450__abs_rect__timenorm_upsample__L_PrimMotor', ...
    'sub18__bp20_200__hilbert_demod__locked_FlxS__R_PrimMotor',  ...
    % add more here ...
};

% ── Paths and setup (same as run_CMC_conditions.m) ────────────────────────
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\LSL\xdf-Matlab-master'));

EEGLAB_path      = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';
data_path        = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path  = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
rawdata_path     = [data_path, '0_source_data\'];
rawEEGLAB_path   = [data_path, '2_raw-EEGLAB\'];
epoched_EEG_path = [data_path, '5_single-subject-EEG-analysis\', ...
                    'timewarp_test\Epoched_data'];
results_path     = [data_path, '5_single-subject-EEG-analysis\CMC_conditions\'];
stats_path       = [results_path, 'stats\'];
if ~exist(stats_path, 'dir'), mkdir(stats_path); end

this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

subject_id = 18;
emg_ids    = [2, 3, 4, 5];
emg_labels = {'VL', 'RF', 'GM', 'BF'};
emg_srate  = 2000;
eeg_srate  = 500;
epoch_limits_s  = [-0.5, 3.5];
epoch_limits_ms = epoch_limits_s * 1000;

% ── Parse condition labels into their components ──────────────────────────
% Label format: sub{id}__{filter}__{proc}__{tf_method}__{cluster}
filter_conds = struct( ...
    'label',  {'bp20_450',   'bp20_200',   'hp20',    'hp10',    'raw'}, ...
    'type',   {'bandpass',   'bandpass',   'highpass','highpass','none'}, ...
    'cutoff', {[20 450],     [20 200],     20,         10,        []}    );
proc_conds = struct( ...
    'label', {'abs_rect', 'hilbert_amp', 'hilbert_demod', 'none'});
tf_methods   = {'timenorm_upsample', 'locked_FlxS', 'timewarp_TF_posthoc'};
eeg_clusters = {'L_PrimMotor', 'R_PrimMotor', 'L_ParOcc', 'R_ParOcc'};

% ── Load EEG data and IC activations (same as phase 1) ───────────────────
fprintf('Loading EEG datasets...\n');
cd([rawEEGLAB_path, 'sub-', num2str(subject_id)])
EEG_orig = pop_loadset('filename', ['sub-', num2str(subject_id), '_merged_EEG.set']);
cd(this_path)
cd(epoched_EEG_path)
EEG_epoched_main = pop_loadset('filename', ...
    ['sub-', num2str(subject_id), '_cleaned_with_ICA_epoched.set']);
cd(this_path)

EEG_times    = EEG_epoched_main.times;
tlimits_tw   = [EEG_epoched_main.xmin, EEG_epoched_main.xmax] * 1000;
tw_latencies = EEG_epoched_main.timewarp.latencies;

output  = runs_concatenated(subject_id, rawdata_path);
EMG_raw = double(output.All_EMG(emg_ids, :));

% IC cluster activations
load('Subjects_ICs_in_clusters.mat')
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
    EEG_IC.(cfield) = squeeze(EEG_epoched_main.icaact(ic_id, :, :));
end

% ── Epoch timing metadata (same once-only extraction) ────────────────────
[b_meta, a_meta] = butter(4, [20 450] / (emg_srate/2), 'bandpass');
EMG_meta = filtfilt(b_meta, a_meta, double(EMG_raw'))';
eeg_emptyset = structfun(@(x) [], EEG_orig, 'UniformOutput', false);
EMG_tmpl = eeg_emptyset;
EMG_tmpl.data = double(EMG_meta);
EMG_tmpl.nbchan = size(EMG_meta, 1); EMG_tmpl.pnts = size(EMG_meta, 2);
EMG_tmpl.trials = 1; EMG_tmpl.srate = emg_srate;
EMG_tmpl.xmin = 0; EMG_tmpl.xmax = (EMG_tmpl.pnts-1)/emg_srate;
EMG_tmpl.setname = 'EMG_continuous'; EMG_tmpl.filename = ''; EMG_tmpl.filepath = '';
for ch = 1:EMG_tmpl.nbchan, EMG_tmpl.chanlocs(ch).labels = emg_labels{ch}; end
EMG_tmpl = eeg_checkset(EMG_tmpl);

[EMG_ev, ~] = import_same_event_on_EMG_stream( ...
    EMG_tmpl, EEG_orig, output, subject_id, processing_path, data_path);
for e = 1:length(EMG_ev.event)
    nums = str2double(strsplit(EMG_ev.event(e).desc, '_')).';
    if strcmp(EMG_ev.event(e).type, 'boundary')
        EMG_ev.event(e).trial = 'none'; EMG_ev.event(e).cond = 'none';
    else
        EMG_ev.event(e).trial = nums(end); EMG_ev.event(e).cond = nums(2);
    end
end
EMG_ep_src = pop_selectevent(EMG_ev, 'type', {'FlxS','ExtS','ExtE'}, 'deleteevents', 'on');
EMG_epoched_meta = pop_epoch(EMG_ep_src, {'FlxS'}, epoch_limits_s, ...
    'newname', 'meta', 'epochinfo', 'yes');
EMG_epoched_meta = eeg_checkset(EMG_epoched_meta);

[matched_vector, ~] = find_matched_epochs(EEG_epoched_main, EMG_epoched_meta);
EMG_time_meta = EMG_epoched_meta.times;
N_EMG_smp     = EMG_epoched_meta.pnts;
delta_t_ms    = NaN(length(matched_vector), 1);
epoch_start_s = NaN(size(EMG_epoched_meta.epoch));

for i = 1:length(matched_vector)
    if isnan(matched_vector(i)), continue; end
    ep = EMG_epoched_meta.epoch(matched_vector(i));
    FlxS_idx = find(strcmpi(ep.eventtype,'FlxS') & cell2mat(ep.eventlatency)==0, 1);
    if isempty(FlxS_idx), continue; end
    t_last = ep.eventinit_time_lastEMGpoint{FlxS_idx};
    t_FlxS = ep.eventinit_time{FlxS_idx};
    delta_t_ms(i) = 1000*(t_FlxS - t_last);
    epoch_start_s(matched_vector(i)) = t_FlxS + epoch_limits_s(1);
end
epoch_start_smp = round(epoch_start_s * emg_srate) + 1;
epoch_end_smp   = epoch_start_smp + N_EMG_smp - 1;

timewarp_indices_all = round(eeg_lat2point(reshape(tw_latencies,1,[]), ...
    1, eeg_srate, tlimits_tw, 1e-3));
timewarp_indices_all = reshape(timewarp_indices_all, size(tw_latencies));

tf_params.srate = eeg_srate; tf_params.freqs = [3 60];
tf_params.nfreqs = 60; tf_params.cycles = [2 0.8];
tf_params.freqscale = 'log'; tf_params.ntimesout = 200;

cmc_opts.nShuffles   = 500;
cmc_opts.alpha_pix   = 0.05;
cmc_opts.alpha_clust = 0.05;
cmc_opts.verbose     = true;
cmc_opts.plot        = true;

% ── Run stats on selected conditions ─────────────────────────────────────
fprintf('\nRunning cluster permutation on %d selected conditions...\n\n', ...
    numel(selected_labels));

for s = 1:numel(selected_labels)
    lbl = selected_labels{s};
    fprintf('[%d/%d] %s\n', s, numel(selected_labels), lbl);

    % ── Parse label → condition components ───────────────────────────────
    parts    = strsplit(lbl, '__');
    % parts: {sub_id, filter, proc, tf_method, cluster}
    fi = find(strcmp({filter_conds.label}, parts{2}));
    pi = find(strcmp({proc_conds.label},   parts{3}));
    ti = find(strcmp(tf_methods,           parts{4}));
    ci = find(strcmp(eeg_clusters,         parts{5}));

    if any([isempty(fi), isempty(pi), isempty(ti), isempty(ci)])
        warning('Could not parse label: %s — skipping.', lbl);
        continue;
    end

    fig_path = fullfile(stats_path, [lbl, '_clusterperm.png']);
    if exist(fig_path, 'file')
        fprintf('  Already done — skipping.\n');
        continue;
    end

    % ── Recompute EMG epochs for this filter / proc combo ────────────────
    EMG_filt = apply_emg_filter(EMG_raw, emg_srate, filter_conds(fi));
    EMG_proc = apply_emg_processing(EMG_filt, proc_conds(pi).label);

    N_EEG = length(EEG_epoched_main.epoch);
    EMG_intrp = NaN(size(EMG_raw, 1), length(EEG_times), N_EEG);
    for i = 1:N_EEG
        mi = matched_vector(i);
        if isnan(mi) || isnan(delta_t_ms(i)), continue; end
        s1 = epoch_start_smp(mi); s2 = epoch_end_smp(mi);
        if s1 < 1 || s2 > size(EMG_proc, 2), continue; end
        emg_ep = EMG_proc(:, s1:s2);
        emg_t  = EMG_time_meta - delta_t_ms(i);
        EMG_intrp(:, :, i) = interp1(emg_t', emg_ep', EEG_times', 'linear')';
    end
    nan_eps = squeeze(any(any(isnan(EMG_intrp), 1), 2));
    keep    = ~nan_eps;

    EMG_epochs_CMC = EMG_intrp(:, :, keep);
    tw_idx         = timewarp_indices_all(keep, :);
    X_epochs       = EEG_IC.(eeg_clusters{ci})(:, keep);
    Y_epochs       = EMG_epochs_CMC;

    % ── Recompute TF ─────────────────────────────────────────────────────
    try
        switch tf_methods{ti}
            case 'timenorm_upsample'
                [TF_X, TF_Y, time_axis, freqs] = compute_TF_timenorm( ...
                    X_epochs, Y_epochs, tw_idx, tf_params);
            case 'locked_FlxS'
                [TF_X, TF_Y, time_axis, freqs] = compute_TF_locked( ...
                    X_epochs, Y_epochs, tf_params, epoch_limits_ms);
            case 'timewarp_TF_posthoc'
                [TF_X, TF_Y, time_axis, freqs] = compute_TF_timewarp_posthoc( ...
                    X_epochs, Y_epochs, tw_idx, tf_params, epoch_limits_ms);
        end

        % ── Cluster permutation ───────────────────────────────────────────
        cmc_opts.title         = strrep(lbl, '_', '\_');
        cmc_opts.muscle_labels = emg_labels;

        [~, ~, ~, cluster_info, fig_h] = compute_and_plot_CMC_cluster_perm( ...
            TF_X, TF_Y, time_axis, freqs, cmc_opts);

        if ~isempty(fig_h) && isvalid(fig_h)
            try
                exportgraphics(fig_h, fig_path, 'Resolution', 150);
            catch
                saveas(fig_h, fig_path, 'png');
            end
            close(fig_h);
            fprintf('  Saved: %s\n', fig_path);
        end

        % Print cluster summary
        for m = 1:numel(cluster_info)
            if isempty(cluster_info{m}), continue; end
            fprintf('  Muscle %d: %d significant clusters\n', m, numel(cluster_info{m}));
            for c_i = 1:numel(cluster_info{m})
                ci_info = cluster_info{m}(c_i);
                fprintf('    cluster %d: %.1f–%.1f Hz, %.1f–%.1f%%, p=%.3f\n', ...
                    c_i, ci_info.freq_range(1), ci_info.freq_range(2), ...
                    ci_info.time_range(1), ci_info.time_range(2), ci_info.p_value);
            end
        end

    catch ME
        warning('Stats failed for %s:\n  %s', lbl, ME.message);
    end
end

fprintf('\nDone. Stats figures saved to:\n  %s\n', stats_path);