%% run_CMC_reepoch.m
%
% Event-locked cortico-muscular coherence — no time normalisation.
% Re-epochs around FlxS / ExtS / ExtE using the CONTINUOUS ICA EEG,
% giving a symmetric [-1, 1] s window for all three events.
%
% Why continuous EEG:
%   EEG_epoched_main was created with a 500 ms pre-epoch buffer around FlxS.
%   Using it for sub-epoching limits FlxS to [-0.5, 1.5] s.
%   The continuous ICA dataset has no such constraint — any symmetric
%   window is possible as long as it falls within the recording bounds.
%
% EEG source  : sub-{id}_cleaned_with_ICA.set  (continuous, ICA weights stored)
%               icaact computed on-the-fly from icaweights * icasphere * data
% EMG filter  : band-pass 20-450 Hz
% EMG proc    : (1) Hilbert amplitude  |H{x}|
%               (2) Hilbert demodulation  cos(phi(t))
% EEG clusters: L/R primary motor, L/R parieto-occipital

clc; clear;

%% ── Subject list ─────────────────────────────────────────────────────────
subject_ids = [5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18];
emg_ids     = [2, 3, 4, 5];
n_subjects  = numel(subject_ids);

%% ── Paths ────────────────────────────────────────────────────────────────
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\LSL\xdf-Matlab-master'));

EEGLAB_path      = 'D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0';
data_path        = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
processing_path  = 'D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\';
rawdata_path     = [data_path, '0_source_data\'];
rawEEGLAB_path   = [data_path, '2_raw-EEGLAB\'];
epoched_EEG_path = [data_path, '5_single-subject-EEG-analysis\', ...
                    'timewarp_test\Epoched_data'];
cont_ICA_path    = [data_path, '5_single-subject-EEG-analysis\', ...
                    'timewarp_test\'];             % continuous ICA EEG
cmc_folder_path  = [processing_path, 'Group_Level_PostProcessing\', ...
                    'Final_paper_plot_generation\_NHB\', ...
                    'Corticomuscular_Coherence\run_all_cases\'];
results_path     = [cmc_folder_path, ...
                    'CMC_results\re_epoch_notimewarp\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

%% ── Initialise EEGLAB ────────────────────────────────────────────────────
this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

%% ── Processing conditions ────────────────────────────────────────────────
proc_conditions = struct( ...
    'label', {'hilbert_amp',              'hilbert_demod'         }, ...
    'desc',  {'Hilbert amplitude |H{x}|', 'Hilbert demod. cos(phi)'});
n_proc = numel(proc_conditions);

%% ── EEG cluster definitions ──────────────────────────────────────────────
% eeg_clusters = {'L_PrimMotor', 'R_PrimMotor', 'L_ParOcc', 'R_ParOcc'};
% cluster_defs = { ...
%     'Left_Prim_Motor',        'L_PrimMotor'; ...
%     'Right_Prim_Motor',       'R_PrimMotor'; ...
%     'Left_Parieto_Occipital', 'L_ParOcc';    ...
%     'Right_Parieto_Occipital','R_ParOcc'     };
eeg_clusters = {'L_PrimMotor'};
cluster_defs = { ...
    'Left_Prim_Motor',        'L_PrimMotor'};
n_clusters = numel(eeg_clusters);

%% ── TF and CMC parameters ────────────────────────────────────────────────
emg_labels     = {'VL', 'RF', 'GM', 'BF'};
emg_srate      = 2000;
eeg_srate      = 500;
epoch_limits_s = [-0.5, 3.5];   % original FlxS-locked epoch window

tf_params.srate     = eeg_srate;
tf_params.freqs     = [3 60];
tf_params.nfreqs    = 60;
tf_params.cycles    = [2 0.8];
tf_params.freqscale = 'log';
tf_params.ntimesout = 200;
tf_params.verbose   = false;

cmc_opts.plot          = true;
cmc_opts.colormap      = 'turbo';
cmc_opts.muscle_labels = emg_labels;

%% ── Re-epoch event definitions ───────────────────────────────────────────
% Column 1: event label
% Column 2: column index in EEG_epoched_main.timewarp.latencies
%           (1 = FlxS at 0 ms, 2 = ExtS, 3 = ExtE)
% Column 3: epoch window [tmin tmax] in seconds
%
% All events use a symmetric [-1, 1] s window because we extract from the
% continuous EEG, which has no pre-epoch buffer constraint.
reepoch_defs = { ...
    'FlxS',  1,  [-1.0  1.0]; ...
    'ExtS',  2,  [-1.0  1.0]; ...
    'ExtE',  3,  [-1.0  1.0]; ...
};
n_reepoch = size(reepoch_defs, 1);

%% ── Load IC cluster assignments ──────────────────────────────────────────
load('Subjects_ICs_in_clusters.mat')   % -> SUBJECTS_ICS

%% ── Group-level storage ──────────────────────────────────────────────────
all_coh             = cell(n_subjects, n_clusters, n_proc, n_reepoch);
freqs_ref           = [];
reepoch_time_ax_ref = cell(n_reepoch, 1);

%% ══════════════════════════════════════════════════════════════════════════
%  SUBJECT LOOP
% ══════════════════════════════════════════════════════════════════════════
for si = 1:n_subjects

    subject_id = subject_ids(si);
    nChan      = numel(emg_ids);
    fprintf('\n════════  Subject %d  (%d / %d)  ════════\n', ...
        subject_id, si, n_subjects);

    %% ── Load EEG datasets ────────────────────────────────────────────────
    try
        % Raw merged EEG — needed for import_same_event_on_EMG_stream
        cd([rawEEGLAB_path, 'sub-', num2str(subject_id)])
        EEG_orig = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_merged_EEG.set']);
        cd(this_path)

        % Epoched EEG — needed for timewarp latencies and urevent lookup
        cd(epoched_EEG_path)
        EEG_epoched_main = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_cleaned_with_ICA_epoched.set']);
        cd(this_path)

        % Continuous ICA EEG — used to extract EEG epochs with no window constraint
        cd(cont_ICA_path)
        EEG_cont = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_cleaned_with_ICA.set']);
        cd(this_path)

    catch ME
        warning('Subject %d: EEG load failed — %s', subject_id, ME.message);
        continue
    end

    tw_latencies = EEG_epoched_main.timewarp.latencies; % [nEpochs x nEvents] ms
    N_EEG        = length(EEG_epoched_main.epoch);

    %% ── Load raw EMG ─────────────────────────────────────────────────────
    try
        output  = runs_concatenated(subject_id, rawdata_path);
        EMG_raw = double(output.All_EMG(emg_ids, :));
        clear output
    catch ME
        warning('Subject %d: EMG load failed — %s', subject_id, ME.message);
        continue
    end

    %% ── Identify clusters present and collect needed IC indices ──────────
    cluster_ic_ids  = NaN(1, n_clusters);   % IC index per cluster slot
    cluster_present = false(1, n_clusters);
    ic_ids_needed   = [];

    for c = 1:size(cluster_defs, 1)
        cname_ic = cluster_defs{c, 1};
        cfield   = cluster_defs{c, 2};
        ci_slot  = find(strcmp(eeg_clusters, cfield));
        cidx     = find(strcmpi(SUBJECTS_ICS(:, 1), cname_ic));
        sidx     = find(SUBJECTS_ICS{cidx, 2}.Subjects == subject_id - 4);
        if isempty(sidx)
            fprintf('  Cluster %-25s not present.\n', cname_ic);
            continue
        end
        ic_id                   = SUBJECTS_ICS{cidx, 2}.ICs(sidx);
        cluster_ic_ids(ci_slot) = ic_id;
        cluster_present(ci_slot)= true;
        ic_ids_needed           = [ic_ids_needed, ic_id]; %#ok<AGROW>
    end
    ic_ids_needed = unique(ic_ids_needed);

    if ~any(cluster_present)
        warning('Subject %d: not in any cluster — skipping.', subject_id);
        continue
    end
    fprintf('  Clusters: %s\n', strjoin(eeg_clusters(cluster_present), ', '));

    %% ── Compute ICA activations for needed ICs only ──────────────────────
    % icaact is NOT stored in the .set file — compute on the fly.
    % Only the rows (ICs) needed for our clusters are computed,
    % keeping memory use proportional to n_clusters rather than n_ICs_total.
    fprintf('  Computing ICA activations for %d IC(s)...\n', numel(ic_ids_needed));
    W           = EEG_cont.icaweights(ic_ids_needed, :) * EEG_cont.icasphere;
    ic_cont_all = W * double(EEG_cont.data);   % [n_needed x nSamples_cont]
    n_cont_smp  = size(ic_cont_all, 2);
    clear EEG_cont W   % free memory — continuous data is large

    %% ── Pre-compute continuous FlxS sample indices from urevent ──────────
    % EEG_epoched_main.urevent stores the original sample indices (in the
    % continuous dataset sub-{id}_cleaned_with_ICA.set) for every event.
    % We use these to locate FlxS in the continuous ICA activation vector,
    % then offset by tw_latencies to find ExtS and ExtE.
    event_epoch_vec = [EEG_epoched_main.event.epoch];
    event_type_vec  = {EEG_epoched_main.event.type};
    FlxS_cont_sample = NaN(N_EEG, 1);

    for i = 1:N_EEG
        ev_in_ep  = find(event_epoch_vec == i);
        flxs_mask = strcmpi(event_type_vec(ev_in_ep), 'FlxS');
        flxs_evs  = ev_in_ep(flxs_mask);
        if isempty(flxs_evs), continue; end
        ur_idx = EEG_epoched_main.event(flxs_evs(1)).urevent;
        if isempty(ur_idx) || ur_idx > length(EEG_epoched_main.urevent)
            continue
        end
        FlxS_cont_sample(i) = round(EEG_epoched_main.urevent(ur_idx).latency);
    end

    %% ── Epoch timing alignment for EMG (once per subject) ────────────────
    [b_tmp, a_tmp] = butter(4, [20 450] / (emg_srate/2), 'bandpass');
    EMG_filt_tmp   = filtfilt(b_tmp, a_tmp, double(EMG_raw'))';

    eeg_emptyset       = structfun(@(x) [], EEG_orig, 'UniformOutput', false);
    EMG_tmpl           = eeg_emptyset;
    EMG_tmpl.data      = double(EMG_filt_tmp);
    EMG_tmpl.nbchan    = nChan;
    EMG_tmpl.pnts      = size(EMG_filt_tmp, 2);
    EMG_tmpl.trials    = 1;
    EMG_tmpl.srate     = emg_srate;
    EMG_tmpl.xmin      = 0;
    EMG_tmpl.xmax      = (EMG_tmpl.pnts - 1) / emg_srate;
    EMG_tmpl.setname   = 'EMG';
    EMG_tmpl.filename  = '';
    EMG_tmpl.filepath  = '';
    for ch = 1:nChan
        EMG_tmpl.chanlocs(ch).labels = emg_labels{ch};
    end
    EMG_tmpl = eeg_checkset(EMG_tmpl);

    output_timing = runs_concatenated(subject_id, rawdata_path);
    [EMG_ev, ~]   = import_same_event_on_EMG_stream( ...
        EMG_tmpl, EEG_orig, output_timing, ...
        subject_id, processing_path, data_path);
    clear EMG_tmpl EMG_filt_tmp output_timing

    for e = 1:length(EMG_ev.event)
        nums = str2double(strsplit(EMG_ev.event(e).desc, '_')).';
        if strcmp(EMG_ev.event(e).type, 'boundary')
            EMG_ev.event(e).trial = 'none';
            EMG_ev.event(e).cond  = 'none';
        else
            EMG_ev.event(e).trial = nums(end);
            EMG_ev.event(e).cond  = nums(2);
        end
    end

    EMG_ev_sel       = pop_selectevent(EMG_ev, ...
        'type', {'FlxS','ExtS','ExtE'}, 'deleteevents', 'on');
    EMG_epoched_meta = pop_epoch(EMG_ev_sel, {'FlxS'}, epoch_limits_s, ...
        'newname', 'meta', 'epochinfo', 'yes');
    EMG_epoched_meta = eeg_checkset(EMG_epoched_meta);
    clear EMG_ev EMG_ev_sel

    [matched_vector, ~] = find_matched_epochs(EEG_epoched_main, EMG_epoched_meta);

    EMG_time_meta   = EMG_epoched_meta.times;
    N_EMG_smp       = EMG_epoched_meta.pnts;
    delta_t_ms      = NaN(length(matched_vector), 1);
    FlxS_abs_time_s = NaN(size(EMG_epoched_meta.epoch));

    for i = 1:length(matched_vector)
        if isnan(matched_vector(i)), continue; end
        ep       = EMG_epoched_meta.epoch(matched_vector(i));
        FlxS_idx = find(strcmpi(ep.eventtype,'FlxS') & ...
            cell2mat(ep.eventlatency) == 0, 1);
        if isempty(FlxS_idx), continue; end
        t_last  = ep.eventinit_time_lastEMGpoint{FlxS_idx};
        t_FlxS  = ep.eventinit_time{FlxS_idx};
        delta_t_ms(i)                      = 1000 * (t_FlxS - t_last);
        FlxS_abs_time_s(matched_vector(i)) = t_FlxS;
    end
    clear EMG_epoched_meta

    % Valid epochs: matched EMG, valid delta_t, and FlxS found in continuous EEG
    keep           = ~isnan(matched_vector) & ...
                     ~isnan(delta_t_ms)     & ...
                     ~isnan(FlxS_cont_sample);
    valid_orig_idx = find(keep);
    fprintf('  Valid epochs: %d / %d\n', sum(keep), N_EEG);

    %% ── Bandpass filter (shared by all processing conditions) ────────────
    [b_filt, a_filt] = butter(4, [20 450] / (emg_srate/2), 'bandpass');
    EMG_filt         = filtfilt(b_filt, a_filt, double(EMG_raw'))';
    clear EMG_raw

    %% ── Processing condition loop ─────────────────────────────────────────
    for pi = 1:n_proc
        proc_label = proc_conditions(pi).label;
        fprintf('\n  Processing: %s\n', proc_conditions(pi).desc);

        z = hilbert(double(EMG_filt'))';
        switch proc_label
            case 'hilbert_amp',   EMG_proc = abs(z);
            case 'hilbert_demod', EMG_proc = real(z ./ abs(z));
        end
        clear z

        %% ── Re-epoch event loop ──────────────────────────────────────────
        for ev = 1:n_reepoch
            ev_name  = reepoch_defs{ev, 1};
            ev_col   = reepoch_defs{ev, 2};
            win_s    = reepoch_defs{ev, 3};
            win_ms   = win_s * 1000;
            n_smp_re = round(diff(win_s) * eeg_srate) + 1;
            t_re_ms  = linspace(win_ms(1), win_ms(2), n_smp_re);

            % Sample offsets from FlxS for EEG extraction
            eeg_offset_s1 = round(win_s(1) * eeg_srate);
            eeg_offset_s2 = round(win_s(2) * eeg_srate);

            fprintf('\n    Event: %s  [%g %g] s\n', ev_name, win_s(1), win_s(2));

            n_valid = length(valid_orig_idx);

            %% ── Extract EMG re-epochs (cluster-independent) ──────────────
            EMG_reepoch = NaN(nChan, n_smp_re, n_valid);

            for k = 1:n_valid
                orig_idx = valid_orig_idx(k);
                mi       = matched_vector(orig_idx);

                ev_lat_ms = tw_latencies(orig_idx, ev_col);
                if isnan(ev_lat_ms), continue; end

                ev_abs_s = FlxS_abs_time_s(mi) + ev_lat_ms / 1000;
                
                % Add a small margin on each side to absorb the delta_t shift
                % (delta_t_ms < 0.5 ms at 2000 Hz, so 5 ms is more than enough)
                margin_s = 0.005;
                
                s1 = round((ev_abs_s + win_s(1) - margin_s) * emg_srate) + 1;
                s2 = round((ev_abs_s + win_s(2) + margin_s) * emg_srate);
                if s1 < 1 || s2 > size(EMG_proc, 2), continue; end
                
                emg_win = EMG_proc(:, s1:s2);
                
                % Time axis reflects the wider window, including the margin
                emg_t = linspace((win_ms(1) - margin_s*1000), ...
                                 (win_ms(2) + margin_s*1000), size(emg_win, 2));
                emg_t = emg_t - delta_t_ms(orig_idx);
                
                EMG_reepoch(:, :, k) = interp1(emg_t, emg_win', t_re_ms, 'linear', NaN)';

            end

            nan_re    = squeeze(any(any(isnan(EMG_reepoch), 1), 2));
            keep_re   = ~nan_re;
            EMG_re_ok = EMG_reepoch(:, :, keep_re);
            orig_ok   = valid_orig_idx(keep_re);
            clear EMG_reepoch

            if sum(keep_re) < 10
                fprintf('    Too few valid EMG epochs (%d) — skipping.\n', ...
                    sum(keep_re));
                continue
            end
            fprintf('    Valid EMG epochs: %d\n', sum(keep_re));

            %% ── EEG cluster loop ─────────────────────────────────────────
            for ci = 1:n_clusters
                cname = eeg_clusters{ci};
                if ~cluster_present(ci), continue; end

                fig_path = fullfile(results_path, ...
                    sprintf('sub%02d__%s__reepoch_%s__%s.png', ...
                    subject_id, proc_label, ev_name, cname));

                if exist(fig_path, 'file')
                    fprintf('    %s | %s — exists, skipping.\n', ev_name, cname);
                    mat_path = fullfile(results_path, ...
                        sprintf('sub%02d__coh_reepoch.mat', subject_id));
                    if exist(mat_path, 'file')
                        tmp = load(mat_path, 'coh_sub');
                        all_coh{si, ci, pi, ev} = tmp.coh_sub{ci, pi, ev};
                    end
                    continue
                end

                %% ── Extract EEG re-epochs from continuous ICA activation ──
                % For each valid epoch, find the event's sample in the
                % continuous EEG and extract the window directly.
                ic_id  = cluster_ic_ids(ci);
                ic_row = find(ic_ids_needed == ic_id);

                X_reepoch = NaN(n_smp_re, length(orig_ok));

                for k = 1:length(orig_ok)
                    orig_idx  = orig_ok(k);
                    ev_lat_ms = tw_latencies(orig_idx, ev_col);
                    if isnan(ev_lat_ms), continue; end

                    % Continuous sample of this event:
                    %   FlxS sample + offset for ExtS/ExtE (0 for FlxS itself)
                    ev_cont_smp = FlxS_cont_sample(orig_idx) + ...
                        round(ev_lat_ms * eeg_srate / 1000);

                    s1 = ev_cont_smp + eeg_offset_s1;
                    s2 = ev_cont_smp + eeg_offset_s2;
                    if s1 < 1 || s2 > n_cont_smp, continue; end

                    X_reepoch(:, k) = ic_cont_all(ic_row, s1:s2)';
                end

                % Drop epochs where EEG extraction also failed
                nan_eeg = any(isnan(X_reepoch), 1);
                X_valid = X_reepoch(:, ~nan_eeg);
                Y_valid = EMG_re_ok(:, :, ~nan_eeg);

                if size(X_valid, 2) < 10
                    warning('Subject %d | %s | %s | %s: only %d valid epochs.', ...
                        subject_id, proc_label, ev_name, cname, size(X_valid,2));
                    continue
                end

                %% ── Compute TF and CMC ───────────────────────────────────
                try
                    [TF_X, TF_Y, time_axis, freqs] = compute_TF_locked( ...
                        X_valid, Y_valid, tf_params, win_ms);

                    if isempty(reepoch_time_ax_ref{ev})
                        reepoch_time_ax_ref{ev} = time_axis;
                    end
                    if isempty(freqs_ref)
                        freqs_ref = freqs;
                    end

                    cmc_opts.title        = sprintf('sub%02d  %s  %s  %s  (%d epochs)', ...
                        subject_id, strrep(proc_label,'_','\_'), ...
                        ev_name, strrep(cname,'_',' '), size(X_valid, 2));
                    cmc_opts.event_label  = ev_name;
                    
                    [coh, fig_h] = compute_and_plot_CMC_eventlocked( ...
                        TF_X, TF_Y, time_axis, freqs, cmc_opts);
                    clear TF_X TF_Y

                    all_coh{si, ci, pi, ev} = coh;

                    if ~isempty(fig_h) && isvalid(fig_h)
                        try
                            exportgraphics(fig_h, fig_path, 'Resolution', 150);
                        catch
                            saveas(fig_h, fig_path, 'png');
                        end
                        close(fig_h);
                    end
                    fprintf('    %s | %s — done (%d epochs).\n', ...
                        ev_name, cname, size(X_valid, 2));

                catch ME
                    warning('Subject %d | %s | %s | %s: %s', ...
                        subject_id, proc_label, ev_name, cname, ME.message);
                end

            end % cluster loop

            clear EMG_re_ok

        end % event loop

        clear EMG_proc

    end % processing loop

    clear EMG_filt ic_cont_all

    %% ── Save per-subject CMC matrices ─────────────────────────────────────
    coh_sub  = squeeze(all_coh(si, :, :, :));
    mat_path = fullfile(results_path, ...
        sprintf('sub%02d__coh_reepoch.mat', subject_id));
    save(mat_path, 'coh_sub', 'freqs_ref', 'reepoch_time_ax_ref', ...
        'subject_id', '-v7.3');
    fprintf('\n  Subject %d — saved.\n', subject_id);

end % subject loop

fprintf('\nAll subjects done.\n');

%% ══════════════════════════════════════════════════════════════════════════
%  RELOAD saved per-subject matrices (handles interrupted runs)
% ══════════════════════════════════════════════════════════════════════════
fprintf('\nReloading saved per-subject matrices...\n');
for si = 1:n_subjects
    subject_id = subject_ids(si);
    mat_path   = fullfile(results_path, ...
        sprintf('sub%02d__coh_reepoch.mat', subject_id));
    if ~exist(mat_path, 'file')
        fprintf('  sub%02d — no .mat file found.\n', subject_id);
        continue
    end
    tmp = load(mat_path, 'coh_sub', 'freqs_ref', 'reepoch_time_ax_ref');
    all_coh(si, :, :, :) = tmp.coh_sub;
    if isempty(freqs_ref) && ~isempty(tmp.freqs_ref)
        freqs_ref           = tmp.freqs_ref;
        reepoch_time_ax_ref = tmp.reepoch_time_ax_ref;
    end
    fprintf('  sub%02d — loaded.\n', subject_id);
end

if isempty(freqs_ref)
    error('No valid .mat files found — cannot compute group averages.');
end

%% ══════════════════════════════════════════════════════════════════════════
%  GROUP-LEVEL AVERAGE CMC
% ══════════════════════════════════════════════════════════════════════════
fprintf('\nComputing group-average CMC figures...\n');

for ev = 1:n_reepoch
    ev_name   = reepoch_defs{ev, 1};
    time_axis = reepoch_time_ax_ref{ev};
    if isempty(time_axis)
        warning('No time axis for event %s — skipping.', ev_name); continue
    end

    for pi = 1:n_proc
        proc_label = proc_conditions(pi).label;

        for ci = 1:n_clusters
            cname = eeg_clusters{ci};

            coh_stack  = [];
            valid_subs = [];
            for si = 1:n_subjects
                if isempty(all_coh{si, ci, pi, ev}), continue; end
                coh_stack  = cat(4, coh_stack, all_coh{si, ci, pi, ev});
                valid_subs = [valid_subs, subject_ids(si)]; %#ok<AGROW>
            end

            if isempty(coh_stack)
                warning('No valid subjects: %s | %s | %s.', ...
                    ev_name, proc_label, cname); continue
            end

            n_valid  = numel(valid_subs);
            coh_mean = mean(coh_stack, 4, 'omitnan');
            coh_sem  = std(coh_stack, 0, 4, 'omitnan') / sqrt(n_valid);

            fprintf('  %s | %s | %s: n = %d  (IDs: %s)\n', ...
                ev_name, proc_label, cname, n_valid, num2str(valid_subs));

            save(fullfile(results_path, ...
                sprintf('group__coh__reepoch_%s__%s__%s.mat', ...
                ev_name, proc_label, cname)), ...
                'coh_mean', 'coh_sem', 'time_axis', 'freqs_ref', ...
                'valid_subs', '-v7.3');

            group_title = sprintf('%s  |  %s  |  %s   (n = %d)', ...
                ev_name, strrep(proc_label,'_','\_'), ...
                strrep(cname,'_',' '), n_valid);
            fig_h = plot_group_CMC(coh_mean, time_axis, freqs_ref, ...
                emg_labels, group_title, ev_name);
            fig_path = fullfile(results_path, ...
                sprintf('group__CMC__reepoch_%s__%s__%s.png', ...
                ev_name, proc_label, cname));
            try; exportgraphics(fig_h, fig_path, 'Resolution', 150);
            catch; saveas(fig_h, fig_path, 'png'); end
            close(fig_h);
            fprintf('  Saved: %s\n', fig_path);

        end
    end
end

fprintf('\nFinished.\n');


%% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function [matched_vector, condition_vector] = find_matched_epochs( ...
        EEG_epoched_main, EMG_epoched)
    eventtypes_EMG   = {EMG_epoched.epoch.eventtype};
    eventlatency_EMG = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventlatency}, 'UniformOutput', false);
    FlxS_type_lgl = cellfun(@(x) strcmpi(x,'FlxS'), eventtypes_EMG, ...
        'UniformOutput', false);
    FlxS_lat_lgl  = cellfun(@(x) x == 0, eventlatency_EMG, ...
        'UniformOutput', false);
    FlxS_idx_EMG  = cellfun(@(x,y) find(x & y), FlxS_type_lgl, ...
        FlxS_lat_lgl, 'UniformOutput', false);
    more_one = cellfun(@(x) numel(x) > 1, FlxS_idx_EMG);

    init_time_cell = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventinit_time_EEG}, 'UniformOutput', false);
    init_time_EMG  = cellfun(@(x,y) x(y(1)), init_time_cell, FlxS_idx_EMG, ...
        'UniformOutput', false);
    init_time_EMG  = cell2mat(init_time_EMG);

    N_EEG = length(EEG_epoched_main.epoch);
    matched_vector   = NaN(N_EEG, 1);
    condition_vector = zeros(N_EEG, 1);
    for i = 1:N_EEG
        evtypes = EEG_epoched_main.epoch(i).eventtype;
        evlat   = cell2mat(EEG_epoched_main.epoch(i).eventlatency);
        FlxS_i  = find(strcmpi(evtypes,'FlxS') & evlat == 0, 1);
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


function fig_h = plot_group_CMC(coh_mean, time_axis, freqs, ...
        muscle_labels, fig_title, ev_name)
    nMuscles   = size(coh_mean, 3);
    nRows      = ceil(sqrt(nMuscles));
    nCols      = ceil(nMuscles / nRows);
    cmax       = prctile(coh_mean(:), 95);
    if isnan(cmax) || cmax == 0, cmax = 1; end
    freq_ticks = [4 8 14 30 60];

    fig_h = figure('Color', 'w', 'Name', fig_title);
    for m = 1:nMuscles
        ax = subplot(nRows, nCols, m);
        contourf(time_axis, freqs, coh_mean(:,:,m), 200, 'LineColor','none');
        set(ax, 'yscale','log', 'ydir','norm', ...
            'ylim', [freqs(1) freqs(end)] + [0 0.01], ...
            'ytick', freq_ticks, 'box', 'on', 'FontSize', 12);
        colormap(ax, 'turbo');
        caxis([0 cmax]);
        hold on
        % xline(0, 'w--', 'LineWidth', 2, 'Label', ev_name, ...
        %     'LabelVerticalAlignment', 'bottom', 'FontSize', 10);
        xline(0, 'k--', 'LineWidth', 2);
        hold off
        title(muscle_labels{m}, 'FontWeight', 'bold');
        xlabel('Time relative to event (ms)');
        ylabel('Frequency (Hz)');
    end
    cb = colorbar(subplot(nRows, nCols, nMuscles), ...
        'Position', [0.93 0.13 0.015 0.755]);
    ylabel(cb, 'MSC');
    sgtitle(fig_title, 'FontSize', 15, 'FontWeight', 'normal');
end