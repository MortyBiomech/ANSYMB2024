%% run_CMC_allsubjects.m
%
% Cortico-muscular coherence — all subjects, two EMG processing conditions:
%   EMG filter    : band-pass 20-450 Hz
%   EMG processing: (1) Hilbert amplitude  |H{x}| = A(t)
%                   (2) Hilbert demodulation  cos(phi(t))
%   TF method     : upsample to max SOI, cycle-% axis (compute_TF_timenorm)
%   EEG clusters  : L/R primary motor, L/R parieto-occipital
%
% Outputs per subject : one PNG figure per processing condition x EEG cluster
% Output at group level: one PNG per processing condition x cluster (mean MSC)

clc; clear;

%% ── Subject list ─────────────────────────────────────────────────────────
subject_ids = [5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18];
emg_ids     = [2, 3, 4, 5];   % same for all subjects
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
cmc_folder_path  = [processing_path, 'Group_Level_PostProcessing\', ...
                    'Final_paper_plot_generation\_NHB\', ...
                    'Corticomuscular_Coherence\run_all_cases\'];
results_path     = [cmc_folder_path, 'CMC_results\bp20450_hilbert_upsample\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

%% ── Initialise EEGLAB ────────────────────────────────────────────────────
this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

%% ── Processing conditions ────────────────────────────────────────────────
proc_conditions = struct( ...
    'label', {'hilbert_amp',              'hilbert_demod'       }, ...
    'desc',  {'Hilbert amplitude |H{x}|', 'Hilbert demod. cos(phi)'});
n_proc = numel(proc_conditions);

%% ── EEG cluster definitions ──────────────────────────────────────────────
eeg_clusters = {'L_PrimMotor', 'R_PrimMotor', 'L_ParOcc', 'R_ParOcc'};
cluster_defs = { ...
    'Left_Prim_Motor',        'L_PrimMotor'; ...
    'Right_Prim_Motor',       'R_PrimMotor'; ...
    'Left_Parieto_Occipital', 'L_ParOcc';    ...
    'Right_Parieto_Occipital','R_ParOcc'     };
n_clusters = numel(eeg_clusters);

%% ── TF and CMC parameters ────────────────────────────────────────────────
emg_labels     = {'VL', 'RF', 'GM', 'BF'};
emg_srate      = 2000;
eeg_srate      = 500;
epoch_limits_s = [-0.5, 3.5];

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

%% ── Load cluster IC assignments ──────────────────────────────────────────
load('Subjects_ICs_in_clusters.mat')   % -> SUBJECTS_ICS

%% ── Group-level storage ──────────────────────────────────────────────────
% all_coh{si, ci, pi} = [nFreqs x nTimes x nMuscles] or []
all_coh          = cell(n_subjects, n_clusters, n_proc);
percent_grid_ref = [];
freqs_ref        = [];

%% ══════════════════════════════════════════════════════════════════════════
%  SUBJECT LOOP
% ══════════════════════════════════════════════════════════════════════════
for si = 10

    subject_id = subject_ids(si);
    nChan      = numel(emg_ids);

    fprintf('\n════════  Subject %d  (%d / %d)  ════════\n', ...
        subject_id, si, n_subjects);

    %% ── Load EEG ─────────────────────────────────────────────────────────
    try
        cd([rawEEGLAB_path, 'sub-', num2str(subject_id)])
        EEG_orig = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_merged_EEG.set']);
        cd(this_path)

        cd(epoched_EEG_path)
        EEG_epoched_main = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_cleaned_with_ICA_epoched.set']);
        cd(this_path)
    catch ME
        warning('Subject %d: EEG load failed — %s', subject_id, ME.message);
        continue
    end

    EEG_times    = EEG_epoched_main.times;
    tlimits_tw   = [EEG_epoched_main.xmin EEG_epoched_main.xmax] * 1000;
    tw_latencies = EEG_epoched_main.timewarp.latencies;

    %% ── Load raw EMG ─────────────────────────────────────────────────────
    try
        output  = runs_concatenated(subject_id, rawdata_path);
        EMG_raw = double(output.All_EMG(emg_ids, :));
        clear output
    catch ME
        warning('Subject %d: EMG load failed — %s', subject_id, ME.message);
        continue
    end

    %% ── Extract IC activations (skip missing clusters gracefully) ─────────
    EEG_IC          = struct();
    cluster_present = false(1, n_clusters);

    for c = 1:size(cluster_defs, 1)
        cname_ic = cluster_defs{c, 1};
        cfield   = cluster_defs{c, 2};
        cidx     = find(strcmpi(SUBJECTS_ICS(:, 1), cname_ic));
        sidx     = find(SUBJECTS_ICS{cidx, 2}.Subjects == subject_id - 4);
        if isempty(sidx)
            fprintf('  Cluster %-25s not present for this subject.\n', cname_ic);
            continue
        end
        ic_id              = SUBJECTS_ICS{cidx, 2}.ICs(sidx);
        EEG_IC.(cfield)    = squeeze(EEG_epoched_main.icaact(ic_id, :, :));
        cluster_present(c) = true;
    end

    if ~any(cluster_present)
        warning('Subject %d: not present in any cluster — skipping.', subject_id);
        continue
    end
    fprintf('  Clusters available: %s\n', strjoin(eeg_clusters(cluster_present), ', '));

    %% ── Epoch timing alignment (once per subject) ─────────────────────────
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

    EMG_time_meta = EMG_epoched_meta.times;
    N_EMG_smp     = EMG_epoched_meta.pnts;
    delta_t_ms    = NaN(length(matched_vector), 1);
    epoch_start_s = NaN(size(EMG_epoched_meta.epoch));

    for i = 1:length(matched_vector)
        if isnan(matched_vector(i)), continue; end
        ep       = EMG_epoched_meta.epoch(matched_vector(i));
        FlxS_idx = find(strcmpi(ep.eventtype,'FlxS') & ...
            cell2mat(ep.eventlatency) == 0, 1);
        if isempty(FlxS_idx), continue; end
        t_last = ep.eventinit_time_lastEMGpoint{FlxS_idx};
        t_FlxS = ep.eventinit_time{FlxS_idx};
        delta_t_ms(i)    = 1000 * (t_FlxS - t_last);
        epoch_start_s(matched_vector(i)) = t_FlxS + epoch_limits_s(1);
    end
    clear EMG_epoched_meta

    epoch_start_smp = round(epoch_start_s * emg_srate) + 1;
    epoch_end_smp   = epoch_start_smp + N_EMG_smp - 1;

    timewarp_indices = round(eeg_lat2point( ...
        reshape(tw_latencies, 1, []), 1, eeg_srate, tlimits_tw, 1e-3));
    timewarp_indices = reshape(timewarp_indices, size(tw_latencies));

    %% ── Bandpass filter (shared by both processing conditions) ───────────
    [b_filt, a_filt] = butter(4, [20 450] / (emg_srate/2), 'bandpass');
    EMG_filt         = filtfilt(b_filt, a_filt, double(EMG_raw'))';
    clear EMG_raw

    %% ── Processing condition loop ─────────────────────────────────────────
    for pi = 1:n_proc
        proc_label = proc_conditions(pi).label;
        fprintf('\n  [%s]\n', proc_conditions(pi).desc);

        %% ── Apply processing ─────────────────────────────────────────────
        z = hilbert(double(EMG_filt'))';
        switch proc_label
            case 'hilbert_amp'
                EMG_proc = abs(z);
            case 'hilbert_demod'
                EMG_proc = real(z ./ abs(z));
        end
        clear z

        %% ── Interpolate epochs onto EEG time grid ────────────────────────
        N_EEG     = length(EEG_epoched_main.epoch);
        EMG_intrp = NaN(nChan, length(EEG_times), N_EEG);

        for i = 1:N_EEG
            mi = matched_vector(i);
            if isnan(mi) || isnan(delta_t_ms(i)), continue; end
            s1 = epoch_start_smp(mi);
            s2 = epoch_end_smp(mi);
            if s1 < 1 || s2 > size(EMG_proc, 2), continue; end
            emg_ep = EMG_proc(:, s1:s2);
            emg_t  = EMG_time_meta - delta_t_ms(i);
            EMG_intrp(:, :, i) = ...
                interp1(emg_t', emg_ep', EEG_times', 'linear')';
        end
        clear EMG_proc

        nan_eps        = squeeze(any(any(isnan(EMG_intrp), 1), 2));
        keep           = ~nan_eps;
        EMG_epochs_CMC = EMG_intrp(:, :, keep);
        clear EMG_intrp
        tw_idx = timewarp_indices(keep, :);
        fprintf('  Valid epochs: %d / %d\n', sum(keep), N_EEG);

        %% ── CMC per EEG cluster ──────────────────────────────────────────
        for ci = 1:n_clusters
            cname = eeg_clusters{ci};

            if ~isfield(EEG_IC, cname)
                fprintf('  %s | %s — no IC, skipping.\n', proc_label, cname);
                continue
            end

            fig_path = fullfile(results_path, ...
                sprintf('sub%02d__%s__%s.png', subject_id, proc_label, cname));

            % if exist(fig_path, 'file')
            %     fprintf('  %s | %s — already saved, skipping.\n', proc_label, cname);
            %     mat_path = fullfile(results_path, ...
            %         sprintf('sub%02d__coh.mat', subject_id));
            %     if exist(mat_path, 'file')
            %         tmp = load(mat_path, 'coh_sub');
            %         all_coh{si, ci, pi} = tmp.coh_sub{ci, pi};
            %     end
            %     continue
            % end

            X_epochs = EEG_IC.(cname)(:, keep);

            try
                [TF_X, TF_Y, percent_grid, freqs] = compute_TF_timenorm( ...
                    X_epochs, EMG_epochs_CMC, tw_idx, tf_params);

                if isempty(percent_grid_ref)
                    percent_grid_ref = percent_grid;
                    freqs_ref        = freqs;
                end

                cmc_opts.title = sprintf('sub%02d  %s  %s', subject_id, ...
                    strrep(proc_label,'_','\_'), strrep(cname,'_',' '));

                [coh, fig_h] = compute_and_plot_CMC( ...
                    TF_X, TF_Y, percent_grid, freqs, cmc_opts);
                clear TF_X TF_Y

                all_coh{si, ci, pi} = coh;

                if ~isempty(fig_h) && isvalid(fig_h)
                    try
                        exportgraphics(fig_h, fig_path, 'Resolution', 150);
                    catch
                        saveas(fig_h, fig_path, 'png');
                    end
                    close(fig_h);
                end
                fprintf(' sub-%d: %s | %s — done.\n', subject_id ,proc_label, cname);

            catch ME
                warning('Subject %d | %s | %s: %s', ...
                    subject_id, proc_label, cname, ME.message);
            end

        end % cluster loop

        clear EMG_epochs_CMC

    end % processing loop

    clear EMG_filt

    %% ── Save per-subject CMC matrices ─────────────────────────────────────
    coh_sub  = squeeze(all_coh(si, :, :));   % {n_clusters x n_proc}
    mat_path = fullfile(results_path, sprintf('sub%02d__coh.mat', subject_id));
    save(mat_path, 'coh_sub', 'percent_grid_ref', 'freqs_ref', ...
        'subject_id', '-v7.3');
    fprintf('\n  Subject %d — saved.\n', subject_id);

end % subject loop

fprintf('\nAll subjects done. Computing group averages...\n');

%% ══════════════════════════════════════════════════════════════════════════
%  GROUP-LEVEL AVERAGE CMC
% ══════════════════════════════════════════════════════════════════════════

%% ── Reload saved per-subject CMC matrices ────────────────────────────────
fprintf('\nLoading saved per-subject CMC matrices...\n');
for si = 1:n_subjects
    subject_id = subject_ids(si);
    mat_path   = fullfile(results_path, sprintf('sub%02d__coh.mat', subject_id));

    if ~exist(mat_path, 'file')
        fprintf('  sub%02d — no .mat file found, skipping.\n', subject_id);
        continue
    end

    tmp = load(mat_path, 'coh_sub', 'percent_grid_ref', 'freqs_ref');
    all_coh(si, :, :) = tmp.coh_sub;

    if isempty(percent_grid_ref) && isfield(tmp, 'percent_grid_ref') ...
            && ~isempty(tmp.percent_grid_ref)
        percent_grid_ref = tmp.percent_grid_ref;
        freqs_ref        = tmp.freqs_ref;
    end

    fprintf('  sub%02d — loaded.\n', subject_id);
end

if isempty(percent_grid_ref)
    error('No valid subject .mat files found. Cannot compute group average.');
end

%% ── Compute and save group average ──────────────────────────────────────
fprintf('\nComputing group-average CMC figures...\n');

for pi = 1:n_proc
    proc_label = proc_conditions(pi).label;

    for ci = 1:n_clusters
        cname = eeg_clusters{ci};

        coh_stack  = [];
        valid_subs = [];
        for si = 1:n_subjects
            if isempty(all_coh{si, ci, pi}), continue; end
            coh_stack  = cat(4, coh_stack, all_coh{si, ci, pi});
            valid_subs = [valid_subs, subject_ids(si)]; %#ok<AGROW>
        end

        if isempty(coh_stack)
            warning('No valid subjects for %s | %s.', proc_label, cname);
            continue
        end

        n_valid  = numel(valid_subs);
        coh_mean = mean(coh_stack, 4, 'omitnan');
        coh_sem  = std(coh_stack, 0, 4, 'omitnan') / sqrt(n_valid);

        fprintf('  %s | %s: n = %d subjects  (IDs: %s)\n', ...
            proc_label, cname, n_valid, num2str(valid_subs));

        save(fullfile(results_path, ...
            sprintf('group__coh__%s__%s.mat', proc_label, cname)), ...
            'coh_mean', 'coh_sem', 'percent_grid_ref', 'freqs_ref', ...
            'valid_subs', '-v7.3');

        group_title = sprintf('%s  |  %s   (n = %d)', ...
            strrep(proc_label,'_','\_'), strrep(cname,'_',' '), n_valid);
        fig_h = plot_group_CMC(coh_mean, percent_grid_ref, freqs_ref, ...
            emg_labels, group_title);
        fig_path = fullfile(results_path, ...
            sprintf('group__CMC__%s__%s.png', proc_label, cname));
        try
            exportgraphics(fig_h, fig_path, 'Resolution', 150);
        catch
            saveas(fig_h, fig_path, 'png');
        end
        close(fig_h);
        fprintf('  Saved: %s\n', fig_path);

    end
end

fprintf('\nFinished.\n');


%% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function [matched_vector, condition_vector] = find_matched_epochs( ...
        EEG_epoched_main, EMG_epoched)
% Match EEG and EMG epochs by FlxS event init_time.
    eventtypes_EMG   = {EMG_epoched.epoch.eventtype};
    eventlatency_EMG = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventlatency}, 'UniformOutput', false);
    FlxS_type_lgl = cellfun(@(x) strcmpi(x,'FlxS'), ...
        eventtypes_EMG, 'UniformOutput', false);
    FlxS_lat_lgl  = cellfun(@(x) x == 0, ...
        eventlatency_EMG, 'UniformOutput', false);
    FlxS_idx_EMG  = cellfun(@(x,y) find(x & y), ...
        FlxS_type_lgl, FlxS_lat_lgl, 'UniformOutput', false);
    more_one = cellfun(@(x) numel(x) > 1, FlxS_idx_EMG);

    init_time_cell = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventinit_time_EEG}, 'UniformOutput', false);
    init_time_EMG  = cellfun(@(x,y) x(y(1)), ...
        init_time_cell, FlxS_idx_EMG, 'UniformOutput', false);
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


function fig_h = plot_group_CMC(coh_mean, percent_grid, freqs, ...
        muscle_labels, fig_title)
% 2x2 contourf of mean MSC across subjects.
%   coh_mean : [nFreqs x nTimes x nMuscles]
    nMuscles   = size(coh_mean, 3);
    nRows      = ceil(sqrt(nMuscles));
    nCols      = ceil(nMuscles / nRows);
    cmax       = prctile(coh_mean(:), 95);
    if isnan(cmax) || cmax == 0, cmax = 1; end
    freq_ticks = [4 8 14 30 60];

    fig_h = figure('Color', 'w', 'Name', fig_title);
    for m = 1:nMuscles
        ax = subplot(nRows, nCols, m);
        contourf(percent_grid, freqs, coh_mean(:,:,m), 200, 'LineColor','none');
        set(ax, 'yscale','log', 'ydir','norm', ...
            'ylim', [freqs(1) freqs(end)] + [0 0.01], ...
            'xlim', [-5 105], 'ytick', freq_ticks, ...
            'box', 'on', 'FontSize', 12);
        colormap(ax, 'turbo');
        caxis([0 cmax]);
        hold on
        xline(0,   'w--', 'LineWidth', 1.5);
        xline(100, 'w--', 'LineWidth', 1.5);
        hold off
        title(muscle_labels{m}, 'FontWeight', 'bold');
        xlabel('Cycle (%)');
        ylabel('Frequency (Hz)');
    end
    cb = colorbar(subplot(nRows, nCols, nMuscles), ...
        'Position', [0.93 0.13 0.015 0.755]);
    ylabel(cb, 'MSC');
    sgtitle(fig_title, 'FontSize', 15, 'FontWeight', 'normal');
end