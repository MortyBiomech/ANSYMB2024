%% run_CMC_synergy.m
%
% CMC analysis using group muscle synergies (k=3, across-subject NMF).
% Two TF approaches run per subject:
%
%   Section A — Upsample / cycle-% (compute_TF_timenorm)
%               H_epochs loaded directly from synergy files.
%               No EMG preprocessing needed — fast.
%
%   Section B — Re-epoch locked to FlxS / ExtS / ExtE (compute_TF_locked)
%               Synergy activations computed on-the-fly from continuous EMG
%               using the pseudo-inverse of the shared W:  H = max(0, pinv(W)*EMG)
%               EEG extracted from the continuous ICA dataset.
%               Same ±1 s symmetric window for all three events.
%
% Prerequisites
%   run_synergy_extraction_group.m must have completed successfully (k=3).
%   Shared W and per-subject H_epochs must be saved in synergy_path.

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
cont_ICA_path    = [data_path, '5_single-subject-EEG-analysis\timewarp_test\'];
cmc_folder_path  = [processing_path, 'Group_Level_PostProcessing\', ...
                    'Final_paper_plot_generation\_NHB\', ...
                    'Corticomuscular_Coherence\run_all_cases\'];
synergy_path     = [cmc_folder_path, 'CMC_results\synergy_extraction_group\'];
results_path     = [cmc_folder_path, 'CMC_results\synergy_CMC_k3\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

%% ── Load shared synergy weights W ───────────────────────────────────────
k = 3;
grp = load(fullfile(synergy_path, sprintf('group__synergies_k%d.mat', k)));
W            = grp.W;               % [4 x 3]
syn_labels   = grp.syn_labels;      % {'extension','flexion','co-contraction'}
W_pinv       = pinv(W);             % [3 x 4]  used in Section B

fprintf('Shared W loaded  (k=%d  synergies: %s)\n', k, ...
    strjoin(syn_labels, ' | '));

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
load('Subjects_ICs_in_clusters.mat')

%% ── TF / CMC parameters ─────────────────────────────────────────────────
emg_labels     = {'VL','RF','GM','BF'};
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
cmc_opts.muscle_labels = syn_labels;   % synergy names on plot

%% ── Re-epoch definitions (Section B) ────────────────────────────────────
reepoch_defs = { ...
    'FlxS',  1,  [-1.0  1.0]; ...
    'ExtS',  2,  [-1.0  1.0]; ...
    'ExtE',  3,  [-1.0  1.0]; ...
};
n_reepoch = size(reepoch_defs, 1);
margin_s  = 0.005;   % 5 ms margin to prevent interp1 boundary NaN

%% ── Group storage ────────────────────────────────────────────────────────
all_coh_A       = cell(n_subjects, n_clusters);        % Section A
all_coh_B       = cell(n_subjects, n_clusters, n_reepoch); % Section B
pct_grid_ref    = [];
freqs_ref       = [];
reepoch_tax_ref = cell(n_reepoch, 1);

%% ══════════════════════════════════════════════════════════════════════════
%  SUBJECT LOOP
% ══════════════════════════════════════════════════════════════════════════
for si = 2:n_subjects

    subject_id = subject_ids(si);
    fprintf('\n════════  Subject %d  (%d/%d)  ════════\n', ...
        subject_id, si, n_subjects);

    %% ── Load synergy file ────────────────────────────────────────────────
    syn_file = fullfile(synergy_path, ...
        sprintf('sub%02d__synergies_group_k%d.mat', subject_id, k));
    if ~exist(syn_file, 'file')
        warning('Subject %d: synergy file not found — skipping.', subject_id);
        continue
    end
    syn = load(syn_file, 'H_epochs', 'keep');
    H_epochs = syn.H_epochs;   % [k x n_smp x n_valid]
    keep     = syn.keep;       % [N_EEG x 1] logical

    n_valid = size(H_epochs, 3);
    fprintf('  Valid epochs from synergy file: %d\n', n_valid);

    %% ── Load EEG datasets ────────────────────────────────────────────────
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
    tlimits_tw   = [EEG_epoched_main.xmin EEG_epoched_main.xmax]*1000;
    tw_latencies = EEG_epoched_main.timewarp.latencies;

    %% ── IC activations ───────────────────────────────────────────────────
    EEG_IC          = struct();
    cluster_present = false(1, n_clusters);
    cluster_ic_ids  = NaN(1, n_clusters);
    ic_ids_needed   = [];

    for c = 1:size(cluster_defs,1)
        cname_ic = cluster_defs{c,1};  cfield = cluster_defs{c,2};
        ci_slot  = find(strcmp(eeg_clusters, cfield));
        cidx     = find(strcmpi(SUBJECTS_ICS(:,1), cname_ic));
        sidx     = find(SUBJECTS_ICS{cidx,2}.Subjects == subject_id-4);
        if isempty(sidx)
            fprintf('  Cluster %-25s not present.\n', cname_ic); continue
        end
        ic_id = SUBJECTS_ICS{cidx,2}.ICs(sidx);
        EEG_IC.(cfield)      = squeeze(EEG_epoched_main.icaact(ic_id,:,:));
        cluster_present(ci_slot) = true;
        cluster_ic_ids(ci_slot)  = ic_id;
        ic_ids_needed = [ic_ids_needed, ic_id]; %#ok<AGROW>
    end
    if ~any(cluster_present)
        warning('Subject %d: no clusters found.', subject_id); continue
    end

    %% ── Timewarp indices for Section A ───────────────────────────────────
    tw_indices = round(eeg_lat2point( ...
        reshape(tw_latencies,1,[]), 1, eeg_srate, tlimits_tw, 1e-3));
    tw_indices = reshape(tw_indices, size(tw_latencies));
    tw_idx_A   = tw_indices(keep, :);   % filtered by keep mask

    %% ════════════════════════════════════════════════════════════════════
    %  SECTION A — Upsample / cycle-% CMC
    %  Y_epochs = H_epochs from file (already on EEG time grid)
    %% ════════════════════════════════════════════════════════════════════
    fprintf('\n  [Section A] Upsample CMC\n');

    for ci = 1:n_clusters
        cname = eeg_clusters{ci};
        if ~cluster_present(ci), continue; end

        fig_path = fullfile(results_path, ...
            sprintf('A__sub%02d__%s.png', subject_id, cname));
        if exist(fig_path,'file')
            fprintf('  %s — exists, skipping.\n', cname); continue
        end

        X_epochs = EEG_IC.(cname)(:, keep);   % [n_smp x n_valid]
        Y_epochs = H_epochs;                   % [k x n_smp x n_valid]

        try
            [TF_X, TF_Y, pct_grid, freqs] = compute_TF_timenorm( ...
                X_epochs, Y_epochs, tw_idx_A, tf_params);

            if isempty(pct_grid_ref)
                pct_grid_ref = pct_grid; freqs_ref = freqs;
            end

            cmc_opts.title = sprintf('sub%02d  %s  %s', ...
                subject_id, strrep(cname,'_',' '), '(synergies)');
            [coh, fig_h] = compute_and_plot_CMC( ...
                TF_X, TF_Y, pct_grid, freqs, cmc_opts);
            clear TF_X TF_Y

            all_coh_A{si, ci} = coh;
            if ~isempty(fig_h) && isvalid(fig_h)
                try; exportgraphics(fig_h, fig_path, 'Resolution',150);
                catch; saveas(fig_h, fig_path, 'png'); end
                close(fig_h);
            end
            fprintf('  A | %s — done.\n', cname);
        catch ME
            warning('A | sub%d | %s: %s', subject_id, cname, ME.message);
        end
    end

    %% ════════════════════════════════════════════════════════════════════
    %  SECTION B — Re-epoch CMC
    %  Synergy activations: H = max(0, pinv(W) * EMG_proc_window)
    %  EEG: extracted from continuous ICA dataset (±1 s, no constraint)
    %% ════════════════════════════════════════════════════════════════════
    fprintf('\n  [Section B] Re-epoch CMC\n');

    %% ── Load continuous ICA EEG ──────────────────────────────────────────
    try
        cd(cont_ICA_path)
        EEG_cont = pop_loadset('filename', ...
            ['sub-', num2str(subject_id), '_cleaned_with_ICA.set']);
        cd(this_path)
    catch ME
        warning('Subject %d: continuous EEG load failed — %s', ...
            subject_id, ME.message); goto_next = true;
    end

    % Compute IC activations for needed ICs only
    ic_ids_needed = unique(ic_ids_needed);
    W_ica = EEG_cont.icaweights(ic_ids_needed,:) * EEG_cont.icasphere;
    ic_cont_all = W_ica * double(EEG_cont.data);  % [n_ic x n_cont_smp]
    n_cont_smp  = size(ic_cont_all, 2);
    clear EEG_cont W_ica

    %% ── FlxS sample indices in continuous EEG (via urevent) ─────────────
    N_EEG = length(EEG_epoched_main.epoch);
    ev_epoch_vec = [EEG_epoched_main.event.epoch];
    ev_type_vec  = {EEG_epoched_main.event.type};
    FlxS_cont_smp = NaN(N_EEG, 1);

    for i = 1:N_EEG
        ev_in   = find(ev_epoch_vec == i);
        flxs_ev = ev_in(strcmpi(ev_type_vec(ev_in),'FlxS'));
        if isempty(flxs_ev), continue; end
        ur = EEG_epoched_main.event(flxs_ev(1)).urevent;
        if isempty(ur) || ur > length(EEG_epoched_main.urevent), continue; end
        FlxS_cont_smp(i) = round(EEG_epoched_main.urevent(ur).latency);
    end

    %% ── Load raw EMG and timing alignment ────────────────────────────────
    try
        output_r  = runs_concatenated(subject_id, rawdata_path);
        EMG_raw_r = double(output_r.All_EMG(emg_ids,:));
        clear output_r
    catch ME
        warning('Subject %d: EMG load failed — %s', subject_id, ME.message)
        continue
    end

    [b_tmp,a_tmp] = butter(4,[20 450]/(emg_srate/2),'bandpass');
    EMG_filt_tmp  = filtfilt(b_tmp,a_tmp,double(EMG_raw_r'))';

    eeg_emptyset = structfun(@(x)[], EEG_orig,'UniformOutput',false);
    EMG_tmpl = eeg_emptyset;
    EMG_tmpl.data = double(EMG_filt_tmp); EMG_tmpl.nbchan = numel(emg_ids);
    EMG_tmpl.pnts = size(EMG_filt_tmp,2); EMG_tmpl.trials = 1;
    EMG_tmpl.srate = emg_srate; EMG_tmpl.xmin = 0;
    EMG_tmpl.xmax = (EMG_tmpl.pnts-1)/emg_srate;
    EMG_tmpl.setname='EMG'; EMG_tmpl.filename=''; EMG_tmpl.filepath='';
    for ch = 1:numel(emg_ids)
        EMG_tmpl.chanlocs(ch).labels = emg_labels{ch};
    end
    EMG_tmpl = eeg_checkset(EMG_tmpl);

    out_t = runs_concatenated(subject_id, rawdata_path);
    [EMG_ev,~] = import_same_event_on_EMG_stream( ...
        EMG_tmpl, EEG_orig, out_t, subject_id, processing_path, data_path);
    clear EMG_tmpl EMG_filt_tmp out_t

    for e = 1:length(EMG_ev.event)
        nums = str2double(strsplit(EMG_ev.event(e).desc,'_')).';
        if strcmp(EMG_ev.event(e).type,'boundary')
            EMG_ev.event(e).trial='none'; EMG_ev.event(e).cond='none';
        else
            EMG_ev.event(e).trial=nums(end); EMG_ev.event(e).cond=nums(2);
        end
    end
    EMG_ev_sel = pop_selectevent(EMG_ev,...
        'type',{'FlxS','ExtS','ExtE'},'deleteevents','on');
    EMG_meta = pop_epoch(EMG_ev_sel,{'FlxS'},epoch_limits_s,...
        'newname','meta','epochinfo','yes');
    EMG_meta = eeg_checkset(EMG_meta);
    clear EMG_ev EMG_ev_sel

    [matched_vector,~] = find_matched_epochs(EEG_epoched_main, EMG_meta);
    delta_t_ms      = NaN(length(matched_vector),1);
    FlxS_abs_time_s = NaN(size(EMG_meta.epoch));

    for i = 1:length(matched_vector)
        if isnan(matched_vector(i)), continue; end
        ep = EMG_meta.epoch(matched_vector(i));
        fi = find(strcmpi(ep.eventtype,'FlxS') & ...
            cell2mat(ep.eventlatency)==0, 1);
        if isempty(fi), continue; end
        t_last = ep.eventinit_time_lastEMGpoint{fi};
        t_FlxS = ep.eventinit_time{fi};
        delta_t_ms(i)                      = 1000*(t_FlxS-t_last);
        FlxS_abs_time_s(matched_vector(i)) = t_FlxS;
    end
    clear EMG_meta

    % keep_B: valid for Section B (also requires continuous EEG sample)
    keep_B       = keep & ~isnan(FlxS_cont_smp) & ~isnan(delta_t_ms);
    valid_idx_B  = find(keep_B);

    % Bandpass → Hilbert amplitude → continuous EMG_proc
    [b_f,a_f]  = butter(4,[20 450]/(emg_srate/2),'bandpass');
    EMG_filt   = filtfilt(b_f,a_f,double(EMG_raw_r'))';
    z          = hilbert(double(EMG_filt'))';
    EMG_proc   = abs(z);
    clear EMG_raw_r EMG_filt z

    %% ── Re-epoch event loop ──────────────────────────────────────────────
    for ev = 1:n_reepoch
        ev_name  = reepoch_defs{ev,1};
        ev_col   = reepoch_defs{ev,2};
        win_s    = reepoch_defs{ev,3};
        win_ms   = win_s * 1000;
        n_smp_re = round(diff(win_s)*eeg_srate) + 1;
        t_re_ms  = linspace(win_ms(1), win_ms(2), n_smp_re);
        eeg_s1_off = round(win_s(1)*eeg_srate);
        eeg_s2_off = round(win_s(2)*eeg_srate);

        fprintf('  B | %s\n', ev_name);
        n_v = length(valid_idx_B);

        %% ── Synergy H re-epochs (from continuous EMG via pinv(W)) ────────
        H_reepoch = NaN(k, n_smp_re, n_v);

        for kk = 1:n_v
            oi = valid_idx_B(kk);
            mi = matched_vector(oi);
            ev_lat_ms = tw_latencies(oi, ev_col);
            if isnan(ev_lat_ms), continue; end

            ev_abs_s = FlxS_abs_time_s(mi) + ev_lat_ms/1000;
            s1 = round((ev_abs_s + win_s(1) - margin_s)*emg_srate) + 1;
            s2 = round((ev_abs_s + win_s(2) + margin_s)*emg_srate);
            if s1 < 1 || s2 > size(EMG_proc,2), continue; end

            emg_win = EMG_proc(:, s1:s2);
            emg_t   = linspace(win_ms(1)-margin_s*1000, ...
                               win_ms(2)+margin_s*1000, size(emg_win,2));
            emg_t   = emg_t - delta_t_ms(oi);

            % Synergy activation via pseudo-inverse
            H_win_raw = max(0, W_pinv * emg_win);        % [k x n_emg]
            H_reepoch(:, :, kk) = interp1(emg_t, H_win_raw', ...
                t_re_ms, 'linear', NaN)';
        end

        nan_re   = squeeze(any(any(isnan(H_reepoch),1),2));
        keep_re  = ~nan_re;
        H_re_ok  = H_reepoch(:,:,keep_re);
        orig_ok  = valid_idx_B(keep_re);
        clear H_reepoch

        if sum(keep_re) < 10
            fprintf('  Too few valid epochs for %s (%d) — skipping.\n', ...
                ev_name, sum(keep_re)); continue
        end
        fprintf('  B | %s: %d valid epochs\n', ev_name, sum(keep_re));

        %% ── EEG cluster loop ─────────────────────────────────────────────
        for ci = 1:n_clusters
            cname = eeg_clusters{ci};
            if ~cluster_present(ci), continue; end

            fig_path = fullfile(results_path, ...
                sprintf('B__sub%02d__%s__reepoch_%s.png', ...
                subject_id, cname, ev_name));
            if exist(fig_path,'file'), continue; end

            ic_id  = cluster_ic_ids(ci);
            ic_row = find(ic_ids_needed == ic_id);

            X_reepoch = NaN(n_smp_re, length(orig_ok));
            for kk = 1:length(orig_ok)
                oi = orig_ok(kk);
                ev_lat_ms = tw_latencies(oi, ev_col);
                if isnan(ev_lat_ms) || isnan(FlxS_cont_smp(oi)), continue; end
                ev_cont = FlxS_cont_smp(oi) + round(ev_lat_ms*eeg_srate/1000);
                s1 = ev_cont + eeg_s1_off;
                s2 = ev_cont + eeg_s2_off;
                if s1 < 1 || s2 > n_cont_smp, continue; end
                X_reepoch(:, kk) = ic_cont_all(ic_row, s1:s2)';
            end

            nan_eeg = any(isnan(X_reepoch),1);
            X_valid = X_reepoch(:, ~nan_eeg);
            Y_valid = H_re_ok(:, :, ~nan_eeg);

            if size(X_valid,2) < 10
                warning('B | sub%d | %s | %s: only %d epochs', ...
                    subject_id, cname, ev_name, size(X_valid,2)); continue
            end

            try
                [TF_X,TF_Y,tax,freqs] = compute_TF_locked( ...
                    X_valid, Y_valid, tf_params, win_ms);

                if isempty(reepoch_tax_ref{ev})
                    reepoch_tax_ref{ev} = tax;
                    if isempty(freqs_ref), freqs_ref = freqs; end
                end

                cmc_ev             = cmc_opts;
                cmc_ev.title       = sprintf('sub%02d  %s  %s  %s  (%d ep)', ...
                    subject_id, strrep(cname,'_',' '), ev_name, ...
                    '(synergies)', size(X_valid,2));
                cmc_ev.event_label = ev_name;

                [coh_b, fig_h] = compute_and_plot_CMC_eventlocked( ...
                    TF_X, TF_Y, tax, freqs, cmc_ev);
                clear TF_X TF_Y

                all_coh_B{si, ci, ev} = coh_b;
                if ~isempty(fig_h) && isvalid(fig_h)
                    try; exportgraphics(fig_h,fig_path,'Resolution',150);
                    catch; saveas(fig_h,fig_path,'png'); end
                    close(fig_h);
                end
                fprintf('  B | %s | %s — done.\n', ev_name, cname);
            catch ME
                warning('B | sub%d | %s | %s: %s', ...
                    subject_id, cname, ev_name, ME.message);
            end
        end % cluster loop
        clear H_re_ok
    end % event loop

    clear EMG_proc ic_cont_all

    %% ── Save per-subject CMC ─────────────────────────────────────────────
    coh_A = all_coh_A(si,:);
    coh_B = all_coh_B(si,:,:);
    save(fullfile(results_path, sprintf('sub%02d__coh_synergy.mat',subject_id)), ...
        'coh_A','coh_B','pct_grid_ref','freqs_ref','reepoch_tax_ref', ...
        'subject_id','-v7.3');
    fprintf('\n  Subject %d saved.\n', subject_id);

end % subject loop

fprintf('\nAll subjects done. Computing group averages...\n');

%% ══════════════════════════════════════════════════════════════════════════
%  GROUP AVERAGES
% ══════════════════════════════════════════════════════════════════════════

%% ── Section A: cycle-% ───────────────────────────────────────────────────
for ci = 1:n_clusters
    cname = eeg_clusters{ci};
    coh_stack = []; valid_subs = [];
    for si = 1:n_subjects
        if isempty(all_coh_A{si,ci}), continue; end
        coh_stack  = cat(4, coh_stack, all_coh_A{si,ci});
        valid_subs = [valid_subs, subject_ids(si)]; %#ok<AGROW>
    end
    if isempty(coh_stack), continue; end
    n_v = numel(valid_subs);
    coh_mean = mean(coh_stack,4,'omitnan');
    coh_sem  = std(coh_stack,0,4,'omitnan')/sqrt(n_v);
    save(fullfile(results_path,sprintf('group__coh_A__%s.mat',cname)), ...
        'coh_mean','coh_sem','pct_grid_ref','freqs_ref','valid_subs','-v7.3');
    ttl = sprintf('Synergy CMC (cycle-%%)\n%s   n=%d', ...
        strrep(cname,'_',' '), n_v);
    fig_h = plot_group_CMC(coh_mean, pct_grid_ref, freqs_ref, syn_labels, ttl);
    fp = fullfile(results_path, sprintf('group__A__%s.png',cname));
    try; exportgraphics(fig_h,fp,'Resolution',150); catch; saveas(fig_h,fp,'png'); end
    close(fig_h);
    fprintf('  A | %s: n=%d\n', cname, n_v);
end

%% ── Section B: re-epoch ──────────────────────────────────────────────────
for ev = 1:n_reepoch
    ev_name = reepoch_defs{ev,1};
    tax     = reepoch_tax_ref{ev};
    if isempty(tax), continue; end

    for ci = 1:n_clusters
        cname = eeg_clusters{ci};
        coh_stack = []; valid_subs = [];
        for si = 1:n_subjects
            if isempty(all_coh_B{si,ci,ev}), continue; end
            coh_stack  = cat(4, coh_stack, all_coh_B{si,ci,ev});
            valid_subs = [valid_subs, subject_ids(si)]; %#ok<AGROW>
        end
        if isempty(coh_stack), continue; end
        n_v = numel(valid_subs);
        coh_mean = mean(coh_stack,4,'omitnan');
        save(fullfile(results_path, ...
            sprintf('group__coh_B__%s__%s.mat',ev_name,cname)), ...
            'coh_mean','tax','freqs_ref','valid_subs','-v7.3');
        ttl = sprintf('Synergy CMC  %s | %s   n=%d', ...
            ev_name, strrep(cname,'_',' '), n_v);
        cmc_grp             = cmc_opts;
        cmc_grp.title       = ttl;
        cmc_grp.event_label = ev_name;
        fig_h = plot_group_CMC_event(coh_mean, tax, freqs_ref, syn_labels, ttl);
        fp = fullfile(results_path, ...
            sprintf('group__B__%s__%s.png',ev_name,cname));
        try; exportgraphics(fig_h,fp,'Resolution',150); catch; saveas(fig_h,fp,'png'); end
        close(fig_h);
        fprintf('  B | %s | %s: n=%d\n', ev_name, cname, n_v);
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
        {EMG_epoched.epoch.eventlatency},'UniformOutput',false);
    FlxS_type_lgl = cellfun(@(x) strcmpi(x,'FlxS'), eventtypes_EMG, ...
        'UniformOutput',false);
    FlxS_lat_lgl  = cellfun(@(x) x==0, eventlatency_EMG,'UniformOutput',false);
    FlxS_idx_EMG  = cellfun(@(x,y) find(x & y), FlxS_type_lgl, ...
        FlxS_lat_lgl,'UniformOutput',false);
    more_one = cellfun(@(x) numel(x)>1, FlxS_idx_EMG);
    init_time_cell = cellfun(@(x) cell2mat(x), ...
        {EMG_epoched.epoch.eventinit_time_EEG},'UniformOutput',false);
    init_time_EMG  = cellfun(@(x,y) x(y(1)), init_time_cell, FlxS_idx_EMG, ...
        'UniformOutput',false);
    init_time_EMG  = cell2mat(init_time_EMG);
    N_EEG = length(EEG_epoched_main.epoch);
    matched_vector   = NaN(N_EEG,1);
    condition_vector = zeros(N_EEG,1);
    for i = 1:N_EEG
        evtypes = EEG_epoched_main.epoch(i).eventtype;
        evlat   = cell2mat(EEG_epoched_main.epoch(i).eventlatency);
        FlxS_i  = find(strcmpi(evtypes,'FlxS') & evlat==0,1);
        if isempty(FlxS_i), continue; end
        t_EEG = EEG_epoched_main.epoch(i).eventinit_time{FlxS_i};
        hit   = find(abs(init_time_EMG - t_EEG) < 1e-10);
        if numel(hit)>1 || any(more_one(hit)) || isempty(hit)
            matched_vector(i) = NaN;
        else
            matched_vector(i)   = hit;
            condition_vector(i) = str2double( ...
                EEG_epoched_main.epoch(i).eventcond{FlxS_i});
        end
    end
end


function fig_h = plot_group_CMC(coh_mean, pct_grid, freqs, ...
        syn_labels, fig_title)
    nS = size(coh_mean,3);
    nR = ceil(sqrt(nS)); nC = ceil(nS/nR);
    cmax = prctile(coh_mean(:),95); if isnan(cmax)||cmax==0, cmax=1; end
    fig_h = figure('Color','w','Name',fig_title);
    for s = 1:nS
        ax = subplot(nR,nC,s);
        contourf(pct_grid, freqs, coh_mean(:,:,s), 200,'LineColor','none');
        set(ax,'yscale','log','ydir','norm', ...
            'ylim',[freqs(1) freqs(end)]+[0 .01],'xlim',[-5 105], ...
            'ytick',[4 8 14 30 60],'box','on','FontSize',12);
        colormap(ax,'turbo'); caxis([0 cmax]);
        hold on; xline(0,'w--','LineWidth',1.5); xline(100,'w--','LineWidth',1.5); hold off
        title(syn_labels{s},'FontWeight','bold'); xlabel('Cycle (%)'); ylabel('Frequency (Hz)');
    end
    cb = colorbar(subplot(nR,nC,nS),'Position',[0.93 0.13 0.015 0.755]);
    ylabel(cb,'MSC'); sgtitle(fig_title,'FontSize',14,'FontWeight','normal');
end


function fig_h = plot_group_CMC_event(coh_mean, tax, freqs, ...
        syn_labels, fig_title)
    nS = size(coh_mean,3);
    nR = ceil(sqrt(nS)); nC = ceil(nS/nR);
    cmax = prctile(coh_mean(:),95); if isnan(cmax)||cmax==0, cmax=1; end
    fig_h = figure('Color','w','Name',fig_title);
    for s = 1:nS
        ax = subplot(nR,nC,s);
        contourf(tax, freqs, coh_mean(:,:,s), 200,'LineColor','none');
        set(ax,'yscale','log','ydir','norm', ...
            'ylim',[freqs(1) freqs(end)]+[0 .01], ...
            'xlim',[tax(1) tax(end)],'ytick',[4 8 14 30 60],'box','on','FontSize',12);
        colormap(ax,'turbo'); caxis([0 cmax]);
        hold on; xline(0,'w--','LineWidth',2); hold off
        title(syn_labels{s},'FontWeight','bold');
        xlabel('Time relative to event (ms)'); ylabel('Frequency (Hz)');
    end
    cb = colorbar(subplot(nR,nC,nS),'Position',[0.93 0.13 0.015 0.755]);
    ylabel(cb,'MSC'); sgtitle(fig_title,'FontSize',14,'FontWeight','normal');
end