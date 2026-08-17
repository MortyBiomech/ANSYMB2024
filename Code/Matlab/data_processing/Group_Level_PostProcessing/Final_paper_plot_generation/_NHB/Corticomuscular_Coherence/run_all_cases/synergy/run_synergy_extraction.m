%% run_synergy_extraction.m
%
% Extract muscle synergies per subject using Non-negative Matrix Factorization
% (NMF) applied to the concatenated epoch matrix.
%
% EMG representation : band-pass 20-450 Hz → Hilbert amplitude  |H{x}|
%                      (non-negative by construction, required for NMF)
% Synergy number     : determined automatically by variance accounted for
%                      (VAF) criterion:
%                        global VAF >= 90%  AND  per-muscle VAF >= 75%
% NMF stability      : 20 random restarts per k; best solution kept
% Epoch set          : same valid epochs used in run_CMC_allsubjects.m,
%                      interpolated onto the EEG time grid (500 Hz)
%
% Outputs per subject (saved to results_path):
%   sub{id}__synergies.mat
%       .W          [4 x k_opt]            synergy weight matrix (normalised)
%       .H_epochs   [k_opt x n_smp x n_ep] per-epoch activation profiles
%       .k_opt      scalar                 number of synergies selected
%       .VAF_global [1 x 4]               global VAF for k = 1..4
%       .VAF_muscle [4 x 4]               per-muscle VAF for k = 1..4
%       .keep       [N_EEG x 1] logical   valid epoch mask (matches CMC scripts)
%       .emg_labels {'VL','RF','GM','BF'}
%       .subject_id scalar
%
% Figures saved:
%   sub{id}__synergy_VAF.png    VAF curves
%   sub{id}__synergy_weights.png W bar plots per synergy

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
cmc_folder_path  = [processing_path, 'Group_Level_PostProcessing\', ...
                    'Final_paper_plot_generation\_NHB\', ...
                    'Corticomuscular_Coherence\run_all_cases\'];
results_path     = [cmc_folder_path, 'CMC_results\synergy_extraction\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

%% ── Fixed parameters ─────────────────────────────────────────────────────
emg_labels     = {'VL', 'RF', 'GM', 'BF'};
emg_srate      = 2000;
eeg_srate      = 500;
epoch_limits_s = [-0.5, 3.5];
nChan          = numel(emg_ids);

% NMF parameters
k_max         = nChan;    % max synergies = number of muscles
n_replicates  = 20;       % random restarts for NMF stability
VAF_global_th = 0.85;     % 90% global VAF threshold
VAF_muscle_th = 0.75;     % 75% per-muscle VAF threshold

%% ══════════════════════════════════════════════════════════════════════════
%  SUBJECT LOOP
% ══════════════════════════════════════════════════════════════════════════
for si = 1:n_subjects

    subject_id = subject_ids(si);
    fprintf('\n════════  Subject %d  (%d / %d)  ════════\n', ...
        subject_id, si, n_subjects);

    % mat_path = fullfile(results_path, ...
    %     sprintf('sub%02d__synergies.mat', subject_id));
    mat_path = fullfile(results_path, ...
        sprintf('sub%02d__synergies_S3.mat', subject_id));
    if exist(mat_path, 'file')
        fprintf('  Already extracted — skipping.\n');
        continue
    end

    %% ── Load EEG (needed for timing alignment only) ──────────────────────
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
    EEG_times = EEG_epoched_main.times;

    %% ── Load raw EMG ─────────────────────────────────────────────────────
    try
        output  = runs_concatenated(subject_id, rawdata_path);
        EMG_raw = double(output.All_EMG(emg_ids, :));
        clear output
    catch ME
        warning('Subject %d: EMG load failed — %s', subject_id, ME.message);
        continue
    end

    %% ── Epoch timing alignment ───────────────────────────────────────────
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
    epoch_start_s   = NaN(size(EMG_epoched_meta.epoch));

    for i = 1:length(matched_vector)
        if isnan(matched_vector(i)), continue; end
        ep       = EMG_epoched_meta.epoch(matched_vector(i));
        FlxS_idx = find(strcmpi(ep.eventtype,'FlxS') & ...
            cell2mat(ep.eventlatency) == 0, 1);
        if isempty(FlxS_idx), continue; end
        t_last = ep.eventinit_time_lastEMGpoint{FlxS_idx};
        t_FlxS = ep.eventinit_time{FlxS_idx};
        delta_t_ms(i)                    = 1000 * (t_FlxS - t_last);
        epoch_start_s(matched_vector(i)) = t_FlxS + epoch_limits_s(1);
    end
    clear EMG_epoched_meta

    epoch_start_smp = round(epoch_start_s * emg_srate) + 1;
    epoch_end_smp   = epoch_start_smp + N_EMG_smp - 1;

    %% ── Apply EMG preprocessing: bp 20-450 Hz → Hilbert amplitude ────────
    [b, a]   = butter(4, [20 450] / (emg_srate/2), 'bandpass');
    EMG_filt = filtfilt(b, a, double(EMG_raw'))';
    z        = hilbert(double(EMG_filt'))';
    EMG_proc = abs(z);
    clear EMG_raw EMG_filt z

    %% ── Interpolate epochs onto EEG time grid ────────────────────────────
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
        EMG_intrp(:, :, i) = interp1(emg_t', emg_ep', EEG_times', 'linear')';
    end
    clear EMG_proc

    nan_eps        = squeeze(any(any(isnan(EMG_intrp), 1), 2));
    keep           = ~nan_eps;
    EMG_epochs_CMC = EMG_intrp(:, :, keep);
    clear EMG_intrp

    n_valid = sum(keep);
    n_smp   = size(EMG_epochs_CMC, 2);
    fprintf('  Valid epochs for NMF: %d\n', n_valid);

    if n_valid < 10
        warning('Subject %d: too few valid epochs (%d) — skipping.', ...
            subject_id, n_valid);
        continue
    end

    %% ── Concatenate epochs: [4 x (n_smp * n_valid)] ─────────────────────
    % Concatenate along the time dimension so NMF sees one long matrix.
    % H_concat can then be reshaped back to [k x n_smp x n_valid].
    EMG_concat = reshape(EMG_epochs_CMC, nChan, n_smp * n_valid);

    %% ── NMF for k = 1 : k_max ───────────────────────────────────────────
    fprintf('  Running NMF (k = 1 to %d, %d replicates each)...\n', ...
        k_max, n_replicates);

    VAF_global = zeros(1, k_max);
    VAF_muscle = zeros(nChan, k_max);
    W_cell     = cell(1, k_max);
    H_cell     = cell(1, k_max);

    SS_total = sum(EMG_concat(:).^2);

    for k = 3 %1:k_max
        [W_k, H_k] = nnmf(EMG_concat, k, ...
            'replicates', n_replicates, ...
            'algorithm',  'als');

        % Normalise W columns to unit norm (scale absorbed into H)
        for col = 1:k
            col_norm = norm(W_k(:, col));
            if col_norm > 0
                W_k(:, col) = W_k(:, col) / col_norm;
                H_k(col, :) = H_k(col, :) * col_norm;
            end
        end

        W_cell{k} = W_k;
        H_cell{k} = H_k;

        % Global VAF
        residual        = EMG_concat - W_k * H_k;
        VAF_global(k)   = 1 - sum(residual(:).^2) / SS_total;

        % Per-muscle VAF
        for m = 1:nChan
            ss_m           = sum(EMG_concat(m,:).^2);
            VAF_muscle(m,k)= 1 - sum(residual(m,:).^2) / ss_m;
        end

        fprintf('    k = %d : global VAF = %.1f%%   per-muscle VAF = %s\n', ...
            k, VAF_global(k)*100, ...
            strjoin(arrayfun(@(v) sprintf('%.0f%%', v*100), ...
            VAF_muscle(:,k)', 'UniformOutput', false), ' | '));
    end

    %% ── Select optimal k ─────────────────────────────────────────────────
    % % Minimum k satisfying: global VAF >= 90% AND all per-muscle VAF >= 75%
    % k_opt = [];
    % for k = 1:k_max
    %     % global VAF criterion only
    %     if VAF_global(k) >= VAF_global_th
    %         k_opt = k;
    %         break
    %     end
    % end
    % if isempty(k_opt)
    %     k_opt = k_max;
    %     warning('Subject %d: VAF threshold not reached — using k = %d.', ...
    %         subject_id, k_opt);
    % end

    k_opt = 3;   % fixed by biomechanical rationale

    fprintf('  Selected k = %d  (global VAF = %.1f%%)\n', ...
        k_opt, VAF_global(k_opt)*100);

    W_opt      = W_cell{k_opt};
    H_opt_flat = H_cell{k_opt};

    %% ── Reshape H back into per-epoch structure ──────────────────────────
    % H_opt_flat: [k_opt x (n_smp * n_valid)]
    % → H_epochs:  [k_opt x n_smp x n_valid]
    H_epochs = reshape(H_opt_flat, k_opt, n_smp, n_valid);

    %% ── Save synergy data ────────────────────────────────────────────────
    synergy_data.W          = W_opt;
    synergy_data.H_epochs   = H_epochs;
    synergy_data.k_opt      = k_opt;
    synergy_data.VAF_global = VAF_global;
    synergy_data.VAF_muscle = VAF_muscle;
    synergy_data.keep       = keep;
    synergy_data.emg_labels = emg_labels;
    synergy_data.subject_id = subject_id;

    save(mat_path, 'synergy_data', '-v7.3');
    fprintf('  Synergy data saved.\n');

    % %% ── Figure 1: VAF curves ─────────────────────────────────────────────
    % fig_vaf = figure('Color', 'w', 'Name', ...
    %     sprintf('sub%02d VAF', subject_id));
    % 
    % subplot(1, 2, 1); hold on
    % plot(1:k_max, VAF_global*100, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
    % yline(VAF_global_th*100, 'r--', 'LineWidth', 1.5, ...
    %     'Label', '90% threshold', 'HandleVisibility', 'off');
    % xline(k_opt, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % xlabel('Number of synergies  k');
    % ylabel('Global VAF (%)');
    % title('Global variance accounted for');
    % set(gca, 'XTick', 1:k_max, 'FontSize', 12, 'Box', 'on');
    % ylim([0 105])
    % 
    % subplot(1, 2, 2); hold on
    % colors = lines(nChan);
    % for m = 1:nChan
    %     plot(1:k_max, VAF_muscle(m,:)*100, 'o-', ...
    %         'Color', colors(m,:), 'LineWidth', 1.5, 'MarkerSize', 7, ...
    %         'DisplayName', emg_labels{m});
    % end
    % yline(VAF_muscle_th*100, 'r--', 'LineWidth', 1.5, ...
    %     'Label', '75% threshold', 'HandleVisibility', 'off');
    % xline(k_opt, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % xlabel('Number of synergies  k');
    % ylabel('Per-muscle VAF (%)');
    % title('Per-muscle variance accounted for');
    % legend('Location', 'southeast');
    % set(gca, 'XTick', 1:k_max, 'FontSize', 12, 'Box', 'on');
    % ylim([0 105])
    % 
    % sgtitle(sprintf('Sub %d — Synergy VAF   (selected k = %d)', ...
    %     subject_id, k_opt), 'FontSize', 14);
    % 
    % fig_vaf_path = fullfile(results_path, ...
    %     sprintf('sub%02d__synergy_VAF.png', subject_id));
    % try; exportgraphics(fig_vaf, fig_vaf_path, 'Resolution', 150);
    % catch; saveas(fig_vaf, fig_vaf_path, 'png'); end
    % close(fig_vaf);

    %% ── Figure 2: Synergy weight matrix W ───────────────────────────────
    fig_w = figure('Color', 'w', 'Name', ...
        sprintf('sub%02d Synergy Weights', subject_id));

    for k = 1:k_opt
        subplot(1, k_opt, k)
        barh(W_opt(:, k), 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
        set(gca, 'YTick', 1:nChan, 'YTickLabel', emg_labels, ...
            'FontSize', 12, 'Box', 'on');
        xlabel('Weight');
        xlim([0, max(W_opt(:)) * 1.15]);
        title(sprintf('Synergy %d', k), 'FontWeight', 'bold');
        if k == 1, ylabel('Muscle'); end
    end

    sgtitle(sprintf('Sub %d — Synergy weights W', subject_id), 'FontSize', 14);

    % fig_w_path = fullfile(results_path, ...
    %     sprintf('sub%02d__synergy_weights.png', subject_id));
    fig_w_path = fullfile(results_path, ...
        sprintf('sub%02d__synergy_weights_S3.png', subject_id));
    try; exportgraphics(fig_w, fig_w_path, 'Resolution', 150);
    catch; saveas(fig_w, fig_w_path, 'png'); end
    close(fig_w);

end % subject loop

fprintf('\nSynergy extraction complete.\nResults in: %s\n', results_path);


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