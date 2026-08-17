%% run_synergy_extraction_group.m
%
% Across-subject NMF for muscle synergy extraction — k = 2 and k = 3.
%
% k = 2  expected synergies:
%         extension  (VL + RF dominant)
%         flexion    (BF + GM dominant)
%
% k = 3  expected synergies:
%         extension      (VL + RF dominant)
%         flexion        (BF + GM dominant)
%         co-contraction (all muscles, intermediate extension score)
%
% Sorting logic (works for any k):
%   extension score = (VL+RF weights) − (BF+GM weights)
%   sort descending:  highest → extension
%                     lowest  → flexion
%                     middle  → co-contraction (k=3 only)
%
% Since W is shared across subjects, sorting is consistent everywhere —
% no per-subject matching is needed.
%
% Outputs per k value:
%   group__synergies_k{k}.mat
%   sub{id}__synergies_group_k{k}.mat
%   group__synergy_W_k{k}.png
%   group__VAF_validation_k{k}.png
%   group__VAF_comparison.png          (k=2 vs k=3 side by side)

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
results_path     = [cmc_folder_path, ...
                    'CMC_results\synergy_extraction_group\'];
if ~exist(results_path, 'dir'), mkdir(results_path); end

this_path = pwd;
cd(EEGLAB_path)
if ~exist('ALLCOM', 'var'), eeglab; end
cd(this_path)

%% ── Fixed parameters ─────────────────────────────────────────────────────
emg_labels   = {'VL', 'RF', 'GM', 'BF'};
n_muscles    = numel(emg_ids);
emg_srate    = 2000;
eeg_srate    = 500;
epoch_limits_s = [-0.5, 3.5];
k_values     = [2, 3];          % run NMF for both
n_replicates = 20;

% Synergy labels per k value
synergy_labels_k = { ...
    {'extension', 'flexion'}, ...                           % k=2
    {'extension', 'flexion', 'co-contraction'} };          % k=3

% Colours: extension=blue, flexion=red, co-contraction=green
synergy_colors = [0.2 0.55 0.85; 0.85 0.3 0.2; 0.2 0.7 0.35];

%% ══════════════════════════════════════════════════════════════════════════
%  PHASE 1 — Extract valid EMG epochs (done once, shared by all k values)
% ══════════════════════════════════════════════════════════════════════════
fprintf('Phase 1: Extracting epochs from all subjects...\n\n');

EMG_epoch_cell = cell(n_subjects, 1);
keep_cell      = cell(n_subjects, 1);
n_smp_per_sub  = zeros(n_subjects, 1);
n_ep_per_sub   = zeros(n_subjects, 1);

for si = 1:n_subjects

    subject_id = subject_ids(si);
    fprintf('  Subject %d (%d/%d)...\n', subject_id, si, n_subjects);

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

    try
        output  = runs_concatenated(subject_id, rawdata_path);
        EMG_raw = double(output.All_EMG(emg_ids, :));
        clear output
    catch ME
        warning('Subject %d: EMG load failed — %s', subject_id, ME.message);
        continue
    end

    % Timing alignment
    [b_tmp, a_tmp] = butter(4, [20 450]/(emg_srate/2), 'bandpass');
    EMG_filt_tmp   = filtfilt(b_tmp, a_tmp, double(EMG_raw'))';

    eeg_emptyset = structfun(@(x) [], EEG_orig, 'UniformOutput', false);
    EMG_tmpl = eeg_emptyset;
    EMG_tmpl.data = double(EMG_filt_tmp); EMG_tmpl.nbchan = n_muscles;
    EMG_tmpl.pnts = size(EMG_filt_tmp,2); EMG_tmpl.trials = 1;
    EMG_tmpl.srate = emg_srate; EMG_tmpl.xmin = 0;
    EMG_tmpl.xmax = (EMG_tmpl.pnts-1)/emg_srate;
    EMG_tmpl.setname = 'EMG'; EMG_tmpl.filename = ''; EMG_tmpl.filepath = '';
    for ch = 1:n_muscles
        EMG_tmpl.chanlocs(ch).labels = emg_labels{ch};
    end
    EMG_tmpl = eeg_checkset(EMG_tmpl);

    output_timing = runs_concatenated(subject_id, rawdata_path);
    [EMG_ev, ~] = import_same_event_on_EMG_stream( ...
        EMG_tmpl, EEG_orig, output_timing, subject_id, processing_path, data_path);
    clear EMG_tmpl EMG_filt_tmp output_timing

    for e = 1:length(EMG_ev.event)
        nums = str2double(strsplit(EMG_ev.event(e).desc,'_')).';
        if strcmp(EMG_ev.event(e).type,'boundary')
            EMG_ev.event(e).trial = 'none'; EMG_ev.event(e).cond = 'none';
        else
            EMG_ev.event(e).trial = nums(end); EMG_ev.event(e).cond = nums(2);
        end
    end

    EMG_ev_sel = pop_selectevent(EMG_ev, ...
        'type',{'FlxS','ExtS','ExtE'}, 'deleteevents','on');
    EMG_epoched_meta = pop_epoch(EMG_ev_sel, {'FlxS'}, epoch_limits_s, ...
        'newname','meta','epochinfo','yes');
    EMG_epoched_meta = eeg_checkset(EMG_epoched_meta);
    clear EMG_ev EMG_ev_sel

    [matched_vector, ~] = find_matched_epochs(EEG_epoched_main, EMG_epoched_meta);

    EMG_time_meta = EMG_epoched_meta.times;
    N_EMG_smp     = EMG_epoched_meta.pnts;
    delta_t_ms    = NaN(length(matched_vector), 1);
    epoch_start_s = NaN(size(EMG_epoched_meta.epoch));

    for i = 1:length(matched_vector)
        if isnan(matched_vector(i)), continue; end
        ep = EMG_epoched_meta.epoch(matched_vector(i));
        FlxS_idx = find(strcmpi(ep.eventtype,'FlxS') & ...
            cell2mat(ep.eventlatency)==0, 1);
        if isempty(FlxS_idx), continue; end
        t_last = ep.eventinit_time_lastEMGpoint{FlxS_idx};
        t_FlxS = ep.eventinit_time{FlxS_idx};
        delta_t_ms(i)                    = 1000*(t_FlxS - t_last);
        epoch_start_s(matched_vector(i)) = t_FlxS + epoch_limits_s(1);
    end
    clear EMG_epoched_meta

    epoch_start_smp = round(epoch_start_s * emg_srate) + 1;
    epoch_end_smp   = epoch_start_smp + N_EMG_smp - 1;

    % Bandpass → Hilbert amplitude
    [b, a]   = butter(4, [20 450]/(emg_srate/2), 'bandpass');
    EMG_filt = filtfilt(b, a, double(EMG_raw'))';
    z        = hilbert(double(EMG_filt'))';
    EMG_proc = abs(z);
    clear EMG_raw EMG_filt z

    % Interpolate onto EEG time grid
    N_EEG     = length(EEG_epoched_main.epoch);
    EMG_intrp = NaN(n_muscles, length(EEG_times), N_EEG);
    for i = 1:N_EEG
        mi = matched_vector(i);
        if isnan(mi) || isnan(delta_t_ms(i)), continue; end
        s1 = epoch_start_smp(mi); s2 = epoch_end_smp(mi);
        if s1 < 1 || s2 > size(EMG_proc,2), continue; end
        emg_ep = EMG_proc(:, s1:s2);
        emg_t  = EMG_time_meta - delta_t_ms(i);
        EMG_intrp(:,:,i) = interp1(emg_t', emg_ep', EEG_times', 'linear')';
    end
    clear EMG_proc


    nan_eps = squeeze(any(any(isnan(EMG_intrp),1),2));
    keep    = ~nan_eps;

    % Maximum normalisation per muscle
    % Use the 99th percentile instead of the absolute max
    % to protect against single-sample noise spikes
    EMG_valid = EMG_intrp(:, :, keep);   % [4 x n_smp x n_valid]
    
    for m = 1:n_muscles
        muscle_vals = EMG_valid(m, :, :);
        scale_m     = prctile(muscle_vals(:), 99);
        if scale_m > 0
            EMG_valid(m, :, :) = EMG_valid(m, :, :) / scale_m;
        end
    end
    
    EMG_epoch_cell{si} = EMG_valid;
    keep_cell{si}      = keep;
    n_smp_per_sub(si)  = size(EMG_intrp, 2);
    n_ep_per_sub(si)   = sum(keep);

    clear EMG_intrp

    fprintf('    Valid epochs: %d\n', n_ep_per_sub(si));


    
end

%% ══════════════════════════════════════════════════════════════════════════
%  PHASE 2 — Build concatenated matrix (shared across all k values)
% ══════════════════════════════════════════════════════════════════════════
fprintf('\nPhase 2: Building concatenated EMG matrix...\n');

valid_subs = find(~cellfun(@isempty, EMG_epoch_cell));
col_ranges = zeros(n_subjects, 2);
EMG_concat = [];

for si = 1:numel(valid_subs)
    n_smp = n_smp_per_sub(si);
    n_ep  = n_ep_per_sub(si);
    flat  = reshape(EMG_epoch_cell{si}, n_muscles, n_smp * n_ep);
    col_start = size(EMG_concat, 2) + 1;
    EMG_concat = [EMG_concat, flat]; %#ok<AGROW>
    col_ranges(si, :) = [col_start, size(EMG_concat, 2)];
end

fprintf('  Concatenated matrix: [%d x %d]  (%d subjects)\n', ...
    size(EMG_concat,1), size(EMG_concat,2), numel(valid_subs));

%% ══════════════════════════════════════════════════════════════════════════
%  PHASE 3 — NMF loop over k values
% ══════════════════════════════════════════════════════════════════════════
VAF_global_all = NaN(n_subjects, numel(k_values));
VAF_muscle_all = NaN(n_muscles, n_subjects, numel(k_values));
W_all          = cell(numel(k_values), 1);

for ki = 1:numel(k_values)

    k = k_values(ki);
    syn_labels = synergy_labels_k{ki};
    fprintf('\n══ k = %d ══════════════════════════════════════════\n', k);

    %% ── Run NMF ──────────────────────────────────────────────────────────
    fprintf('  Running NMF (%d replicates)...\n', n_replicates);
    [W, H_concat] = nnmf(EMG_concat, k, ...
        'replicates', n_replicates, 'algorithm', 'als');

    % Normalise columns
    for col = 1:k
        cn = norm(W(:,col));
        if cn > 0
            W(:,col)        = W(:,col) / cn;
            H_concat(col,:) = H_concat(col,:) * cn;
        end
    end

    %% ── Sort synergies ───────────────────────────────────────────────────
    % Extension score = (VL+RF) - (BF+GM)
    % Sort descending: highest → extension, lowest → flexion, middle → co-contraction
    ext_score         = (W(1,:) + W(2,:)) - (W(3,:) + W(4,:));
    [~, score_order]  = sort(ext_score, 'descend');
    % Map sorted positions to synergy slots: 1=ext, end=flx, middle=cocon
    sort_order        = [score_order(1), score_order(end), score_order(2:end-1)];
    sort_order        = sort_order(1:k);   % trim to k (handles k=2: no middle)

    W        = W(:, sort_order);
    H_concat = H_concat(sort_order, :);
    W_all{ki} = W;

    fprintf('  Extension score per synergy after sorting:\n');
    for s = 1:k
        fprintf('    %s : %.3f  (VL=%.2f RF=%.2f GM=%.2f BF=%.2f)\n', ...
            syn_labels{s}, ext_score(sort_order(s)), ...
            W(1,s), W(2,s), W(3,s), W(4,s));
    end

    %% ── Slice H and save per subject ─────────────────────────────────────
    fprintf('  Saving per-subject activations...\n');
    for si = 1:numel(valid_subs)
        subject_id = subject_ids(si);
        n_smp = n_smp_per_sub(si);
        n_ep  = n_ep_per_sub(si);
        H_flat   = H_concat(:, col_ranges(si,1):col_ranges(si,2));
        H_epochs = reshape(H_flat, k, n_smp, n_ep);
        keep     = keep_cell{si};

        mat_path = fullfile(results_path, ...
            sprintf('sub%02d__synergies_group_k%d.mat', subject_id, k));
        save(mat_path, 'H_epochs', 'keep', 'subject_id', ...
            'W', 'k', 'syn_labels', 'emg_labels', '-v7.3');
    end

    %% ── Validate: per-subject VAF ────────────────────────────────────────
    fprintf('  Computing per-subject VAF...\n');
    for si = 1:numel(valid_subs)
        n_smp  = n_smp_per_sub(si);
        n_ep   = n_ep_per_sub(si);
        E_flat = reshape(EMG_epoch_cell{si}, n_muscles, n_smp * n_ep);
        H_flat = H_concat(:, col_ranges(si,1):col_ranges(si,2));
        resid  = E_flat - W * H_flat;

        VAF_global_all(si, ki) = 1 - sum(resid(:).^2)/sum(E_flat(:).^2);
        for m = 1:n_muscles
            ss_m = sum(E_flat(m,:).^2);
            VAF_muscle_all(m, si, ki) = 1 - sum(resid(m,:).^2)/ss_m;
        end

        fprintf('    sub%02d: global = %.1f%%   per-muscle: %s\n', ...
            subject_ids(si), VAF_global_all(si,ki)*100, ...
            strjoin(arrayfun(@(v) sprintf('%.0f%%',v*100), ...
            VAF_muscle_all(:,si,ki)', 'UniformOutput',false), ' | '));
    end

    %% ── Save group summary ───────────────────────────────────────────────
    VAF_global = VAF_global_all(:, ki);
    VAF_muscle = VAF_muscle_all(:, :, ki);
    save(fullfile(results_path, sprintf('group__synergies_k%d.mat', k)), ...
        'W', 'syn_labels', 'emg_labels', 'VAF_global', 'VAF_muscle', ...
        'subject_ids', 'k', '-v7.3');

    %% ── Figure: synergy weights W ────────────────────────────────────────
    fig_w = figure('Color','w','Name',sprintf('W k=%d',k));
    for s = 1:k
        subplot(1, k, s)
        barh(1:n_muscles, W(:,s), 0.5, ...
            'FaceColor', synergy_colors(s,:), 'EdgeColor','none');
        set(gca, 'YTick',1:n_muscles, 'YTickLabel',emg_labels, ...
            'XLim',[0 1.1], 'FontSize',13, 'Box','on');
        xlabel('Weight');
        title(syn_labels{s}, 'FontWeight','bold', 'FontSize',14);
        if s == 1, ylabel('Muscle'); end
    end
    sgtitle(sprintf('Group synergy weights W   (k = %d)', k), ...
        'FontSize',14, 'FontWeight','normal');
    try; exportgraphics(fig_w, fullfile(results_path, ...
        sprintf('group__synergy_W_k%d.png',k)), 'Resolution',150);
    catch; saveas(fig_w, fullfile(results_path, ...
        sprintf('group__synergy_W_k%d.png',k)),'png'); 
    end
    close(fig_w);

    %% ── Figure: per-subject VAF validation ───────────────────────────────
    fig_vaf = figure('Color','w','Name',sprintf('VAF k=%d',k));
    valid_mask = ~isnan(VAF_global);
    x_idx      = find(valid_mask);
    x_labels   = arrayfun(@(x) sprintf('sub%02d',x), ...
        subject_ids(valid_mask), 'UniformOutput',false);

    subplot(1,2,1)
    bar(x_idx, VAF_global(valid_mask)*100, 0.6, ...
        'FaceColor',[0.4 0.7 0.4], 'EdgeColor','none');
    yline(90,'r--','LineWidth',1.5,'HandleVisibility','off','Label','90%');
    set(gca,'XTick',x_idx,'XTickLabel',x_labels, ...
        'XTickLabelRotation',45,'YLim',[0 105],'FontSize',11,'Box','on');
    ylabel('Global VAF (%)');
    title(sprintf('Global VAF  (k=%d)',k), 'FontWeight','bold');

    subplot(1,2,2); hold on
    colors_m = lines(n_muscles);
    for m = 1:n_muscles
        plot(x_idx, VAF_muscle_all(m,valid_mask,ki)*100, 'o-', ...
            'Color',colors_m(m,:),'LineWidth',1.5,'MarkerSize',6, ...
            'DisplayName',emg_labels{m});
    end
    yline(75,'r--','LineWidth',1.5,'HandleVisibility','off','Label','75%');
    set(gca,'XTick',x_idx,'XTickLabel',x_labels, ...
        'XTickLabelRotation',45,'YLim',[0 105],'FontSize',11,'Box','on');
    ylabel('Per-muscle VAF (%)');
    title(sprintf('Per-muscle VAF  (k=%d)',k),'FontWeight','bold');
    legend('Location','southeast','FontSize',10);

    sgtitle(sprintf('Shared W validation — k = %d', k), ...
        'FontSize',14,'FontWeight','normal');
    try; exportgraphics(fig_vaf, fullfile(results_path, ...
        sprintf('group__VAF_validation_k%d.png',k)),'Resolution',150);
    catch; saveas(fig_vaf, fullfile(results_path, ...
        sprintf('group__VAF_validation_k%d.png',k)),'png'); 
    end
    close(fig_vaf);

end % k loop

%% ══════════════════════════════════════════════════════════════════════════
%  PHASE 4 — Comparison figure: k=2 vs k=3 VAF and W side by side
% ══════════════════════════════════════════════════════════════════════════
fprintf('\nPhase 4: Generating comparison figures...\n');

valid_mask = ~isnan(VAF_global_all(:,1));
x_idx      = find(valid_mask);
x_labels   = arrayfun(@(x) sprintf('sub%02d',x), ...
    subject_ids(valid_mask), 'UniformOutput',false);

%% ── VAF comparison: global k=2 vs k=3 per subject ────────────────────────
fig_cmp = figure('Color','w','Name','VAF comparison k=2 vs k=3');

subplot(1,2,1); hold on
bar_data = [VAF_global_all(valid_mask,1)*100, ...
            VAF_global_all(valid_mask,2)*100];
b = bar(x_idx, bar_data, 0.7);
b(1).FaceColor = [0.4 0.55 0.8];   % k=2
b(2).FaceColor = [0.85 0.55 0.2];  % k=3
yline(90,'r--','LineWidth',1.5,'HandleVisibility','off','Label','90%');
set(gca,'XTick',x_idx,'XTickLabel',x_labels, ...
    'XTickLabelRotation',45,'YLim',[0 105],'FontSize',11,'Box','on');
legend({'k=2','k=3'},'Location','southeast');
ylabel('Global VAF (%)');
title('Global VAF comparison','FontWeight','bold');

% VAF gain from adding the 3rd synergy
subplot(1,2,2)
gain = (VAF_global_all(valid_mask,2) - VAF_global_all(valid_mask,1))*100;
bar(x_idx, gain, 0.6, 'FaceColor',[0.3 0.75 0.45], 'EdgeColor','none');
set(gca,'XTick',x_idx,'XTickLabel',x_labels, ...
    'XTickLabelRotation',45,'YLim',[0, max(gain)*1.3+1], ...
    'FontSize',11,'Box','on');
ylabel('\DeltaVAF (%)');
title('VAF gain: k=3 over k=2','FontWeight','bold');
yline(mean(gain),'k--','LineWidth',1.5,'Label', ...
    sprintf('mean = %.1f%%', mean(gain)),'HandleVisibility','off');

sgtitle('Group NMF: k=2 vs k=3', 'FontSize',14,'FontWeight','normal');
try; exportgraphics(fig_cmp, fullfile(results_path,'group__VAF_comparison.png'), ...
    'Resolution',150);
catch; saveas(fig_cmp, fullfile(results_path,'group__VAF_comparison.png'),'png'); 
end
close(fig_cmp);

%% ── W comparison: all synergies from k=2 and k=3 on one figure ───────────
fig_wall = figure('Color','w','Name','W comparison all k');
k_tot = sum(k_values);   % 2+3=5 subplots total
sp_idx = 0;

for ki = 1:numel(k_values)
    k          = k_values(ki);
    W          = W_all{ki};
    syn_labels = synergy_labels_k{ki};
    for s = 1:k
        sp_idx = sp_idx + 1;
        subplot(numel(k_values), max(k_values), ...
            (ki-1)*max(k_values) + s)
        barh(1:n_muscles, W(:,s), 0.5, ...
            'FaceColor', synergy_colors(s,:), 'EdgeColor','none');
        set(gca,'YTick',1:n_muscles,'YTickLabel',emg_labels, ...
            'XLim',[0 1.1],'FontSize',11,'Box','on');
        xlabel('Weight');
        title(sprintf('k=%d  %s', k, syn_labels{s}), ...
            'FontWeight','bold','FontSize',12);
        if s == 1, ylabel('Muscle'); end
    end
end

sgtitle('Group synergy weights W: k=2 (top) vs k=3 (bottom)', ...
    'FontSize',14,'FontWeight','normal');
try; exportgraphics(fig_wall, fullfile(results_path,'group__W_comparison.png'), ...
    'Resolution',150);
catch; saveas(fig_wall, fullfile(results_path,'group__W_comparison.png'),'png'); 
end
close(fig_wall);

fprintf('\nAll done. Results in:\n  %s\n', results_path);


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
    FlxS_lat_lgl  = cellfun(@(x) x==0, eventlatency_EMG, ...
        'UniformOutput',false);
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
        FlxS_i  = find(strcmpi(evtypes,'FlxS') & evlat==0, 1);
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