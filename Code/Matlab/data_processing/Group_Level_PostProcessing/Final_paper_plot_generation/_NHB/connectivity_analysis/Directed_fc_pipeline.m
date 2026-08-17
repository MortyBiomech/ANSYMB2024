%% =========================================================================
%  DIRECTED FUNCTIONAL CONNECTIVITY PIPELINE (PDC/DTF)
%  Knee-Tracking MoBI EEG Study — Single Subject Analysis
%
%  Pipeline overview:
%   1. Load subject data and extract IC cluster memberships
%   2. Extract IC activations for 6 clusters
%   3. Timewarp epochs to normalized movement cycle
%   4. Fit MVAR model (with BIC model order selection)
%   5. Compute PDC & DTF per frequency band
%   6. Epoch-averaged connectivity matrix
%   7. Time-resolved connectivity across normalized cycle (sliding window)
%   8. Phase-randomized surrogate testing for significance
%   9. Visualization
%
%  Requirements:
%   - EEGLAB (with STUDY loaded or individual .set file)
%   - BSMART toolbox  OR  custom MVAR functions (included below)
%   - Signal Processing Toolbox
%
%  Author note: Adjust paths and structure field names to match your setup.
% =========================================================================

clear; clc; close all;

%% =========================================================================
%  SECTION 0 — USER SETTINGS  (edit these)
% =========================================================================

% --- Paths ---
eeglab_path     = '/path/to/eeglab';          % path to EEGLAB root
data_path       = '/path/to/data';            % folder containing .set files
subject_id      = 'sub-18';                   % subject to analyze

% --- EEG parameters ---
srate           = 500;                        % sampling rate (Hz)
epoch_tmin_ms   = -500;                       % epoch start (ms)
epoch_tmax_ms   = 3500;                       % epoch end   (ms)
epoch_tmin_s    = epoch_tmin_ms / 1000;
epoch_tmax_s    = epoch_tmax_ms / 1000;
n_samples       = (epoch_tmax_ms - epoch_tmin_ms) / 1000 * srate; % = 2000

% --- Condition labels (as stored in EEG.epoch event types) ---
cond_values     = {1, 3, 6};                  % spring stiffness levels
cond_labels     = {'Low','Medium','High'};

% --- Event names within each epoch ---
event_Flex_Start = 'Flexion_Start';           % time-locking event (t=0)
event_Flex_End   = 'Flexion_End';
event_Ext_End    = 'Extension_End';

% --- Cluster definition ---
%  cluster_ICs: struct with fields named after clusters.
%  Each field is a struct with subject IDs as fields.
%  Value = IC index (integer) or [] if subject has no IC in that cluster.
%
%  Example structure (replace with your actual cluster struct):
%
%  cluster_ICs.L_PrimMotor.sub18  = 3;
%  cluster_ICs.R_PrimMotor.sub18  = 5;
%  cluster_ICs.L_PreMotor.sub18   = [];   % no IC for this subject
%  etc.
%
%  --> Adjust 'cluster_names' and how you index your struct below.

cluster_names = {'L_PrimMotor', 'R_PrimMotor', ...
                 'L_PreMotor',  'R_PreMotor',  ...
                 'L_ParOcc',    'R_ParOcc'};

% Subject key for indexing your cluster struct (adjust naming convention)
subj_key = 'sub18';   % field name as it appears in your cluster struct

% --- MVAR / PDC settings ---
max_model_order = 20;          % maximum model order to test for BIC
freq_bands = struct(...
    'Delta',   [1   4], ...
    'Theta',   [4   8], ...
    'Alpha',   [8  13], ...
    'Beta',    [13 30], ...
    'Gamma',   [30 45]);
freq_resolution = 500;         % number of frequency bins (0 to srate/2)
n_surrogates    = 200;         % number of phase-randomized surrogates
alpha_level     = 0.05;        % significance threshold

% --- Timewarping settings ---
n_norm_samples  = 200;         % normalized cycle length (time points)
% Landmark positions in normalized cycle (will be computed from data median)
% [1, mid1, mid2, n_norm_samples] = [Flex_Start, Flex_End, Ext_End, epoch_end]

% --- Sliding window settings (for time-resolved PDC) ---
win_samples     = 100;         % window length in normalized samples (~50% cycle)
win_step        = 10;          % step size in normalized samples

%% =========================================================================
%  SECTION 1 — SETUP & LOAD DATA
% =========================================================================

% Add EEGLAB to path
addpath(eeglab_path);
eeglab nogui;

% Load subject EEG file
% Assumes one .set file per subject containing all conditions as epochs
set_file = fullfile(data_path, subject_id, [subject_id '_preprocessed_ICA.set']);
fprintf('Loading %s ...\n', set_file);
EEG = pop_loadset(set_file);

fprintf('Loaded: %d ICs, %d epochs, %d samples/epoch, srate=%d Hz\n', ...
    EEG.nbchan, EEG.trials, EEG.pnts, EEG.srate);

%% =========================================================================
%  SECTION 2 — EXTRACT IC INDICES FOR THIS SUBJECT'S CLUSTERS
% =========================================================================

% cluster_ICs should already be loaded in your workspace.
% If it's saved in a .mat file, load it here:
% load(fullfile(data_path, 'cluster_ICs.mat'));  % uncomment and adjust

n_clusters = length(cluster_names);
ic_indices  = nan(1, n_clusters);  % IC index for each cluster (NaN = absent)

for c = 1:n_clusters
    cname = cluster_names{c};
    if isfield(cluster_ICs, cname) && isfield(cluster_ICs.(cname), subj_key)
        ic_val = cluster_ICs.(cname).(subj_key);
        if ~isempty(ic_val)
            ic_indices(c) = ic_val;
        end
    end
end

% Report cluster coverage for this subject
fprintf('\n--- Cluster IC mapping for %s ---\n', subject_id);
valid_clusters = [];
for c = 1:n_clusters
    if ~isnan(ic_indices(c))
        fprintf('  %-25s  IC %d\n', cluster_names{c}, ic_indices(c));
        valid_clusters(end+1) = c; %#ok<AGROW>
    else
        fprintf('  %-25s  [ABSENT]\n', cluster_names{c});
    end
end

% Use only clusters where this subject has an IC
active_cluster_names = cluster_names(valid_clusters);
active_ic_indices    = ic_indices(valid_clusters);
n_active = length(active_ic_indices);

fprintf('\n%d active clusters for analysis.\n', n_active);
if n_active < 2
    error('Need at least 2 clusters for connectivity. Subject has too few ICs.');
end

%% =========================================================================
%  SECTION 3 — EXTRACT EPOCHS PER CONDITION & GET TIMEWARP LATENCIES
% =========================================================================

% EEG.timewarp.latency: [n_epochs x n_landmarks] matrix (in ms from epoch start,
% but stored as sample indices or ms — check your EEGLAB version)
% EEGLAB stores EEG.timewarp.latency as sample indices (1-based) by default.
% We convert to ms below.

% Time vector for the epoch
times_ms = linspace(epoch_tmin_ms, epoch_tmax_ms, n_samples);

% Find condition label for each epoch
epoch_cond = nan(1, EEG.trials);
for ep = 1:EEG.trials
    events_in_epoch = EEG.epoch(ep).eventtype;
    if ~iscell(events_in_epoch)
        events_in_epoch = {events_in_epoch};
    end
    latencies_in_epoch = EEG.epoch(ep).eventlatency;
    if ~iscell(latencies_in_epoch)
        latencies_in_epoch = num2cell(latencies_in_epoch);
    end

    % Find the time-locking event (Flexion_Start at t=0)
    % and extract condition value from its associated field
    for ev = 1:length(events_in_epoch)
        etype = events_in_epoch{ev};
        if ischar(etype) && contains(etype, event_Flex_Start)
            % Condition is encoded in event type string or a separate field
            % Option A: condition encoded in EEG.epoch(ep).eventcondition
            % Option B: condition encoded as numeric in event type
            % --> Adjust based on your actual data structure:
            if isfield(EEG.epoch(ep), 'eventcondition')
                cval = EEG.epoch(ep).eventcondition{ev};
                if iscell(cval), cval = cval{1}; end
                epoch_cond(ep) = str2double(num2str(cval));
            end
            break;
        end
    end
end

% Separate epoch indices per condition
cond_epoch_idx = cell(1, 3);
for ci = 1:3
    cond_epoch_idx{ci} = find(epoch_cond == cond_values{ci});
    fprintf('Condition %s (val=%d): %d epochs\n', ...
        cond_labels{ci}, cond_values{ci}, length(cond_epoch_idx{ci}));
end

%% =========================================================================
%  SECTION 4 — TIMEWARP EPOCHS TO NORMALIZED MOVEMENT CYCLE
% =========================================================================
% EEG.timewarp.latency: [n_epochs x 3] (columns = Flex_Start, Flex_End, Ext_End)
% Values are in ms relative to epoch start (or sample indices — check below).
%
% IMPORTANT: Verify units of EEG.timewarp.latency in your data.
%   - If in ms:      landmark_ms = EEG.timewarp.latency
%   - If in samples: landmark_ms = (EEG.timewarp.latency - 1) / srate * 1000 + epoch_tmin_ms

% --- Convert timewarp latencies to sample indices within epoch ---
% Here we assume EEG.timewarp.latency is in ms relative to epoch start (0 = epoch start)
% Adjust the formula if your data uses a different convention.

tw_latencies_ms = EEG.timewarp.latency;  % [n_epochs x 3]

% Convert ms latencies → sample indices within the epoch (1-based)
% epoch starts at epoch_tmin_ms, so sample index = (lat_ms - epoch_tmin_ms)/1000*srate + 1
tw_samp = round((tw_latencies_ms - epoch_tmin_ms) / 1000 * srate) + 1;
tw_samp = max(1, min(tw_samp, n_samples));  % clamp to valid range

% Compute median landmark positions (in samples) across all epochs
median_landmarks = round(median(tw_samp, 1));  % [1 x 3]
fprintf('\nMedian landmark positions (samples): Flex_Start=%d, Flex_End=%d, Ext_End=%d\n', ...
    median_landmarks(1), median_landmarks(2), median_landmarks(3));

% Define normalized time grid
% Landmarks map to fixed positions in [1 ... n_norm_samples]
% We place them proportionally based on their median positions
norm_lm1 = 1;
norm_lm2 = round((median_landmarks(2) - median_landmarks(1)) / ...
           (median_landmarks(3) - median_landmarks(1)) * (n_norm_samples - 1)) + 1;
norm_lm3 = n_norm_samples;

norm_times = linspace(0, 100, n_norm_samples);  % normalized cycle (0-100%)
fprintf('Normalized landmarks: Flex_Start=%d%%, Flex_End=%d%%, Ext_End=%d%%\n', ...
    norm_times(norm_lm1), round(norm_times(norm_lm2)), norm_times(norm_lm3));

% --- Timewarp function ---
% For each epoch: piecewise linear interpolation from original sample grid
% to normalized grid, using the 3 landmarks as anchors.

function warped = timewarp_epoch(epoch_data, orig_landmarks, norm_lm, n_norm)
    % epoch_data:    [n_channels x n_samples]
    % orig_landmarks: [1 x 3] sample indices in original epoch
    % norm_lm:       [1 x 3] sample indices in normalized epoch
    % n_norm:        scalar, number of output samples

    % Build piecewise mapping: original sample → normalized sample
    orig_pts = [1, orig_landmarks, size(epoch_data, 2)];
    norm_pts = [1, norm_lm,        n_norm];

    % Remove duplicates (in case landmarks coincide)
    [orig_pts, ui] = unique(orig_pts);
    norm_pts = norm_pts(ui);

    % Interpolate each channel
    orig_grid = 1:size(epoch_data, 2);
    norm_grid = interp1(orig_pts, norm_pts, orig_grid, 'linear', 'extrap');
    norm_grid = max(1, min(norm_grid, n_norm));

    new_grid = 1:n_norm;
    warped = zeros(size(epoch_data, 1), n_norm);
    for ch = 1:size(epoch_data, 1)
        warped(ch, :) = interp1(norm_grid, epoch_data(ch, :), new_grid, 'linear', 'extrap');
    end
end

% --- Apply timewarping per condition ---
% Output: warped_data{cond}  = [n_active_clusters x n_norm_samples x n_epochs]

warped_data = cell(1, 3);
for ci = 1:3
    ep_idx = cond_epoch_idx{ci};
    n_ep   = length(ep_idx);
    warped_data{ci} = zeros(n_active, n_norm_samples, n_ep);

    for ep = 1:n_ep
        epoch_idx = ep_idx(ep);
        % Extract IC activations for active clusters [n_active x n_samples]
        raw_epoch = squeeze(EEG.data(active_ic_indices, :, epoch_idx));

        % Get this epoch's landmark positions
        orig_lm = tw_samp(epoch_idx, :);  % [1 x 3]

        % Timewarp
        warped_data{ci}(:, :, ep) = timewarp_epoch(raw_epoch, orig_lm, ...
            [norm_lm1, norm_lm2, norm_lm3], n_norm_samples);
    end
    fprintf('Timewarped %d epochs for condition %s\n', n_ep, cond_labels{ci});
end

%% =========================================================================
%  SECTION 5 — MVAR MODEL ORDER SELECTION (BIC)
% =========================================================================
% We fit the MVAR on all conditions concatenated to find a single optimal
% model order, then use it consistently across conditions.

fprintf('\n--- MVAR Model Order Selection (BIC) ---\n');

% Concatenate all conditions for model order selection
all_warped = cat(3, warped_data{1}, warped_data{2}, warped_data{3});

% Reshape: concatenate epochs along time axis → [n_active x (n_norm * n_epochs)]
% Each epoch is treated as a segment; we do NOT concatenate across epoch boundaries
% (MVAR is fitted per-epoch and coefficients are averaged — Schlögl method)

% For model order selection: fit on a random subset of 50 epochs
rng(42);
n_test_epochs = min(50, size(all_warped, 3));
test_idx      = randperm(size(all_warped, 3), n_test_epochs);

bic_values = nan(1, max_model_order);
for p = 1:max_model_order
    bic_p = 0;
    for ep = 1:n_test_epochs
        seg = all_warped(:, :, test_idx(ep));  % [n_active x n_norm_samples]
        try
            [~, ~, bic_ep] = mvar_fit_bic(seg, p);
            bic_p = bic_p + bic_ep;
        catch
            bic_p = bic_p + inf;
        end
    end
    bic_values(p) = bic_p / n_test_epochs;
    fprintf('  Model order p=%2d  BIC=%.4f\n', p, bic_values(p));
end

[~, optimal_order] = min(bic_values);
fprintf('\nOptimal MVAR model order: p = %d\n', optimal_order);

% Plot BIC curve
figure('Name','MVAR Model Order Selection','NumberTitle','off');
plot(1:max_model_order, bic_values, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xline(optimal_order, 'r--', 'LineWidth', 2, 'Label', sprintf('p=%d', optimal_order));
xlabel('Model Order p'); ylabel('BIC (averaged across epochs)');
title('MVAR Model Order Selection via BIC');
grid on;

%% =========================================================================
%  SECTION 6 — COMPUTE PDC & DTF (EPOCH-AVERAGED CONNECTIVITY MATRIX)
% =========================================================================
% Strategy: fit MVAR per epoch, average AR coefficients across epochs,
% then compute PDC/DTF from the averaged coefficients.
% This is the Schlögl multi-trial approach — valid for quasi-stationary epochs.

p = optimal_order;
n_freqs = freq_resolution;
freqs   = linspace(0, srate/2, n_freqs);  % frequency axis

% Output containers
%   PDC_avg{cond}[from, to, freq]  — epoch-averaged PDC matrix
%   DTF_avg{cond}[from, to, freq]  — epoch-averaged DTF matrix
%   PDC_band{cond}[from, to, band] — band-averaged PDC

PDC_avg  = cell(1, 3);
DTF_avg  = cell(1, 3);
PDC_band = cell(1, 3);
DTF_band = cell(1, 3);

band_names = fieldnames(freq_bands);
n_bands    = length(band_names);

fprintf('\n--- Computing epoch-averaged PDC/DTF ---\n');

for ci = 1:3
    n_ep = size(warped_data{ci}, 3);
    k    = n_active;

    % Accumulate AR coefficient matrices across epochs
    AR_sum = zeros(k, k*p);  % [k x k*p] concatenated AR matrices
    SIG_sum = zeros(k, k);

    valid_ep = 0;
    for ep = 1:n_ep
        seg = warped_data{ci}(:, :, ep);  % [k x n_norm_samples]
        try
            [AR_ep, SIG_ep] = mvar_fit(seg, p);
            AR_sum  = AR_sum  + AR_ep;
            SIG_sum = SIG_sum + SIG_ep;
            valid_ep = valid_ep + 1;
        catch
            % Skip epochs where MVAR fails (rank-deficient)
        end
    end

    AR_avg  = AR_sum  / valid_ep;
    SIG_avg = SIG_sum / valid_ep;

    fprintf('  Condition %s: fitted MVAR on %d/%d epochs\n', ...
        cond_labels{ci}, valid_ep, n_ep);

    % Compute PDC and DTF from averaged AR coefficients
    [PDC_avg{ci}, DTF_avg{ci}] = compute_PDC_DTF(AR_avg, SIG_avg, p, k, freqs, srate);

    % Average over frequency bands
    PDC_band{ci} = zeros(k, k, n_bands);
    DTF_band{ci} = zeros(k, k, n_bands);
    for b = 1:n_bands
        band_range = freq_bands.(band_names{b});
        freq_idx   = freqs >= band_range(1) & freqs <= band_range(2);
        PDC_band{ci}(:,:,b) = mean(PDC_avg{ci}(:,:,freq_idx), 3);
        DTF_band{ci}(:,:,b) = mean(DTF_avg{ci}(:,:,freq_idx), 3);
    end
end

%% =========================================================================
%  SECTION 7 — TIME-RESOLVED PDC (SLIDING WINDOW ACROSS NORMALIZED CYCLE)
% =========================================================================
% For each sliding window position, fit MVAR on the windowed segment
% across all epochs, then compute PDC.
% Output: PDC_timeres{cond}[from, to, freq, time_window]

fprintf('\n--- Computing time-resolved PDC (sliding window) ---\n');

win_starts = 1 : win_step : (n_norm_samples - win_samples + 1);
n_windows  = length(win_starts);
win_centers_pct = (win_starts + win_samples/2 - 1) / n_norm_samples * 100;

PDC_timeres = cell(1, 3);
DTF_timeres = cell(1, 3);

for ci = 1:3
    n_ep = size(warped_data{ci}, 3);
    k    = n_active;
    PDC_timeres{ci} = zeros(k, k, n_freqs, n_windows);
    DTF_timeres{ci} = zeros(k, k, n_freqs, n_windows);

    fprintf('  Condition %s: computing %d windows...\n', cond_labels{ci}, n_windows);
    for w = 1:n_windows
        w_start = win_starts(w);
        w_end   = w_start + win_samples - 1;

        % Fit MVAR on windowed segments, averaged across epochs
        AR_sum  = zeros(k, k*p);
        SIG_sum = zeros(k, k);
        valid_ep = 0;

        for ep = 1:n_ep
            seg = warped_data{ci}(:, w_start:w_end, ep);  % [k x win_samples]
            if size(seg, 2) < k*p + 10
                continue;  % too short for this model order
            end
            try
                [AR_ep, SIG_ep] = mvar_fit(seg, p);
                AR_sum  = AR_sum  + AR_ep;
                SIG_sum = SIG_sum + SIG_ep;
                valid_ep = valid_ep + 1;
            catch
                continue;
            end
        end

        if valid_ep < 10
            continue;  % not enough valid epochs for this window
        end

        AR_w   = AR_sum  / valid_ep;
        SIG_w  = SIG_sum / valid_ep;
        [PDC_w, DTF_w] = compute_PDC_DTF(AR_w, SIG_w, p, k, freqs, srate);
        PDC_timeres{ci}(:,:,:,w) = PDC_w;
        DTF_timeres{ci}(:,:,:,w) = DTF_w;
    end
end

%% =========================================================================
%  SECTION 8 — SURROGATE TESTING (PHASE RANDOMIZATION)
% =========================================================================
% For each condition, generate phase-randomized surrogates of the
% timewarped data, recompute PDC, and build a null distribution.
% Threshold = 95th percentile of surrogate PDC per frequency bin.

fprintf('\n--- Surrogate testing (%d surrogates) ---\n', n_surrogates);

PDC_surr_thresh = cell(1, 3);  % [k x k x n_freqs] — significance threshold

for ci = 1:3
    k    = n_active;
    n_ep = size(warped_data{ci}, 3);
    surr_PDC = zeros(k, k, n_freqs, n_surrogates);

    for s = 1:n_surrogates
        AR_sum  = zeros(k, k*p);
        SIG_sum = zeros(k, k);
        valid_ep = 0;

        for ep = 1:n_ep
            seg = warped_data{ci}(:, :, ep);  % [k x n_norm_samples]
            % Phase randomize: randomize phase of FFT, keep magnitude
            seg_surr = phase_randomize(seg);
            try
                [AR_ep, SIG_ep] = mvar_fit(seg_surr, p);
                AR_sum  = AR_sum  + AR_ep;
                SIG_sum = SIG_sum + SIG_ep;
                valid_ep = valid_ep + 1;
            catch
                continue;
            end
        end

        if valid_ep > 0
            AR_s = AR_sum / valid_ep;
            SIG_s = SIG_sum / valid_ep;
            [PDC_s, ~] = compute_PDC_DTF(AR_s, SIG_s, p, k, freqs, srate);
            surr_PDC(:,:,:,s) = PDC_s;
        end

        if mod(s, 50) == 0
            fprintf('  Condition %s: %d/%d surrogates done\n', ...
                cond_labels{ci}, s, n_surrogates);
        end
    end

    % Threshold at (1-alpha) percentile across surrogates
    PDC_surr_thresh{ci} = prctile(surr_PDC, (1-alpha_level)*100, 4);
    fprintf('  Condition %s: surrogate thresholds computed\n', cond_labels{ci});
end

% Significant PDC mask: PDC > surrogate threshold
PDC_sig = cell(1, 3);
for ci = 1:3
    PDC_sig{ci} = PDC_avg{ci} > PDC_surr_thresh{ci};
end

%% =========================================================================
%  SECTION 9 — VISUALIZATION
% =========================================================================

%% 9A — Epoch-averaged PDC matrices per condition and frequency band
band_to_plot = 'Beta';   % change to any band name
b_idx = find(strcmp(band_names, band_to_plot));

figure('Name', sprintf('PDC Matrices — %s band', band_to_plot), ...
       'Position', [100 100 1400 400], 'NumberTitle', 'off');

clim_max = 0;
for ci = 1:3
    clim_max = max(clim_max, max(PDC_band{ci}(:,:,b_idx), [], 'all'));
end

for ci = 1:3
    subplot(1, 3, ci);
    pdc_mat = PDC_band{ci}(:,:,b_idx);

    % Mask non-significant connections (based on full-spectrum threshold)
    surr_band_thresh = mean(PDC_surr_thresh{ci}(:,:, ...
        freqs >= freq_bands.(band_to_plot)(1) & freqs <= freq_bands.(band_to_plot)(2)), 3);
    sig_mask = pdc_mat > surr_band_thresh;

    imagesc(pdc_mat .* sig_mask);
    colormap(gca, 'hot');
    clim([0 clim_max]);
    colorbar;
    axis square;
    set(gca, 'XTick', 1:n_active, 'XTickLabel', active_cluster_names, ...
             'YTick', 1:n_active, 'YTickLabel', active_cluster_names, ...
             'XTickLabelRotation', 45, 'FontSize', 9);
    title(sprintf('%s (spring=%d)', cond_labels{ci}, cond_values{ci}));
    xlabel('To (target)');
    ylabel('From (source)');
end
sgtitle(sprintf('Epoch-averaged PDC — %s band (significant only, α=%.2f)', ...
    band_to_plot, alpha_level));

%% 9B — Time-resolved PDC for a specific connection across conditions
% Choose a connection of interest, e.g., L_PrimMotor → R_PrimMotor
% Find indices of these clusters in active_cluster_names
from_cluster = 'L_PrimMotor';
to_cluster   = 'R_PrimMotor';
from_idx = find(strcmp(active_cluster_names, from_cluster));
to_idx   = find(strcmp(active_cluster_names, to_cluster));

if ~isempty(from_idx) && ~isempty(to_idx)
    figure('Name', sprintf('Time-resolved PDC: %s → %s', from_cluster, to_cluster), ...
           'Position', [100 560 1400 700], 'NumberTitle', 'off');

    colors = [0.2 0.6 1; 1 0.6 0; 0.8 0.2 0.2];  % blue, orange, red

    for b = 1:n_bands
        bname = band_names{b};
        band_range = freq_bands.(bname);
        f_idx = freqs >= band_range(1) & freqs <= band_range(2);

        subplot(n_bands, 1, b);
        hold on;

        for ci = 1:3
            % Average PDC over frequency band, for this connection, across windows
            pdc_trace = squeeze(mean(PDC_timeres{ci}(from_idx, to_idx, f_idx, :), 3));
            plot(win_centers_pct, pdc_trace, '-', 'Color', colors(ci,:), ...
                'LineWidth', 2, 'DisplayName', cond_labels{ci});
        end

        % Mark landmark events as vertical lines
        xline(norm_times(norm_lm1), 'k--', 'Flex Start', 'LineWidth', 1.5, ...
            'LabelVerticalAlignment', 'bottom');
        xline(norm_times(norm_lm2), 'k:', 'Flex End', 'LineWidth', 1.5, ...
            'LabelVerticalAlignment', 'bottom');
        xline(norm_times(norm_lm3), 'k-.', 'Ext End', 'LineWidth', 1.5, ...
            'LabelVerticalAlignment', 'bottom');

        ylabel('PDC'); title(sprintf('%s band (%d-%d Hz)', bname, band_range(1), band_range(2)));
        if b == 1, legend('Location','best'); end
        if b == n_bands, xlabel('Normalized cycle position (%)'); end
        grid on; xlim([0 100]);
        hold off;
    end
    sgtitle(sprintf('Time-resolved PDC: %s → %s', from_cluster, to_cluster));
else
    warning('Connection %s → %s not found in active clusters.', from_cluster, to_cluster);
end

%% 9C — Full connectivity graph for one condition (Beta band)
figure('Name', 'Connectivity Graph — Beta band', 'NumberTitle', 'off', ...
       'Position', [100 100 800 700]);

ci_plot = 1;  % change to 1/2/3 for Low/Medium/High

% Layout: arrange clusters in approximate anatomical positions
%   L_PrimMotor  R_PrimMotor
%   L_PreMotor   R_PreMotor
%   L_ParOcc     R_ParOcc
layout_pos = [-1 2; 1 2; -1 0; 1 0; -1 -2; 1 -2];
% Trim to active clusters only
active_pos = layout_pos(valid_clusters, :);

pdc_mat  = PDC_band{ci_plot}(:,:,b_idx);
surr_b   = mean(PDC_surr_thresh{ci_plot}(:,:, ...
    freqs >= freq_bands.(band_to_plot)(1) & freqs <= freq_bands.(band_to_plot)(2)), 3);
sig_mask = pdc_mat > surr_b;

hold on;
% Draw arrows for significant connections
for from = 1:n_active
    for to = 1:n_active
        if from ~= to && sig_mask(from, to)
            x_from = active_pos(from, 1);
            y_from = active_pos(from, 2);
            x_to   = active_pos(to,   1);
            y_to   = active_pos(to,   2);
            strength = pdc_mat(from, to);
            lw = 1 + strength * 15;  % scale line width by PDC strength
            quiver(x_from, y_from, x_to-x_from, y_to-y_from, 0, ...
                'Color', [0.8*strength, 0.2, 1-strength], ...
                'LineWidth', lw, 'MaxHeadSize', 0.5);
        end
    end
end

% Draw nodes
for c = 1:n_active
    scatter(active_pos(c,1), active_pos(c,2), 300, 'filled', ...
        'MarkerFaceColor', [0.3 0.6 0.9], 'MarkerEdgeColor', 'k');
    text(active_pos(c,1), active_pos(c,2), active_cluster_names{c}, ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

axis off; axis equal;
title(sprintf('Significant PDC connections — %s — %s band\n(arrow width = PDC strength)', ...
    cond_labels{ci_plot}, band_to_plot));
hold off;

%% =========================================================================
%  SECTION 10 — SAVE RESULTS
% =========================================================================

results.subject_id          = subject_id;
results.cluster_names       = active_cluster_names;
results.ic_indices          = active_ic_indices;
results.cond_labels         = cond_labels;
results.cond_values         = cond_values;
results.freqs               = freqs;
results.freq_bands          = freq_bands;
results.optimal_mvar_order  = optimal_order;
results.PDC_avg             = PDC_avg;
results.DTF_avg             = DTF_avg;
results.PDC_band            = PDC_band;
results.DTF_band            = DTF_band;
results.PDC_timeres         = PDC_timeres;
results.DTF_timeres         = DTF_timeres;
results.PDC_surr_thresh     = PDC_surr_thresh;
results.PDC_sig             = PDC_sig;
results.win_centers_pct     = win_centers_pct;
results.norm_landmark_pct   = [norm_times(norm_lm1), ...
                                norm_times(norm_lm2), ...
                                norm_times(norm_lm3)];
results.warped_data         = warped_data;  % save for inspection
results.bic_values          = bic_values;

save_path = fullfile(data_path, subject_id, [subject_id '_PDC_results.mat']);
save(save_path, 'results', '-v7.3');
fprintf('\nResults saved to: %s\n', save_path);
fprintf('Pipeline complete!\n');

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [AR, SIG] = mvar_fit(data, p)
    % Fit MVAR model of order p to data [k x T] using least squares (Yule-Walker)
    % Returns:
    %   AR  [k x k*p]  — concatenated AR coefficient matrices [A1 A2 ... Ap]
    %   SIG [k x k]    — noise covariance matrix

    [k, T] = size(data);
    if T <= k * p + k
        error('Too few samples for model order p=%d with k=%d channels', p, k);
    end

    % Build regression matrix
    % Y = [data(:, p+1 : T)]   [k x (T-p)]
    % X = [lagged data matrix] [k*p x (T-p)]
    Y = data(:, p+1:T);          % [k x (T-p)]
    n_obs = T - p;

    X = zeros(k*p, n_obs);
    for lag = 1:p
        X((lag-1)*k+1 : lag*k, :) = data(:, p+1-lag : T-lag);
    end

    % Least squares: AR = Y * X' * (X * X')^{-1}
    AR  = Y * X' / (X * X');     % [k x k*p]
    res = Y - AR * X;            % residuals [k x (T-p)]
    SIG = (res * res') / (n_obs - k*p - 1);  % noise covariance
end

function [AR, SIG, bic] = mvar_fit_bic(data, p)
    % Fit MVAR and compute BIC
    [k, T] = size(data);
    [AR, SIG] = mvar_fit(data, p);
    n_obs = T - p;
    n_params = k^2 * p;

    % Log-likelihood of multivariate Gaussian
    log_det_SIG = log(det(SIG));
    log_lik = -n_obs/2 * (k * log(2*pi) + log_det_SIG + k);
    bic = -2 * log_lik + n_params * log(n_obs);
end

function [PDC, DTF] = compute_PDC_DTF(AR, SIG, p, k, freqs, srate)
    % Compute PDC and DTF from MVAR coefficients
    %
    % AR:    [k x k*p]  concatenated AR matrices
    % SIG:   [k x k]    noise covariance
    % p:     model order
    % k:     number of channels
    % freqs: frequency vector (Hz)
    % srate: sampling rate (Hz)
    %
    % PDC[i,j,f] = direct causal influence FROM j TO i at frequency f
    % DTF[i,j,f] = total causal influence FROM j TO i at frequency f

    n_freqs = length(freqs);
    PDC = zeros(k, k, n_freqs);
    DTF = zeros(k, k, n_freqs);

    % Reshape AR into cell array of matrices {A1, A2, ..., Ap}
    A = cell(1, p);
    for lag = 1:p
        A{lag} = AR(:, (lag-1)*k+1 : lag*k);  % [k x k]
    end

    for fi = 1:n_freqs
        f  = freqs(fi);
        z  = exp(-1j * 2 * pi * f / srate);  % unit circle

        % Spectral matrix of AR coefficients: A(f) = I - sum_lag A_lag * z^-lag
        Af = eye(k);
        for lag = 1:p
            Af = Af - A{lag} * z^(-lag);
        end

        % Transfer matrix H(f) = A(f)^{-1}
        H  = inv(Af);  % [k x k]

        % --- PDC (normalized by column norm of A(f)) ---
        % PDC_{ij}(f) = |A_{ij}(f)|^2 / sum_m |A_{mj}(f)|^2
        % Interpretation: fraction of j's output going directly to i
        for j = 1:k
            col_norm_sq = sum(abs(Af(:,j)).^2);
            for i = 1:k
                if i ~= j
                    PDC(i,j,fi) = abs(Af(i,j))^2 / col_norm_sq;
                end
            end
        end

        % --- DTF (normalized by row norm of H(f)) ---
        % DTF_{ij}(f) = |H_{ij}(f)|^2 / sum_m |H_{im}(f)|^2
        % Interpretation: fraction of i's input coming from j (all paths)
        for i = 1:k
            row_norm_sq = sum(abs(H(i,:)).^2);
            for j = 1:k
                if i ~= j
                    DTF(i,j,fi) = abs(H(i,j))^2 / row_norm_sq;
                end
            end
        end
    end
end

function surr = phase_randomize(data)
    % Phase randomization: preserve spectral amplitude, randomize phase
    % data: [k x T]
    [k, T] = size(data);
    surr   = zeros(k, T);
    for ch = 1:k
        xf  = fft(data(ch,:));
        mag = abs(xf);
        pha = angle(xf);

        % Random phase shift (preserving conjugate symmetry for real signal)
        n_half = floor(T/2);
        rand_phase = 2 * pi * rand(1, n_half - 1);  % random phases

        % Build new phase vector
        new_pha = pha;
        new_pha(2 : n_half) = rand_phase;
        new_pha(end : -1 : end - n_half + 2) = -rand_phase;

        surr(ch,:) = real(ifft(mag .* exp(1j * new_pha)));
    end
end