% FAUVET2019_SIMULATION
%
%   Reproduces the "Simulated data" section of:
%
%     Fauvet et al. (2019). "A novel method to generalize time-frequency
%     coherence analysis between EEG or EMG signals during repetitive
%     trials with high intra-subject variability in duration."
%     9th Int. IEEE EMBS Conf. Neural Engineering, pp. 437-440.
%
%   Exact parameters from paper text and Fig. 1 caption:
%     - 30 paired trials, zero mean, -5 dB SNR, Fs = 1000 Hz
%     - SOI: 10 Hz + 30 Hz sinusoids (fully coherent between X and Y)
%     - SOI length range: 1353 to 4789 samples  (Fig. 1 caption)
%     - Pre-SOI flank: uniform random [1.0, 2.1] s of Gaussian white noise
%     - Post-SOI flank: Gaussian white noise, total signal length must
%       not exceed max(SOI_dur) + 2.1 s
%     - Aligned: 1 s before SOI onset + SOI + noise pad
%       Reference length = longest_SOI + 1000 samples = 5789 samples
%     - Time-normalised: upsample aligned signal so SOI = 4789 samples
%       Fixed output length = 7099 points (Fig. 1 caption)
%       Effective srate ranges from 1000 to 3817 Hz (Fig. 1 caption)

clear; clc; close all;
rng(42);   % reproducibility

% =========================================================================
% PARAMETERS  (directly from paper and Fig. 1 caption)
% =========================================================================
Fs          = 1000;        % Hz
nTrials     = 30;
SNR_dB      = -5;
snr_lin     = 10^(SNR_dB/10);   % signal_power / noise_power ≈ 0.3162
noise_std   = 1 / sqrt(snr_lin);   % noise amplitude ≈ 1.778 (signal amp = 1)

f1          = 10;          % Hz — alpha
f2          = 30;          % Hz — beta

% SOI range from Fig. 1 caption
soi_min_n   = 1353;        % samples
soi_max_n   = 4789;        % samples  (l2 — the reference)

% Flanks (paper)
pre_min_n   = 1000;        % 1.0 s
pre_max_n   = 2100;        % 2.1 s
total_max_n = soi_max_n + 2100;   % = 6889 samples

% Aligned reference (paper: "longest SOI plus one second" = 4789 + 1000)
ref_aligned_n = soi_max_n + Fs;   % = 5789 samples

% Time-normalised fixed output length (Fig. 1 caption)
ref_norm_n  = 7099;

% =========================================================================
% STEP 1 — Draw SOI and flank durations
% =========================================================================
% Paper says "normal distribution" but only the range [1353, 4789] is
% given numerically. We draw uniformly within this range so every trial
% stays inside the documented bounds.
soi_n = randi([soi_min_n, soi_max_n], nTrials, 1);   % [nTrials x 1]

% Pre-SOI flank: uniform [1000, 2100] samples
pre_n = randi([pre_min_n, pre_max_n], nTrials, 1);

% Post-SOI: random, total must not exceed total_max_n
avail_post = total_max_n - pre_n - soi_n;
avail_post = max(100, avail_post);                     % minimum 100 ms
post_n = round(rand(nTrials, 1) .* avail_post);

fprintf('=== Fauvet 2019 Simulation ===\n');
fprintf('SOI length  : min=%d  max=%d  mean=%.0f samples\n', ...
    min(soi_n), max(soi_n), mean(soi_n));
fprintf('Pre-flank   : min=%d  max=%d samples\n', min(pre_n), max(pre_n));
fprintf('Total length: min=%d  max=%d samples\n', ...
    min(pre_n+soi_n+post_n), max(pre_n+soi_n+post_n));

% =========================================================================
% STEP 2 — Build RAW signals  (Fig. 1, left column)
%
%   Structure: [pre-flank noise | SOI: sinusoids + noise | post-flank noise]
%   X and Y share identical sinusoidal SOI components (fully coherent)
%   but have independent noise realizations everywhere.
% =========================================================================
X_raw = cell(nTrials, 1);
Y_raw = cell(nTrials, 1);

for t = 1:nTrials
    n_tot = pre_n(t) + soi_n(t) + post_n(t);

    % Independent noise for X and Y
    nx = noise_std * randn(n_tot, 1);
    ny = noise_std * randn(n_tot, 1);

    % Deterministic SOI signal: same for both channels → perfect coherence
    soi_time = (0 : soi_n(t)-1)' / Fs;
    soi_sig  = sin(2*pi*f1*soi_time) + sin(2*pi*f2*soi_time);

    % Add SOI signal into the central segment
    soi_idx = pre_n(t)+1 : pre_n(t)+soi_n(t);
    nx(soi_idx) = nx(soi_idx) + soi_sig;
    ny(soi_idx) = ny(soi_idx) + soi_sig;

    % Zero-mean
    X_raw{t} = nx - mean(nx);
    Y_raw{t} = ny - mean(ny);
end
fprintf('Raw signals built.\n');

% =========================================================================
% STEP 3 — Build ALIGNED signals  (Fig. 1, middle column)
%
%   Paper: "settling the raw signals 1 s prior to the segment of interest"
%   → Extract exactly 1000 samples (1 s) before SOI onset from raw signal.
%   → Append SOI segment.
%   → Pad end with Gaussian white noise to reach ref_aligned_n = 5789.
%
%   In the aligned epoch:
%     samples  1–1000     : 1 s pre-SOI noise
%     samples  1001–(1000+soi_n(t)) : SOI
%     samples  after SOI  : white noise pad
% =========================================================================
X_aligned = zeros(ref_aligned_n, nTrials);
Y_aligned = zeros(ref_aligned_n, nTrials);

soi_onset_aligned_n  = Fs;        % SOI starts at sample 1001 (0-indexed: 1000)
soi_onsets_ms_al     = Fs;        % = 1000 ms (fixed for all trials)
soi_offsets_ms_al    = soi_onsets_ms_al + soi_n / Fs * 1000;  % [nTrials x 1]

for t = 1:nTrials
    rx = X_raw{t};
    ry = Y_raw{t};

    % 1 s pre-SOI: take last 1000 samples before the SOI in raw signal
    pre_start = pre_n(t) - Fs + 1;
    if pre_start < 1
        % pre-flank shorter than 1 s → left-pad with fresh noise
        pad = Fs - pre_n(t);
        pre_x = [noise_std*randn(pad,1); rx(1:pre_n(t))];
        pre_y = [noise_std*randn(pad,1); ry(1:pre_n(t))];
    else
        pre_x = rx(pre_start : pre_n(t));
        pre_y = ry(pre_start : pre_n(t));
    end
    % pre_x / pre_y are now exactly Fs = 1000 samples

    % SOI
    soi_x = rx(pre_n(t)+1 : pre_n(t)+soi_n(t));
    soi_y = ry(pre_n(t)+1 : pre_n(t)+soi_n(t));

    % Concatenate
    seg_x = [pre_x; soi_x];
    seg_y = [pre_y; soi_y];

    % Pad or truncate to ref_aligned_n
    cur = length(seg_x);
    if cur < ref_aligned_n
        pad = ref_aligned_n - cur;
        seg_x = [seg_x; noise_std*randn(pad,1)];
        seg_y = [seg_y; noise_std*randn(pad,1)];
    else
        seg_x = seg_x(1:ref_aligned_n);
        seg_y = seg_y(1:ref_aligned_n);
    end

    X_aligned(:,t) = seg_x;
    Y_aligned(:,t) = seg_y;
end
fprintf('Aligned signals built. ref_aligned_n = %d samples (%.0f ms).\n', ...
    ref_aligned_n, ref_aligned_n/Fs*1000);

% =========================================================================
% STEP 4 — Build TIME-NORMALISED signals  (Fig. 1, right column)
%
%   Paper: upsample each aligned signal so SOI becomes l2=4789 samples.
%   Upsampling ratio for trial t: p/q ≈ l2 / soi_n(t)
%   The full aligned epoch is upsampled at this ratio.
%   Output is truncated/padded to ref_norm_n = 7099 samples.
%
%   Time axis becomes % of movement:
%     0%   = SOI onset  (was at sample 1001 in aligned, now at
%             round(1000 * l2/soi_n(t)) + 1 in normalised)
%     100% = SOI offset (onset + l2 samples)
% =========================================================================
l2 = soi_max_n;   % = 4789

X_norm      = zeros(ref_norm_n, nTrials);
Y_norm      = zeros(ref_norm_n, nTrials);
fs_eff_norm = zeros(nTrials, 1);

for t = 1:nTrials
    l1 = soi_n(t);
    [p, q] = rat(l2 / l1, 1e-4);
    fs_eff_norm(t) = Fs * p / q;

    xu = resample(double(X_aligned(:,t)), p, q);
    yu = resample(double(Y_aligned(:,t)), p, q);

    len = length(xu);
    X_norm(1:min(len,ref_norm_n), t) = xu(1:min(len,ref_norm_n));
    Y_norm(1:min(len,ref_norm_n), t) = yu(1:min(len,ref_norm_n));
end
fprintf('Time-normalised signals built. ref_norm_n = %d samples.\n', ref_norm_n);
fprintf('Effective srate: min=%.0f Hz  max=%.0f Hz\n', ...
    min(fs_eff_norm), max(fs_eff_norm));

% =========================================================================
% STEP 5 — PLOT: reproduce Fig. 1
% =========================================================================
trials_to_show = [1, 2, 3, round(nTrials/2), nTrials];
row_labels     = {'Trial 1','Trial 2','Trial 3', ...
                  sprintf('Trial %d',round(nTrials/2)),'Trial 30'};
nShow = length(trials_to_show);

% x-axis for time-normalised: compute onset position in normalised epoch
% SOI onset in aligned = 1000 samples; after upsampling by l2/l1:
onset_norm_n = zeros(nTrials,1);
for t = 1:nTrials
    onset_norm_n(t) = round(Fs * l2 / soi_n(t));   % upsampled pre-flank
end

figure('Name','Fauvet 2019 — Fig. 1 reproduction', ...
    'Position',[30 30 1500 800], 'Color','w');

for s = 1:nShow
    t = trials_to_show(s);

    % ---- Raw ----
    ax = subplot(nShow, 3, (s-1)*3 + 1);
    t_raw_ms = (0:length(X_raw{t})-1) / Fs * 1000;
    plot(t_raw_ms, X_raw{t}, 'k', 'LineWidth', 0.4);
    xlim([0 max(cellfun(@length,X_raw))/Fs*1000]);
    ylim([-8 8]);  grid off; box off;
    hold on;
    xline(pre_n(t)/Fs*1000,            'r--', 'LineWidth',1);
    xline((pre_n(t)+soi_n(t))/Fs*1000, 'r--', 'LineWidth',1);
    hold off;
    if s==1,    title('raw signals','FontWeight','bold'); end
    if s==nShow, xlabel('Time (ms)'); end
    ylabel(row_labels{s},'FontWeight','bold','FontSize',9);

    % ---- Aligned ----
    subplot(nShow, 3, (s-1)*3 + 2);
    t_al_ms = (0:ref_aligned_n-1) / Fs * 1000;
    plot(t_al_ms, X_aligned(:,t), 'k', 'LineWidth', 0.4);
    xlim([0 (ref_aligned_n-1)/Fs*1000]);
    ylim([-8 8]);  grid off; box off;
    hold on;
    xline(soi_onsets_ms_al,       'r--', 'LineWidth',1);
    xline(soi_offsets_ms_al(t),   'r--', 'LineWidth',1);
    hold off;
    if s==1,    title('aligned signals','FontWeight','bold'); end
    if s==nShow, xlabel('Movement time (ms)'); end

    % ---- Time-normalised ----
    subplot(nShow, 3, (s-1)*3 + 3);
    on  = onset_norm_n(t);
    pct = ((0:ref_norm_n-1) - on) / l2 * 100;
    plot(pct, X_norm(:,t), 'k', 'LineWidth', 0.4);
    xlim([pct(1) pct(end)]);
    ylim([-8 8]);  grid off; box off;
    hold on;
    xline(0,   'r--', 'LineWidth',1);
    xline(100, 'r--', 'LineWidth',1);
    hold off;
    if s==1,    title('time-normalized signals','FontWeight','bold'); end
    if s==nShow, xlabel('Percent of movement realization'); end
end

sgtitle('Fauvet et al. (2019) — Simulation: Fig. 1 reproduction', ...
    'FontWeight','bold','FontSize',13);

% =========================================================================
% STEP 6 — Prepare inputs for my_newcrossf_timenorm
% =========================================================================
% Concatenate aligned trials into single row vectors
x_al = reshape(X_aligned, 1, []);
y_al = reshape(Y_aligned, 1, []);

tlimits_al = [0,  (ref_aligned_n-1)/Fs*1000];   % [0, 5788] ms

% SOI onset/offset vectors (ms) in the aligned epoch
soi_onsets_ms_vec  = repmat(double(soi_onsets_ms_al), nTrials, 1);  % all 1000 ms
soi_offsets_ms_vec = soi_offsets_ms_al;                              % [nTrials x 1]

fprintf('\n=== Ready for coherence computation ===\n');
fprintf('Aligned input:\n');
fprintf('  x_al / y_al : 1 x %d\n', length(x_al));
fprintf('  frame       : %d\n', ref_aligned_n);
fprintf('  tlimits     : [%.0f  %.0f] ms\n', tlimits_al(1), tlimits_al(2));
fprintf('  soi_onsets  : all %.0f ms\n', soi_onsets_ms_al);
fprintf('  soi_offsets : min=%.0f ms   max=%.0f ms\n', ...
    min(soi_offsets_ms_vec), max(soi_offsets_ms_vec));

% =========================================================================
% STEP 7 — Coherence computation  (uncomment to run)
% =========================================================================
% Requires my_newcrossf_timenorm.m on the MATLAB path.
%
% [coherres, mbase, timesout_pct, freqs, ~, Rangle] = ...
%     my_newcrossf_timenorm( ...
%         x_al, y_al,                  ...
%         ref_aligned_n,               ...   % frame
%         tlimits_al,                  ...   % [0, 5788] ms
%         Fs,                          ...   % 1000 Hz
%         [3 0.8],                     ...   % Morlet cycles
%         soi_onsets_ms_vec,           ...   % all 1000 ms
%         soi_offsets_ms_vec,          ...   % per-trial ExtE latency
%         'freqs',    [1 50],          ...
%         'nfreqs',   100,             ...
%         'timesout', 200,             ...
%         'plotamp',  'off',           ...
%         'plotphase','off');
%
% figure;
% imagesc(timesout_pct, freqs, coherres);
% axis xy; colorbar; clim([0 1]);
% xlabel('% of movement realization');
% ylabel('Frequency (Hz)');
% title('Time-frequency coherence — time-normalised (Fauvet 2019)');
% hold on;
% xline(0,   'w--', 'LineWidth', 1.5);   % SOI onset  = 0%
% xline(100, 'w--', 'LineWidth', 1.5);   % SOI offset = 100%
% yline(f1, 'w-',  'LineWidth', 1);      % 10 Hz marker
% yline(f2, 'w-',  'LineWidth', 1);      % 30 Hz marker
% hold off;