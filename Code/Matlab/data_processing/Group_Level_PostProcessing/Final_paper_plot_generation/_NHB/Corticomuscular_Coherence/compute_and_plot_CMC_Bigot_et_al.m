function [coh, sig_mask, lambda, fig_handle] = ...
    compute_and_plot_CMC_Bigot_et_al(TF_X, TF_Y, percent_grid, freqs, ...
                         X_epochs, Y_epochs, timewarp_indices, ...
                         params, opts)
% COMPUTE_AND_PLOT_CMC
%   Magnitude-squared coherence (MSC) between an EEG TF representation
%   and one or more EMG TF representations, with significance masking
%   based on the wavelet cross-spectrum threshold of Bigot et al. (2011),
%   and a 2x2 plot.
%
%   Pipeline (matches Fauvet/Cremoux/Amarantini group, 5-step procedure):
%     1. Compute per-pixel MSC from TF_X and TF_Y (NaN-aware sums).
%     2. Compute |S_xy(f,t)| (magnitude of trial-averaged cross-spectrum).
%     3. Re-upsample the original time-domain epochs (same as in
%        compute_TF_timenorm) and isolate the per-trial SOI segments
%        (all of equal length max_SOI = l2 samples).
%     4. From these matched-length SOI signals, compute the empirical
%        covariance per channel and its largest eigenvalue rho_x, rho_y.
%     5. Compute the Bigot et al. (2011) data-driven threshold
%        lambda_alpha and mask MSC pixels where |S_xy| < lambda_alpha.
%
% Inputs:
%   TF_X             = [nFreqs x nTimes x nEpochs] complex EEG TF
%                      (output of compute_TF_timenorm)
%   TF_Y             = [nFreqs x nTimes x nEpochs x nMuscles] complex EMG TF
%                      (or [nFreqs x nTimes x nEpochs] for one muscle)
%   percent_grid     = [1 x nTimes] cycle-percent time axis
%   freqs            = [1 x nFreqs] frequency vector (Hz)
%   X_epochs         = [nSamples x nEpochs] original EEG epochs
%   Y_epochs         = [nMuscles x nSamples x nEpochs] OR [nSamples x nEpochs]
%   timewarp_indices = [nEpochs x nEvents] event sample indices
%   params (struct):
%       .srate              (Hz, original sampling rate)
%       .soi_event_indices  (default [1 size(timewarp_indices,2)])
%       .resample_tol       (default 1e-8)
%       .verbose            (default true)
%   opts (struct, optional):
%       .alpha          - significance level (default 0.05)
%       .muscle_labels  - default {'VL','RF','GM','BF'} or auto
%       .freq_ticks     - default [4 8 14 30 60]
%       .freq_lim       - default [freqs(1) freqs(end)]
%       .title          - default 'CMC (Bigot threshold)'
%       .colormap       - default 'hot'
%       .plot           - default true
%       .clim           - default [0  prctile(coh(:),90)]
%
% Outputs:
%   coh        = [nFreqs x nTimes x nMuscles] MSC, with non-significant
%                pixels set to NaN
%   sig_mask   = [nFreqs x nTimes x nMuscles] logical mask (true where
%                |S_xy| > lambda_alpha)
%   lambda     = [1 x nMuscles] Bigot threshold per muscle
%   fig_handle = handle to the figure (empty if opts.plot = false)

% --- defaults ---
if nargin < 9, opts = struct(); end
if ~isfield(opts, 'alpha'),       opts.alpha       = 0.05; end
if ~isfield(opts, 'freq_ticks'),  opts.freq_ticks  = [4 8 14 30 60]; end
if ~isfield(opts, 'freq_lim'),    opts.freq_lim    = [freqs(1) freqs(end)]; end
if ~isfield(opts, 'title'),       opts.title       = 'CMC (Bigot threshold)'; end
if ~isfield(opts, 'colormap'),    opts.colormap    = 'hot'; end
if ~isfield(opts, 'plot'),        opts.plot        = true; end

if ~isfield(params, 'soi_event_indices')
    params.soi_event_indices = [1 size(timewarp_indices, 2)];
end
if ~isfield(params, 'resample_tol'), params.resample_tol = 1e-8; end
if ~isfield(params, 'verbose'),      params.verbose      = true; end

soi_evt = params.soi_event_indices;
alpha   = opts.alpha;

% --- handle 3-D vs 4-D TF_Y ---
if ndims(TF_Y) == 3
    nMuscles = 1;
    TF_Y = reshape(TF_Y, size(TF_Y,1), size(TF_Y,2), size(TF_Y,3), 1);
else
    nMuscles = size(TF_Y, 4);
end

% --- handle 2-D vs 3-D Y_epochs ---
Y_is_3d = (ndims(Y_epochs) == 3);
if ~Y_is_3d
    % Promote to [1 x nSamples x nEpochs]
    Y_epochs = reshape(Y_epochs, 1, size(Y_epochs,1), size(Y_epochs,2));
end
assert(size(Y_epochs, 1) == nMuscles, ...
    'Y_epochs muscle count (%d) does not match TF_Y muscle count (%d).', ...
    size(Y_epochs, 1), nMuscles);

% --- muscle labels default ---
if ~isfield(opts, 'muscle_labels') || isempty(opts.muscle_labels)
    if nMuscles == 4
        opts.muscle_labels = {'VL','RF','GM','BF'};
    else
        opts.muscle_labels = arrayfun(@(m) sprintf('Muscle %d', m), ...
            1:nMuscles, 'UniformOutput', false);
    end
end
assert(numel(opts.muscle_labels) == nMuscles, ...
    'opts.muscle_labels must have %d entries.', nMuscles);

% --- shape checks ---
[nF, nT, nE] = size(TF_X);
assert(isequal([nF, nT, nE], size(TF_Y, 1:3)), ...
    'TF_X and TF_Y must agree on [nFreqs, nTimes, nEpochs].');
assert(length(percent_grid) == nT, 'percent_grid length must match nTimes.');
assert(length(freqs) == nF,        'freqs length must match nFreqs.');
assert(size(X_epochs, 2) == nE,    'X_epochs columns must match nEpochs.');
assert(size(Y_epochs, 3) == nE,    'Y_epochs 3rd dim must match nEpochs.');
assert(size(timewarp_indices, 1) == nE, ...
    'timewarp_indices rows must match nEpochs.');

if params.verbose
    fprintf('compute_and_plot_CMC: %d epochs, %d muscles, alpha = %g\n', ...
        nE, nMuscles, alpha);
end

% =========================================================================
% STEP 1 — Compute MSC and cross-spectrum magnitude per muscle
% =========================================================================
coh           = nan(nF, nT, nMuscles);
cross_mag     = nan(nF, nT, nMuscles);
validX_per_ep = ~any(isnan(reshape(TF_X, [], nE)), 1);   % [1 x nE]

for m = 1:nMuscles
    Y_m = TF_Y(:, :, :, m);
    validY_per_ep = ~any(isnan(reshape(Y_m, [], nE)), 1);
    valid_ep_m    = validX_per_ep & validY_per_ep;

    if ~any(valid_ep_m)
        warning('No valid epochs for muscle %d. Output will be NaN.', m);
        continue;
    end

    Xv = TF_X(:, :, valid_ep_m);
    Yv = Y_m(:, :, valid_ep_m);
    n_valid_m = size(Xv, 3);

    cross_sum = sum(Xv .* conj(Yv), 3);     % [nF x nT] complex
    XX_sum    = sum(abs(Xv).^2, 3);
    YY_sum    = sum(abs(Yv).^2, 3);

    denom = XX_sum .* YY_sum;
    denom(denom == 0) = NaN;

    coh(:, :, m)       = abs(cross_sum).^2 ./ denom;
    cross_mag(:, :, m) = abs(cross_sum) / n_valid_m;       % |S_xy| (mean)
end

% =========================================================================
% STEP 2 — Re-upsample epochs and extract SOI segments of equal length
% =========================================================================
SOI_lengths        = timewarp_indices(:, soi_evt(2)) - timewarp_indices(:, soi_evt(1) + 1);
[max_SOI, max_idx] = max(SOI_lengths);
fs_eff_ratio       = max_SOI ./ SOI_lengths;

if params.verbose
    fprintf('  Extracting SOI segments (max_SOI = %d samples)...\n', max_SOI);
end

% Pre-cache (p, q) per epoch
P_all = zeros(nE, 1);
Q_all = zeros(nE, 1);
for ep = 1:nE
    if SOI_lengths(ep) <= 0 || isnan(SOI_lengths(ep)), continue; end
    [P_all(ep), Q_all(ep)] = rat(fs_eff_ratio(ep), params.resample_tol);
end

X_soi = nan(max_SOI, nE);
Y_soi = nan(max_SOI, nE, nMuscles);

for ep = 1:nE
    if P_all(ep) == 0,   continue; end
    if ~validX_per_ep(ep), continue; end

    p     = P_all(ep);
    q     = Q_all(ep);
    e_on  = timewarp_indices(ep, soi_evt(1));
    e_off = timewarp_indices(ep, soi_evt(2));
    if isnan(e_on) || isnan(e_off) || e_off <= e_on, continue; end

    if ep == max_idx
        x_up    = double(X_epochs(:, ep));
        e_on_up = e_on;
    else
        x_up    = resample(double(X_epochs(:, ep)), p, q);
        e_on_up = round(e_on * p / q);
    end
    idx_end = e_on_up + max_SOI - 1;
    if e_on_up < 1 || idx_end > length(x_up), continue; end
    X_soi(:, ep) = x_up(e_on_up : idx_end);

    Y_ep = squeeze(Y_epochs(:, :, ep))';   % [nSamples x nMuscles]
    if ep == max_idx
        y_up = double(Y_ep);
    else
        y_up = resample(double(Y_ep), p, q);
    end
    if size(y_up, 1) >= idx_end
        Y_soi(:, ep, :) = reshape(y_up(e_on_up : idx_end, :), max_SOI, 1, nMuscles);
    end
end

% =========================================================================
% STEP 3 — Empirical covariance and largest eigenvalue per channel
%
% rho^2 = max eig of (1/n) X X' = max eig of (1/n) X' X
% (the second is n x n, much smaller for nSamples >> nTrials)
% =========================================================================
soi_valid_X = ~any(isnan(X_soi), 1);

% --- per-muscle threshold and mask ---
lambda   = nan(1, nMuscles);
sig_mask = false(nF, nT, nMuscles);

for m = 1:nMuscles
    Y_m_soi      = squeeze(Y_soi(:, :, m));
    soi_valid_Ym = ~any(isnan(Y_m_soi), 1);

    Y_m_TF_valid = ~any(isnan(reshape(TF_Y(:,:,:,m), [], nE)), 1);
    valid_m      = soi_valid_X & soi_valid_Ym & Y_m_TF_valid;

    n_m = sum(valid_m);
    if n_m < 2
        warning('Muscle %d: only %d valid epoch(s). Skipping threshold.', m, n_m);
        continue;
    end

    % Centre across trials, then largest eigenvalue via Gram trick
    X_paired = X_soi(:, valid_m);
    X_paired = X_paired - mean(X_paired, 2);
    G_x      = (X_paired' * X_paired) / n_m;
    G_x      = (G_x + G_x') / 2;
    rho_x_m  = sqrt(max(eig(G_x)));

    Y_paired = Y_m_soi(:, valid_m);
    Y_paired = Y_paired - mean(Y_paired, 2);
    G_y      = (Y_paired' * Y_paired) / n_m;
    G_y      = (G_y + G_y') / 2;
    rho_y_m  = sqrt(max(eig(G_y)));

    % --- Bigot threshold ---
    T_len = max_SOI;
    C     = (-log(alpha/2)) / n_m + sqrt(-2*log(alpha/2) / n_m);
    fac   = (1 + sqrt(T_len / n_m))^2;
    lambda_m = rho_x_m * rho_y_m * fac * C;
    lambda(m) = lambda_m;

    if params.verbose
        fprintf('  Muscle %d (%s): n=%d, rho_x=%.3g, rho_y=%.3g, lambda=%.3g\n', ...
            m, opts.muscle_labels{m}, n_m, rho_x_m, rho_y_m, lambda_m);
    end

    % --- Apply mask ---
    mask_m            = cross_mag(:, :, m) > lambda_m;
    sig_mask(:, :, m) = mask_m;
    coh_m             = coh(:, :, m);
    coh_m(~mask_m)    = NaN;
    coh(:, :, m)      = coh_m;
end

% =========================================================================
% STEP 4 — Plot
% =========================================================================
fig_handle = [];
if ~opts.plot
    return;
end

% Shared colour limits — 90th percentile of surviving (significant) MSC
if ~isfield(opts, 'clim') || isempty(opts.clim)
    cmax = prctile(coh(:), 90);
    if isempty(cmax) || isnan(cmax) || cmax == 0, cmax = 1; end
    opts.clim = [0, cmax];
end

% Subplot layout
nRows = ceil(sqrt(nMuscles));
nCols = ceil(nMuscles / nRows);

fig_handle = figure('Name', opts.title, 'Color', 'w', ...
    'Position', [100 100 1100 750]);

ax_handles = gobjects(nMuscles, 1);
for m = 1:nMuscles
    ax_handles(m) = subplot(nRows, nCols, m);

    coh_m       = coh(:, :, m);
    nan_mask    = isnan(coh_m);
    coh_display = coh_m;
    coh_display(nan_mask) = opts.clim(1);

    h_im = imagesc(percent_grid, freqs, coh_display);
    set(h_im, 'AlphaData', ~nan_mask);     % NaN pixels transparent
    set(gca, 'Color', 'w');                % white background

    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     opts.freq_lim, ...
        'ytick',    opts.freq_ticks, ...
        'box',      'on', ...
        'FontSize', 12);
    colormap(gca, opts.colormap);
    caxis(opts.clim);

    hold on;
    xline(0,   'k--', 'LineWidth', 1.5);
    xline(100, 'k--', 'LineWidth', 1.5);
    hold off;

    if isnan(lambda(m))
        title(sprintf('%s (no threshold)', opts.muscle_labels{m}), ...
            'FontWeight', 'bold');
    else
        title(sprintf('%s  (\\lambda = %.2g)', opts.muscle_labels{m}, lambda(m)), ...
            'FontWeight', 'bold');
    end
    xlabel('Cycle (%)');
    ylabel('Frequency (Hz)');
end

% Single shared colorbar
cb = colorbar(ax_handles(end), 'Position', [0.93 0.11 0.015 0.815]);
ylabel(cb, 'MSC');

% Make room for the colorbar
for m = 1:nMuscles
    pos = get(ax_handles(m), 'Position');
    pos(3) = pos(3) * 0.95;
    set(ax_handles(m), 'Position', pos);
end

sgtitle(sprintf('%s   (\\alpha = %g)', opts.title, alpha), ...
    'FontSize', 14, 'FontWeight', 'bold');

end