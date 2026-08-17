function [coh, sig_mask, lambda, fig_handle] = ...
    compute_and_plot_CMC_with_shuffle(TF_X, TF_Y, percent_grid, freqs, opts)
% COMPUTE_AND_PLOT_CMC
%   Magnitude-squared coherence (MSC) between an EEG TF representation
%   and one or more EMG TF representations, with significance masking
%   based on a permutation/shuffle test on the wavelet cross-spectrum,
%   and a 2x2 plot.
%
%   Pipeline:
%     1. Compute per-pixel MSC and trial-averaged cross-spectrum |S_xy|
%        from TF_X and TF_Y (NaN-aware).
%     2. For each muscle: run nShuffles trial-permutations of EMG vs EEG.
%        For each permutation, compute the shuffled |S_xy_shuf(f,t)| map
%        and record its maximum across (f,t). This builds an empirical
%        null distribution of "max |S_xy| under no dependence".
%     3. The threshold lambda is the (1-alpha) percentile of these
%        nShuffles maxima — a single scalar per muscle that gives
%        family-wise alpha control across the whole TF plane.
%     4. Mask MSC pixels where |S_xy_real| <= lambda.
%
% Inputs:
%   TF_X         = [nFreqs x nTimes x nEpochs] complex EEG TF
%   TF_Y         = [nFreqs x nTimes x nEpochs x nMuscles] complex EMG TF
%                  (or [nFreqs x nTimes x nEpochs] for one muscle)
%   percent_grid = [1 x nTimes] cycle-percent time axis
%   freqs        = [1 x nFreqs] frequency vector (Hz)
%   opts (struct, optional):
%       .alpha         - significance level (default 0.05)
%       .nShuffles     - number of shuffle iterations (default 500)
%       .rng_seed      - RNG seed for reproducibility (default 'shuffle')
%       .muscle_labels - default {'VL','RF','GM','BF'} or auto
%       .freq_ticks    - default [4 8 14 30 60]
%       .freq_lim      - default [freqs(1) freqs(end)]
%       .title         - default 'CMC (shuffle threshold)'
%       .colormap      - default 'hot'
%       .plot          - default true
%       .clim          - default [0  prctile(coh(:),90)]
%       .verbose       - default true
%
% Outputs:
%   coh        = [nFreqs x nTimes x nMuscles] MSC, with non-significant
%                pixels set to NaN
%   sig_mask   = [nFreqs x nTimes x nMuscles] logical mask (true where
%                |S_xy_real| > lambda)
%   lambda     = [1 x nMuscles] scalar threshold per muscle (in units
%                of |S_xy|, on the same scale as TF_X .* conj(TF_Y))
%   fig_handle = handle to the figure (empty if opts.plot = false)

% --- defaults ---
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'alpha'),         opts.alpha         = 0.05; end
if ~isfield(opts, 'nShuffles'),     opts.nShuffles     = 500; end
if ~isfield(opts, 'rng_seed'),      opts.rng_seed      = 'shuffle'; end
if ~isfield(opts, 'freq_ticks'),    opts.freq_ticks    = [4 8 14 30 60]; end
if ~isfield(opts, 'freq_lim'),      opts.freq_lim      = [freqs(1) freqs(end)]; end
if ~isfield(opts, 'title'),         opts.title         = 'CMC (shuffle threshold)'; end
if ~isfield(opts, 'colormap'),      opts.colormap      = 'turbo'; end
if ~isfield(opts, 'plot'),          opts.plot          = true; end
if ~isfield(opts, 'verbose'),       opts.verbose       = true; end

alpha     = opts.alpha;
nShuffles = opts.nShuffles;

% --- handle 3-D vs 4-D TF_Y ---
if ndims(TF_Y) == 3
    nMuscles = 1;
    TF_Y = reshape(TF_Y, size(TF_Y,1), size(TF_Y,2), size(TF_Y,3), 1);
else
    nMuscles = size(TF_Y, 4);
end

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

% --- RNG ---
rng(opts.rng_seed);

if opts.verbose
    fprintf('compute_and_plot_CMC: %d epochs, %d muscles\n', nE, nMuscles);
    fprintf('  alpha = %g, nShuffles = %d\n', alpha, nShuffles);
end

% =========================================================================
% STEP 1 — Compute MSC and cross-spectrum magnitude per muscle
% =========================================================================
coh           = nan(nF, nT, nMuscles);
cross_mag     = nan(nF, nT, nMuscles);   % |S_xy| (averaged across trials)
validX_per_ep = ~any(isnan(reshape(TF_X, [], nE)), 1);

% Per-muscle storage of the trial-restricted TF arrays for the shuffle step
TFX_valid_per_m = cell(nMuscles, 1);
TFY_valid_per_m = cell(nMuscles, 1);
nValid_per_m    = zeros(nMuscles, 1);

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

    cross_sum = sum(Xv .* conj(Yv), 3);
    XX_sum    = sum(abs(Xv).^2, 3);
    YY_sum    = sum(abs(Yv).^2, 3);

    denom = XX_sum .* YY_sum;
    denom(denom == 0) = NaN;

    coh(:, :, m)       = abs(cross_sum).^2 ./ denom;
    % cross_mag(:, :, m) = abs(cross_sum) / n_valid_m;

    TFX_valid_per_m{m} = Xv;
    TFY_valid_per_m{m} = Yv;
    nValid_per_m(m)    = n_valid_m;
end

% =========================================================================
% STEP 2 — Shuffle test per muscle
% =========================================================================
lambda   = nan(1, nMuscles);
sig_mask = false(nF, nT, nMuscles);

for m = 1:nMuscles
    n_m = nValid_per_m(m);
    if n_m < 2
        warning('Muscle %d: only %d valid epoch(s). Skipping threshold.', m, n_m);
        continue;
    end

    if opts.verbose
        fprintf('  Muscle %d (%s): %d shuffles on %d trials...\n', ...
            m, opts.muscle_labels{m}, nShuffles, n_m);
    end

    Xv = TFX_valid_per_m{m};   % [nF x nT x n_m]
    Yv = TFY_valid_per_m{m};   % [nF x nT x n_m]

    max_null = zeros(nShuffles, 1);

    for k = 1:nShuffles
        perm = randperm(n_m);
        Yk   = Yv(:, :, perm);

        % cross_sum_k = sum(Xv .* conj(Yk), 3);    % [nF x nT] complex
        % cross_mag_k = abs(cross_sum_k) / n_m;    % [nF x nT] real
        % max_null(k) = max(cross_mag_k(:));

        cross_sum_k = sum(Xv .* conj(Yk), 3);
        XX_sum_k    = sum(abs(Xv).^2, 3);          % unchanged across shuffles, but cheap
        YY_sum_k    = sum(abs(Yk).^2, 3);          % unchanged across shuffles
        msc_k       = abs(cross_sum_k).^2 ./ (XX_sum_k .* YY_sum_k);
        max_null(k) = max(msc_k(:));

        if opts.verbose && rem(k, 100) == 0
            fprintf('    shuffle %d / %d\n', k, nShuffles);
        end
    end

    % --- Threshold: (1 - alpha) percentile of the max distribution ---
    lambda_m  = prctile(max_null, 100 * (1 - alpha));
    lambda(m) = lambda_m;

    if opts.verbose
        fprintf('    lambda = %.4g  (range of null max: %.4g to %.4g)\n', ...
            lambda_m, min(max_null), max(max_null));
    end

    % --- Apply mask ---
    mask_m            = coh(:, :, m) > lambda_m;
    sig_mask(:, :, m) = mask_m;
    coh_m             = coh(:, :, m);
    coh_m(~mask_m)    = NaN;
    coh(:, :, m)      = coh_m;
end

% =========================================================================
% STEP 3 — Plot
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

nRows = ceil(sqrt(nMuscles));
nCols = ceil(nMuscles / nRows);

% fig_handle = figure('Name', opts.title, 'Color', 'w', ...
%     'Position', [100 100 1100 750]);
fig_handle = figure('Name', opts.title, 'Color', 'w');

ax_handles = gobjects(nMuscles, 1);
for m = 1:nMuscles
    ax_handles(m) = subplot(nRows, nCols, m);

    coh_m       = coh(:, :, m);
    nan_mask    = isnan(coh_m);
    coh_display = coh_m;
    coh_display(nan_mask) = opts.clim(1);

    % h_im = contourf(percent_grid, freqs, coh_display, 200, 'LineColor', 'none');
    h_im = imagesc(percent_grid, freqs, coh_display);

    set(h_im, 'AlphaData', ~nan_mask);
    set(gca, 'Color', 'w');

    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     opts.freq_lim  + [0 0.01], ...
        'ytick',    opts.freq_ticks, ...
        'xlim',     [-5 105], ... %percent
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

% cb = colorbar(ax_handles(end), 'Position', [0.93 0.11 0.015 0.815]);
cb = colorbar(ax_handles(end), 'Position', [0.93 0.13 0.015 0.755]);

ylabel(cb, 'MSC');

% for m = 1:nMuscles
%     pos = get(ax_handles(m), 'Position');
%     pos(3) = pos(3) * 0.95;
%     set(ax_handles(m), 'Position', pos);
% end

sgtitle(sprintf('%s   (\\alpha = %g, %d shuffles)', ...
    opts.title, alpha, nShuffles), ...
    'FontSize', 18, 'FontWeight', 'bold');

end