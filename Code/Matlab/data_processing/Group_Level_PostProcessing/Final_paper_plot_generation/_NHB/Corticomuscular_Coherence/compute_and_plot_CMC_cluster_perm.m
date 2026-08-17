function [coh, sig_mask, lambda, cluster_info, fig_handle] = ...
    compute_and_plot_CMC_cluster_perm(TF_X, TF_Y, percent_grid, freqs, opts)
% COMPUTE_AND_PLOT_CMC
%   Magnitude-squared coherence (MSC) between an EEG TF representation
%   and one or more EMG TF representations, with significance masking
%   based on a cluster-based permutation test (Maris & Oostenveld 2007).
%
%   Pipeline (per muscle):
%     1. Compute per-pixel real MSC from TF_X and TF_Y (NaN-aware).
%     2. Run nShuffles trial permutations of EMG vs EEG. For each shuffle,
%        compute the shuffled MSC map and STORE it.
%     3. Build a per-pixel cluster-forming threshold t_pix_map(f,t)
%        as the (1-alpha_pix) percentile of each pixel's null
%        distribution. (Default alpha_pix = 0.05.) If opts.t_pix is a
%        scalar, that value is used everywhere instead.
%     4. In the real MSC map, find pixels above t_pix_map. Group
%        spatially-connected ones into clusters (4-connectivity in
%        freq x time). For each cluster compute its MASS = sum of MSC.
%     5. For each shuffled map, apply the same t_pix_map, find clusters,
%        record the MAX cluster mass across the map -> null distribution
%        of "max cluster mass under no dependence", nShuffles long.
%     6. The cluster-level threshold lambda is the (1-alpha_clust)
%        percentile of these null maxima. (Default alpha_clust = 0.05.)
%     7. Real clusters with mass > lambda are significant. The output
%        coh map keeps only pixels belonging to significant clusters;
%        all other pixels are set to NaN.
%
% Inputs:
%   TF_X         = [nFreqs x nTimes x nEpochs] complex EEG TF
%   TF_Y         = [nFreqs x nTimes x nEpochs x nMuscles] complex EMG TF
%                  (or [nFreqs x nTimes x nEpochs] for one muscle)
%   percent_grid = [1 x nTimes] cycle-percent time axis
%   freqs        = [1 x nFreqs] frequency vector (Hz)
%   opts (struct, optional):
%       .alpha_pix     - per-pixel cluster-forming alpha (default 0.05)
%       .alpha_clust   - cluster-level FWE alpha          (default 0.05)
%       .t_pix         - if scalar, fixed cluster-forming threshold;
%                        if [] (default), use per-pixel data-driven t_pix
%       .nShuffles     - number of permutations (default 500)
%       .connectivity  - 4 or 8 (default 4)
%       .rng_seed      - RNG seed (default 'shuffle')
%       .muscle_labels - default {'VL','RF','GM','BF'} or auto
%       .freq_ticks    - default [4 8 14 30 60]
%       .freq_lim      - default [freqs(1) freqs(end)]
%       .title         - default 'CMC (cluster-based permutation)'
%       .colormap      - default 'hot'
%       .plot          - default true
%       .clim          - default [0  prctile(coh(:),90)]
%       .verbose       - default true
%
% Outputs:
%   coh          = [nFreqs x nTimes x nMuscles] MSC, with non-significant
%                  pixels set to NaN
%   sig_mask     = [nFreqs x nTimes x nMuscles] logical mask (true where
%                  pixel belongs to a significant cluster)
%   lambda       = [1 x nMuscles] cluster-mass threshold per muscle
%   cluster_info = [nMuscles x 1] cell array; cluster_info{m} is a struct
%                  array of significant clusters with fields:
%                    .pixels    - linear indices in [nF x nT]
%                    .mass      - sum of MSC values
%                    .p_value   - empirical p-value from the null
%                    .freq_range- [fmin fmax]
%                    .time_range- [tmin tmax]
%   fig_handle   = handle to the figure (empty if opts.plot = false)

% --- defaults ---
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'alpha_pix'),     opts.alpha_pix     = 0.05; end
if ~isfield(opts, 'alpha_clust'),   opts.alpha_clust   = 0.05; end
if ~isfield(opts, 't_pix'),         opts.t_pix         = []; end
if ~isfield(opts, 'nShuffles'),     opts.nShuffles     = 500; end
if ~isfield(opts, 'connectivity'),  opts.connectivity  = 8; end
if ~isfield(opts, 'rng_seed'),      opts.rng_seed      = 'shuffle'; end
if ~isfield(opts, 'freq_ticks'),    opts.freq_ticks    = [4 8 14 30 60]; end
if ~isfield(opts, 'freq_lim'),      opts.freq_lim      = [freqs(1) freqs(end)]; end
if ~isfield(opts, 'title'),         opts.title         = 'CMC (cluster-based permutation)'; end
if ~isfield(opts, 'colormap'),      opts.colormap      = 'turbo'; end
if ~isfield(opts, 'plot'),          opts.plot          = true; end
if ~isfield(opts, 'verbose'),       opts.verbose       = true; end

assert(opts.connectivity == 4 || opts.connectivity == 8, ...
    'opts.connectivity must be 4 or 8.');

nShuffles   = opts.nShuffles;
alpha_pix   = opts.alpha_pix;
alpha_clust = opts.alpha_clust;

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

rng(opts.rng_seed);

if opts.verbose
    fprintf('compute_and_plot_CMC (cluster-based permutation):\n');
    fprintf('  %d epochs, %d muscles, %d shuffles\n', nE, nMuscles, nShuffles);
    fprintf('  alpha_pix = %g, alpha_clust = %g, connectivity = %d\n', ...
        alpha_pix, alpha_clust, opts.connectivity);
    if isempty(opts.t_pix)
        fprintf('  cluster-forming threshold: per-pixel data-driven\n');
    else
        fprintf('  cluster-forming threshold: fixed at %.4g\n', opts.t_pix);
    end
end

% =========================================================================
% Initialise outputs
% =========================================================================
coh          = nan(nF, nT, nMuscles);
sig_mask     = false(nF, nT, nMuscles);
lambda       = nan(1, nMuscles);
cluster_info = cell(nMuscles, 1);

validX_per_ep = ~any(isnan(reshape(TF_X, [], nE)), 1);   % [1 x nE]

% =========================================================================
% Per-muscle processing
% =========================================================================
for m = 1:nMuscles

    if opts.verbose
        fprintf('\n  Muscle %d (%s):\n', m, opts.muscle_labels{m});
    end

    Y_m = TF_Y(:, :, :, m);
    validY_per_ep = ~any(isnan(reshape(Y_m, [], nE)), 1);
    valid_ep_m    = validX_per_ep & validY_per_ep;

    if ~any(valid_ep_m)
        warning('No valid epochs for muscle %d. Output will be NaN.', m);
        continue;
    end

    Xv = TF_X(:, :, valid_ep_m);
    Yv = Y_m(:, :, valid_ep_m);
    n_m = size(Xv, 3);

    if n_m < 2
        warning('Muscle %d: only %d valid epoch(s). Skipping.', m, n_m);
        continue;
    end

    % --- Real MSC ---
    cross_sum_real = sum(Xv .* conj(Yv), 3);
    XX_sum         = sum(abs(Xv).^2, 3);   % invariant under shuffle
    YY_sum         = sum(abs(Yv).^2, 3);   % invariant under shuffle
    denom          = XX_sum .* YY_sum;
    denom(denom == 0) = NaN;

    msc_real = abs(cross_sum_real).^2 ./ denom;
    coh(:, :, m) = msc_real;

    % --- Shuffle loop: store all null MSC maps ---
    if opts.verbose
        fprintf('    Running %d shuffles...\n', nShuffles);
    end

    null_msc = zeros(nF, nT, nShuffles, 'single');   % single to halve memory
    for k = 1:nShuffles
        perm = randperm(n_m);
        Yk   = Yv(:, :, perm);
        cs_k = sum(Xv .* conj(Yk), 3);
        msc_k = abs(cs_k).^2 ./ denom;
        null_msc(:, :, k) = single(msc_k);

        if opts.verbose && rem(k, 100) == 0
            fprintf('      shuffle %d / %d\n', k, nShuffles);
        end
    end

    % --- Per-pixel cluster-forming threshold ---
    if isempty(opts.t_pix)
        t_pix_map = double(prctile(null_msc, 100*(1 - alpha_pix), 3));   % [nF x nT]
    else
        t_pix_map = opts.t_pix * ones(nF, nT);
    end

    % --- Real clusters and their masses ---
    real_above   = msc_real > t_pix_map;
    [real_clusters, real_masses] = find_clusters(real_above, msc_real, opts.connectivity);
    nClust_real  = numel(real_masses);

    if opts.verbose
        fprintf('    Real map: %d clusters above t_pix\n', nClust_real);
    end

    if nClust_real == 0
        if opts.verbose
            fprintf('    No clusters formed in real map -- nothing to test.\n');
        end
        continue;
    end

    % --- Null distribution of max cluster mass ---
    if opts.verbose
        fprintf('    Computing null distribution of max cluster mass...\n');
    end

    null_max_mass = zeros(nShuffles, 1);
    for k = 1:nShuffles
        msc_k     = double(null_msc(:, :, k));
        above_k   = msc_k > t_pix_map;
        [~, m_k]  = find_clusters(above_k, msc_k, opts.connectivity);
        if isempty(m_k)
            null_max_mass(k) = 0;
        else
            null_max_mass(k) = max(m_k);
        end
    end

    % --- Cluster-mass threshold and significance ---
    lambda_m  = prctile(null_max_mass, 100*(1 - alpha_clust));
    lambda(m) = lambda_m;

    if opts.verbose
        fprintf('    Null max-cluster-mass: range [%.3g, %.3g], lambda = %.3g\n', ...
            min(null_max_mass), max(null_max_mass), lambda_m);
    end

    % --- Identify significant clusters ---
    is_sig    = real_masses > lambda_m;
    sig_idx   = find(is_sig);
    nSig      = numel(sig_idx);

    if opts.verbose
        fprintf('    Significant clusters: %d / %d\n', nSig, nClust_real);
    end

    % --- Build sig_mask and cluster_info ---
    mask_m = false(nF, nT);
    info   = struct('pixels', {}, 'mass', {}, 'p_value', {}, ...
                    'freq_range', {}, 'time_range', {});

    for c_i = 1:nClust_real
        if ~is_sig(c_i), continue; end
        pix = real_clusters{c_i};

        if isempty(pix), continue; end
        mask_m(pix) = true;

        [fi, ti] = ind2sub([nF, nT], pix);
        p_val    = mean(null_max_mass >= real_masses(c_i));

        info(end+1).pixels    = pix;            %#ok<AGROW>
        info(end).mass        = real_masses(c_i);
        info(end).p_value     = p_val;
        info(end).freq_range  = [freqs(min(fi)), freqs(max(fi))];
        info(end).time_range  = [percent_grid(min(ti)), percent_grid(max(ti))];
    end

    sig_mask(:, :, m)  = mask_m;
    cluster_info{m}    = info;

    % --- Mask MSC: keep significant clusters only ---
    coh_m         = coh(:, :, m);
    coh_m(~mask_m) = NaN;
    coh(:, :, m)  = coh_m;

    clear null_msc;   % free memory before next muscle
end

% =========================================================================
% Plot
% =========================================================================
fig_handle = [];
if ~opts.plot
    return;
end

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

    h_im = imagesc(percent_grid, freqs, coh_display);
    set(h_im, 'AlphaData', ~nan_mask);
    set(gca, 'Color', 'w');

    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     opts.freq_lim + [0 0.01], ...
        'xlim',     [-5 105], ...
        'ytick',    opts.freq_ticks, ...
        'box',      'on', ...
        'FontSize', 12);
    colormap(gca, opts.colormap);
    caxis(opts.clim);

    hold on;
    xline(0,   'k--', 'LineWidth', 1.5);
    xline(100, 'k--', 'LineWidth', 1.5);
    hold off;

    nSig = 0;
    if ~isempty(cluster_info{m}), nSig = numel(cluster_info{m}); end
    if isnan(lambda(m))
        title(sprintf('%s (no test)', opts.muscle_labels{m}), ...
            'FontWeight', 'bold');
    else
        title(sprintf('%s  (\\lambda = %.2g, %d sig clusters)', ...
            opts.muscle_labels{m}, lambda(m), nSig), 'FontWeight', 'bold');
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

sgtitle(sprintf('%s   (\\alpha_{pix} = %g, \\alpha_{clust} = %g, %d shuffles)', ...
    opts.title, alpha_pix, alpha_clust, nShuffles), ...
    'FontSize', 18, 'FontWeight', 'normal');

end


% =========================================================================
% LOCAL FUNCTION: find_clusters
% =========================================================================
function [pixels_per_cluster, masses] = find_clusters(binary_map, value_map, conn)
%   binary_map : [nF x nT] logical, true where pixel is above threshold
%   value_map  : [nF x nT] real, the values to sum within each cluster
%   conn       : 4 or 8 (pixel connectivity)
%
%   Returns:
%     pixels_per_cluster : cell array, linear indices of each cluster
%     masses             : vector of cluster masses (sum of value_map)

CC = bwconncomp(binary_map, conn);
pixels_per_cluster = CC.PixelIdxList;
nC = numel(pixels_per_cluster);

masses = zeros(nC, 1);
for c = 1:nC
    masses(c) = sum(value_map(pixels_per_cluster{c}));
end
end