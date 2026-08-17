function [coh, fig_handle] = compute_and_plot_CMC(TF_X, TF_Y, percent_grid, freqs, opts)
% COMPUTE_AND_PLOT_CMC
%   Magnitude-squared coherence (MSC) between an EEG TF representation
%   and one or more EMG TF representations, plus a 2x2 plot.
%
% Inputs:
%   TF_X         = [nFreqs x nTimes x nEpochs] complex TF of EEG (one channel)
%   TF_Y         = [nFreqs x nTimes x nEpochs x nMuscles] complex TF of EMG
%                  (or [nFreqs x nTimes x nEpochs] for a single muscle)
%   percent_grid = [1 x nTimes] cycle-percent time axis
%   freqs        = [1 x nFreqs] frequency vector (Hz)
%   opts (struct, optional):
%       .muscle_labels = cellstr, default {'VL','RF','GM','BF'} (or auto)
%       .freq_ticks    = default [4 8 14 30 60]
%       .freq_lim      = default [freqs(1) freqs(end)]
%       .title         = figure super-title (default 'CMC')
%       .colormap      = colormap name (default 'hot')
%       .plot          = true/false to draw figure (default true)
%       .clim          = [cmin cmax] colour limits, default [0 max(coh(:))]
%
% Outputs:
%   coh        = [nFreqs x nTimes x nMuscles] real-valued MSC, in [0, 1]
%   fig_handle = handle to the figure (empty if opts.plot = false)
%
% Formula (MSC):
%   R^2(f,t) = |sum_i  W_X(f,t,i) * conj(W_Y(f,t,i))|^2
%              / ( sum_i |W_X(f,t,i)|^2  *  sum_i |W_Y(f,t,i)|^2 )
%   Trials with NaN values are excluded (per (f,t) sample).

% --- defaults ---
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'freq_ticks'), opts.freq_ticks = [4 8 14 30 60]; end
if ~isfield(opts, 'freq_lim'),   opts.freq_lim   = [freqs(1) freqs(end)]; end
if ~isfield(opts, 'title'),      opts.title      = 'Cortico-Muscular Coherence (MSC)'; end
if ~isfield(opts, 'colormap'),   opts.colormap   = 'turbo'; end
if ~isfield(opts, 'plot'),       opts.plot       = true; end

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
    'opts.muscle_labels must have %d entries (one per muscle).', nMuscles);

% --- shape checks ---
[nF, nT, nE] = size(TF_X);
assert(isequal([nF, nT, nE], size(TF_Y, 1:3)), ...
    'TF_X and TF_Y must agree on [nFreqs, nTimes, nEpochs].');
assert(length(percent_grid) == nT, 'percent_grid length must match nTimes.');
assert(length(freqs)        == nF, 'freqs length must match nFreqs.');

% =========================================================================
% Compute MSC per muscle
% =========================================================================
coh = zeros(nF, nT, nMuscles);

% Mask of valid (non-NaN) X samples — same across muscles
validX = ~isnan(TF_X);

for m = 1:nMuscles

    Y_m    = TF_Y(:, :, :, m);
    validY = ~isnan(Y_m);
    valid  = validX & validY;        % [nF x nT x nE]

    % Replace NaNs with 0 so they contribute nothing to the sums
    X_use = TF_X;  X_use(~valid) = 0;
    Y_use = Y_m;   Y_use(~valid) = 0;

    % Sums over trials (3rd dim)
    cross_sum = sum(X_use .* conj(Y_use), 3);     % [nF x nT] complex
    XX_sum    = sum(abs(X_use).^2, 3);            % [nF x nT] real
    YY_sum    = sum(abs(Y_use).^2, 3);            % [nF x nT] real

    denom = XX_sum .* YY_sum;
    denom(denom == 0) = NaN;                      % avoid 0/0

    coh(:, :, m) = abs(cross_sum).^2 ./ denom;
    % coh(:, :, m) = abs(cross_sum) ./ sqrt(denom);
end

% =========================================================================
% Plot 2x2 grid
% =========================================================================
fig_handle = [];
if ~opts.plot
    return;
end

% Shared colour limits across all muscles
if ~isfield(opts, 'clim') || isempty(opts.clim)
    cmax = prctile(coh(:), 90);
    if isempty(cmax) || isnan(cmax) || cmax == 0, cmax = 1; end
    opts.clim = [0, cmax];
end

% Determine subplot grid
nRows = ceil(sqrt(nMuscles));
nCols = ceil(nMuscles / nRows);

% fig_handle = figure('Name', opts.title, 'Color', 'w', ...
%     'Position', [100 100 1100 750]);
fig_handle = figure('Name', opts.title, 'Color', 'w');

ax_handles = gobjects(nMuscles, 1);
for m = 1:nMuscles
    ax_handles(m) = subplot(nRows, nCols, m);
    contourf(percent_grid, freqs, coh(:, :, m), 200, 'LineColor', 'none');

    % set(gca, ...
    %     'yscale',   'log', ...
    %     'ydir',     'norm', ...
    %     'ylim',     opts.freq_lim, ...
    %     'ytick',    opts.freq_ticks, ...
    %     'box',      'on', ...
    %     'FontSize', 12);
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
    xline(0,   'w--', 'LineWidth', 1.5);
    xline(100, 'w--', 'LineWidth', 1.5);
    hold off;

    title(opts.muscle_labels{m}, 'FontWeight', 'bold');
    xlabel('Cycle (%)');
    ylabel('Frequency (Hz)');
end

% Single shared colorbar
cb = colorbar(ax_handles(end), 'Position', [0.93 0.13 0.015 0.755]);
% cb = colorbar(ax_handles(end));
ylabel(cb, 'MSC');
% ylabel(cb, 'Coherence');

% Make room for the colorbar by shifting subplots left a touch
% for m = 1:nMuscles
%     pos = get(ax_handles(m), 'Position');
%     pos(3) = pos(3) * 0.95;
%     set(ax_handles(m), 'Position', pos);
% end

sgtitle(opts.title, 'FontSize', 18, 'FontWeight', 'normal');

end