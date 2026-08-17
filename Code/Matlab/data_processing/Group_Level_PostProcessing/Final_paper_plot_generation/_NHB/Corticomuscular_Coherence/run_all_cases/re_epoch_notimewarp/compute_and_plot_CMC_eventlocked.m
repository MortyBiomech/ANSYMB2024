function [coh, fig_handle] = compute_and_plot_CMC_eventlocked( ...
        TF_X, TF_Y, time_axis_ms, freqs, opts)
% COMPUTE_AND_PLOT_CMC_EVENTLOCKED
%   Same MSC computation as compute_and_plot_CMC, but plotting is designed
%   for event-locked (time in ms) rather than cycle-% axes:
%     - x-axis spans the full time_axis_ms range
%     - x-axis label: "Time relative to event (ms)"
%     - a single dashed vertical line at t = 0 (the locking event)
%     - optional event label on the t = 0 marker
%
% Inputs:
%   TF_X         = [nFreqs x nTimes x nEpochs] complex EEG TF
%   TF_Y         = [nFreqs x nTimes x nEpochs x nMuscles] complex EMG TF
%                  (or [nFreqs x nTimes x nEpochs] for a single muscle)
%   time_axis_ms = [1 x nTimes] time axis in milliseconds
%   freqs        = [1 x nFreqs] frequency axis in Hz
%   opts (struct, optional):
%       .muscle_labels = cellstr   (default {'VL','RF','GM','BF'} or auto)
%       .freq_ticks    = default [4 8 14 30 60]
%       .freq_lim      = default [freqs(1) freqs(end)]
%       .title         = figure super-title (default 'CMC — event-locked')
%       .event_label   = string label for the t=0 marker (default 't=0')
%       .colormap      = default 'turbo'
%       .plot          = true/false (default true)
%       .clim          = [cmin cmax] (default auto from 90th percentile)
%
% Outputs:
%   coh        = [nFreqs x nTimes x nMuscles] real-valued MSC in [0, 1]
%   fig_handle = handle to the figure (empty if opts.plot = false)

% --- defaults ---
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'freq_ticks'),   opts.freq_ticks   = [4 8 14 30 60]; end
if ~isfield(opts, 'freq_lim'),     opts.freq_lim     = [freqs(1) freqs(end)]; end
if ~isfield(opts, 'title'),        opts.title        = 'CMC — event-locked'; end
if ~isfield(opts, 'event_label'),  opts.event_label  = 't = 0'; end
if ~isfield(opts, 'colormap'),     opts.colormap     = 'turbo'; end
if ~isfield(opts, 'plot'),         opts.plot         = true; end

% --- handle 3-D vs 4-D TF_Y ---
if ndims(TF_Y) == 3
    nMuscles = 1;
    TF_Y = reshape(TF_Y, size(TF_Y,1), size(TF_Y,2), size(TF_Y,3), 1);
else
    nMuscles = size(TF_Y, 4);
end

% --- muscle labels ---
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
assert(length(time_axis_ms) == nT, 'time_axis_ms length must match nTimes.');
assert(length(freqs)        == nF, 'freqs length must match nFreqs.');

% =========================================================================
%  Compute MSC
% =========================================================================
coh   = zeros(nF, nT, nMuscles);
validX = ~isnan(TF_X);

for m = 1:nMuscles
    Y_m    = TF_Y(:, :, :, m);
    validY = ~isnan(Y_m);
    valid  = validX & validY;

    X_use = TF_X;  X_use(~valid) = 0;
    Y_use = Y_m;   Y_use(~valid) = 0;

    cross_sum = sum(X_use .* conj(Y_use), 3);
    XX_sum    = sum(abs(X_use).^2, 3);
    YY_sum    = sum(abs(Y_use).^2, 3);

    denom = XX_sum .* YY_sum;
    denom(denom == 0) = NaN;

    coh(:, :, m) = abs(cross_sum).^2 ./ denom;
end

% =========================================================================
%  Plot
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

% x-axis limits: exact time axis range — no margin, no white bands
xlim_vals = [time_axis_ms(1), time_axis_ms(end)];

nRows = ceil(sqrt(nMuscles));
nCols = ceil(nMuscles / nRows);

fig_handle = figure('Name', opts.title, 'Color', 'w');

ax_handles = gobjects(nMuscles, 1);
for m = 1:nMuscles
    ax_handles(m) = subplot(nRows, nCols, m);
    contourf(time_axis_ms, freqs, coh(:, :, m), 200, 'LineColor', 'none');

    set(gca, ...
        'yscale',   'log', ...
        'ydir',     'norm', ...
        'ylim',     opts.freq_lim + [0 0.01], ...
        'xlim',     xlim_vals, ...
        'ytick',    opts.freq_ticks, ...
        'box',      'on', ...
        'FontSize', 12);
    colormap(gca, opts.colormap);
    caxis(opts.clim);

    hold on;
    % Single dashed marker at the locking event (t = 0); label in title
    xline(0, 'w--', 'LineWidth', 2);
    hold off;

    title(opts.muscle_labels{m}, 'FontWeight', 'bold');
    xlabel('Time relative to event (ms)');
    ylabel('Frequency (Hz)');
end

cb = colorbar(ax_handles(end), 'Position', [0.93 0.13 0.015 0.755]);
ylabel(cb, 'MSC');

sgtitle(opts.title, 'FontSize', 16, 'FontWeight', 'normal');

end