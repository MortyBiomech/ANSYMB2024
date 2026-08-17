function [TF_X, TF_Y, percent_grid, freqs] = ...
    compute_TF_timewarp_posthoc(X_epochs, Y_epochs, timewarp_indices, params, tlimits_ms)
% COMPUTE_TF_TIMEWARP_POSTHOC
%   Compute wavelet TF at the native sampling rate (no upsampling), then
%   warp each epoch's complex TF matrix to a common cycle-% grid by
%   linear interpolation of real and imaginary parts separately.
%
%   This implements the "calculate each X(t,f,k) & Y(t,f,k), then timewarp
%   the cross-spectrum & auto-spectrum and calculate the CMC" approach from
%   the design diagram — specifically the variant where the TF matrices
%   themselves are normalised to a common %-grid before coherence is computed.
%
%   Compared with compute_TF_timenorm (which warps the signal before TF),
%   this approach preserves the original wavelet scale at each time point
%   but introduces mild interpolation artefacts in the complex plane.
%   Interpolating real and imaginary parts independently (MATLAB's interp1
%   with complex input) is equivalent and numerically stable.
%
% Inputs:
%   X_epochs         : [nSamples x nEpochs]
%   Y_epochs         : [nMuscles x nSamples x nEpochs]  or  [nSamples x nEpochs]
%   timewarp_indices : [nEpochs x nEvents] — sample indices of timewarp events
%                      (column 1 = SOI onset, column end = SOI offset)
%   params           : struct  (.srate, .freqs, .nfreqs, .cycles,
%                               .freqscale, .ntimesout)
%                      Optional fields:
%                        .soi_event_indices — [onset_col, offset_col],
%                                            default [1, size(tw,2)]
%                        .resample_tol     — (unused here, kept for API compat)
%                        .verbose          — default true
%   tlimits_ms       : [tmin tmax] in milliseconds (original epoch window)
%
% Outputs:
%   TF_X        : [nFreqs x nTimesOut x nEpochs]           complex EEG TF
%   TF_Y        : [nFreqs x nTimesOut x nEpochs x nMuscles] complex EMG TF
%   percent_grid: [1 x nTimesOut] cycle-% axis (0 = SOI onset, 100 = SOI offset)
%   freqs       : [1 x nFreqs] frequency axis in Hz

% --- defaults ---
if ~isfield(params, 'verbose'),          params.verbose          = true; end
if ~isfield(params, 'soi_event_indices')
    params.soi_event_indices = [1, size(timewarp_indices, 2)];
end
soi_evt = params.soi_event_indices;

% --- handle 2-D vs 3-D Y_epochs ---
Y_is_3d = (ndims(Y_epochs) == 3);
if Y_is_3d
    [nMuscles, nSamples_Y, nEpochs_Y] = size(Y_epochs);
else
    [nSamples_Y, nEpochs_Y] = size(Y_epochs);
    nMuscles = 1;
    Y_epochs = reshape(Y_epochs, 1, nSamples_Y, nEpochs_Y);
end
[nSamples_X, nEpochs_X] = size(X_epochs);
assert(nEpochs_X == nEpochs_Y,   'X and Y must have the same number of epochs.');
assert(nSamples_X == nSamples_Y, 'X and Y must have the same number of samples.');
nEpochs = nEpochs_X;

% ══════════════════════════════════════════════════════════════════════════
%  STEP 1 — Compute TF for all epochs at native srate (batch mode)
% ══════════════════════════════════════════════════════════════════════════
if params.verbose
    fprintf('compute_TF_timewarp_posthoc: computing TF at native srate...\n');
end

[TF_X_native, freqs, time_grid_ms] = timefreq( ...
    double(X_epochs), params.srate, ...
    'tlimits',   tlimits_ms, ...
    'freqs',     params.freqs, ...
    'nfreqs',    params.nfreqs, ...
    'ntimesout', params.ntimesout, ...
    'freqscale', params.freqscale, ...
    'cycles',    params.cycles);
% TF_X_native: [nFreqs x nTimes_native x nEpochs]

nF        = length(freqs);
nT_native = length(time_grid_ms);
time_grid_ms = time_grid_ms(:)';   % ensure row vector

TF_Y_native = nan(nF, nT_native, nEpochs, nMuscles);
for m = 1:nMuscles
    Y_m = squeeze(Y_epochs(m, :, :));   % [nSamples x nEpochs]
    [tf_m, ~, ~] = timefreq( ...
        double(Y_m), params.srate, ...
        'tlimits',   tlimits_ms, ...
        'freqs',     params.freqs, ...
        'nfreqs',    params.nfreqs, ...
        'ntimesout', params.ntimesout, ...
        'freqscale', params.freqscale, ...
        'cycles',    params.cycles);
    TF_Y_native(:, :, :, m) = tf_m;
end

% ══════════════════════════════════════════════════════════════════════════
%  STEP 2 — Map timewarp event sample indices to TF time-point indices
% ══════════════════════════════════════════════════════════════════════════
% Sample e at epoch time:
%   t_ms(e) = tlimits_ms(1) + (e - 1) / srate * 1000
% Find nearest output TF time point via knnsearch.

srate = params.srate;
t_events_ms = tlimits_ms(1) + (timewarp_indices - 1) / srate * 1000;
% t_events_ms: [nEpochs x nEvents]

tf_event_idx = knnsearch(time_grid_ms', t_events_ms(:));
tf_event_idx = reshape(tf_event_idx, size(timewarp_indices));
% tf_event_idx: [nEpochs x nEvents] — indices into time_grid_ms

% ══════════════════════════════════════════════════════════════════════════
%  STEP 3 — Build the reference cycle-% grid from the longest-SOI epoch
% ══════════════════════════════════════════════════════════════════════════
SOI_tf_lengths = tf_event_idx(:, soi_evt(2)) - tf_event_idx(:, soi_evt(1));
[~, max_idx]   = max(SOI_tf_lengths);

e_on_ref  = tf_event_idx(max_idx, soi_evt(1));
e_off_ref = tf_event_idx(max_idx, soi_evt(2));

% Reference % axis spans the full native TF time grid
soi_pct_ref   = linspace(0, 100, e_off_ref - e_on_ref + 1);
delta_pct_ref = soi_pct_ref(2) - soi_pct_ref(1);
before_ref    = (-delta_pct_ref*(e_on_ref - 1)) : delta_pct_ref : -delta_pct_ref;
after_ref     = (100 + delta_pct_ref) : delta_pct_ref : ...
                (100 + (nT_native - e_off_ref) * delta_pct_ref);
ref_pct_full  = [before_ref, soi_pct_ref, after_ref];
% ref_pct_full has length nT_native — one % value per native TF time point

% Common output % grid (evenly spaced, ntimesout points)
nTout        = params.ntimesout;
percent_grid = linspace(ref_pct_full(1), ref_pct_full(end), nTout);

if params.verbose
    fprintf('  Reference epoch %d: SOI = %d TF time points, delta_pct = %.4f%%\n', ...
        max_idx, e_off_ref - e_on_ref + 1, delta_pct_ref);
    fprintf('  Common %% grid: %.2f%% to %.2f%% (%d points)\n', ...
        percent_grid(1), percent_grid(end), nTout);
end

% ══════════════════════════════════════════════════════════════════════════
%  STEP 4 — Allocate output arrays and warp each epoch
% ══════════════════════════════════════════════════════════════════════════
TF_X = nan(nF, nTout, nEpochs);
if Y_is_3d
    TF_Y = nan(nF, nTout, nEpochs, nMuscles);
else
    TF_Y = nan(nF, nTout, nEpochs);
end

skipped = 0;

for ep = 1:nEpochs

    if params.verbose && rem(ep, 200) == 0
        fprintf('  Warping epoch %d / %d...\n', ep, nEpochs);
    end

    e_on  = tf_event_idx(ep, soi_evt(1));
    e_off = tf_event_idx(ep, soi_evt(2));
    if isnan(e_on) || isnan(e_off) || e_off <= e_on
        skipped = skipped + 1;
        continue;
    end

    % Build per-epoch % axis (length = nT_native, same as native TF grid)
    soi_pct_ep   = linspace(0, 100, e_off - e_on + 1);
    d_ep         = soi_pct_ep(2) - soi_pct_ep(1);
    before_ep    = (-d_ep*(e_on - 1)) : d_ep : -d_ep;
    after_ep     = (100 + d_ep) : d_ep : (100 + (nT_native - e_off) * d_ep);
    pct_ep       = [before_ep, soi_pct_ep, after_ep];
    % pct_ep: [1 x nT_native]

    % Guard against edge cases where pct_ep length drifts by 1 sample
    if length(pct_ep) ~= nT_native
        skipped = skipped + 1;
        continue;
    end

    % Warp X: interpolate complex TF from pct_ep → percent_grid
    % MATLAB interp1 handles complex data by interpolating Re and Im separately.
    TF_ep_X = TF_X_native(:, :, ep);   % [nF x nT_native]
    for f = 1:nF
        TF_X(f, :, ep) = interp1(pct_ep, TF_ep_X(f, :), percent_grid, 'linear', NaN);
    end

    % Warp Y
    for m = 1:nMuscles
        TF_ep_Y = TF_Y_native(:, :, ep, m);   % [nF x nT_native]
        for f = 1:nF
            if Y_is_3d
                TF_Y(f, :, ep, m) = interp1(pct_ep, TF_ep_Y(f, :), percent_grid, 'linear', NaN);
            else
                TF_Y(f, :, ep) = interp1(pct_ep, TF_ep_Y(f, :), percent_grid, 'linear', NaN);
            end
        end
    end

end

if params.verbose
    nValid = sum(~squeeze(isnan(TF_X(1, 1, :))));
    fprintf('  Done. %d / %d epochs warped (%d skipped).\n', ...
        nValid, nEpochs, skipped);
    fprintf('  TF_X: [%s]\n', num2str(size(TF_X)));
end

end
