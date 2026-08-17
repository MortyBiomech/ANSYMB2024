function [TF_X, TF_Y, percent_grid, freqs, fs_eff, mid_event_pct] = ...
    compute_TF_timenorm(X_epochs, Y_epochs, timewarp_indices, params)
% COMPUTE_TF_TIMENORM
%   Time-frequency decomposition of EEG and EMG epochs with Fauvet-style
%   duration normalisation, using a cycle-percent time axis.
%
%   Implements the approach from Morteza Khosrotabar's test 4:
%   each epoch is upsampled so its SOI occupies the same number of samples
%   as the longest-SOI epoch, then timefreq() is called with tlimits
%   expressed in % of cycle (event 1 = 0 %, event 3 = 100 %).
%   All epochs share the same output % grid via the 'timesout' parameter,
%   so the resulting TF matrices have identical [freq x time] dimensions
%   and can be stacked into [freq x time x trials] for coherence.
%
% Inputs:
%   X_epochs         = [nSamples x nEpochs] - first signal (e.g. EEG IC)
%                      Already aligned to a common time grid at params.srate.
%   Y_epochs         = [nMuscles x nSamples x nEpochs] OR [nSamples x nEpochs]
%                      Second signal (e.g. EMG). If 3-D, the function loops
%                      over muscles and returns TF_Y as a 4-D array.
%   timewarp_indices = [nEpochs x nEvents] - per-epoch sample indices of the
%                      timewarp events. Must have at least 2 columns.
%                      Example with 3 events:
%                        column 1: SOI onset event (e.g. FlxS)
%                        column 2: middle event   (e.g. FlexE/ExtS)
%                        column 3: SOI offset event (e.g. ExtE)
%   params (struct)  = .srate       (Hz, original sampling rate)
%                      .tlimits     ([tmin tmax] ms, original epoch limits)
%                      .freqs       ([fmin fmax] Hz)
%                      .nfreqs      (number of frequencies)
%                      .cycles      (Morlet cycles, e.g. [3 0.8])
%                      .freqscale   ('log' or 'linear')
%                      .ntimesout   (number of % grid points)
%                      .soi_event_indices ([2-element vector] which columns
%                                          of timewarp_indices define the
%                                          SOI onset and offset.
%                                          Default: [1 size(timewarp_indices,2)]
%                                          Examples:
%                                            [1 3] = first to last event
%                                            [1 2] = first to middle event
%                                            [2 3] = middle to last event)
%                      .resample_tol (rational approx tolerance, default 1e-8)
%                      .verbose     (true/false, default true)
%
% Outputs:
%   TF_X         = [nFreqs x nTimesOut x nEpochs] complex TF of X
%   TF_Y         = if Y_epochs is 2-D: [nFreqs x nTimesOut x nEpochs]
%                  if Y_epochs is 3-D: [nFreqs x nTimesOut x nEpochs x nMuscles]
%   percent_grid = [1 x nTimesOut] cycle-percent time axis (0 % = SOI onset,
%                  100 % = SOI offset, negative/>100 = guard regions)
%   freqs        = [1 x nFreqs] frequency vector (Hz)
%   fs_eff       = [nEpochs x 1] effective sampling rate of each upsampled
%                  epoch (Hz)
%   mid_event_pct = [nEpochs x 1] % grid value of the middle event for each
%                  epoch, snapped to the nearest percent_grid sample.
%                  ONLY populated when timewarp_indices has exactly 3 events
%                  AND soi_event_indices = [1 3]. Otherwise returned as [].
%                  NaN for skipped epochs.
%
% Notes:
%   - The longest-SOI epoch (where SOI is defined by soi_event_indices)
%     is used as the reference. Its event indices define the % grid
%     used by all epochs.
%   - Epochs with zero or invalid SOI (offset <= onset) are skipped
%     and their TF columns are filled with NaN.
%   - Only the two events specified by soi_event_indices are used for
%     normalisation. Other events in timewarp_indices are ignored
%     (but may be useful for plotting markers afterwards).

% --- defaults ---
if ~isfield(params, 'resample_tol'), params.resample_tol = 1e-8; end
if ~isfield(params, 'verbose'),      params.verbose      = true; end
if ~isfield(params, 'freqscale'),    params.freqscale    = 'log'; end
if ~isfield(params, 'ntimesout'),    params.ntimesout    = 400;   end
if ~isfield(params, 'soi_event_indices')
    params.soi_event_indices = [1 size(timewarp_indices, 2)];
end

srate    = params.srate;
freqs_in = params.freqs;
nfreqs   = params.nfreqs;
cycles   = params.cycles;
nTout    = params.ntimesout;
soi_evt  = params.soi_event_indices;

% Validate soi_event_indices
assert(numel(soi_evt) == 2, ...
    'params.soi_event_indices must be a 2-element vector [onset_col, offset_col].');
assert(soi_evt(1) ~= soi_evt(2), ...
    'params.soi_event_indices: onset and offset columns must differ.');
assert(all(soi_evt >= 1) && all(soi_evt == round(soi_evt)), ...
    'params.soi_event_indices entries must be positive integers.');

% --- decide whether to track the middle event ---
% Only populated when timewarp_indices has exactly 3 events AND the SOI
% spans columns 1 and 3 (i.e. soi_event_indices = [1 3]).
track_mid_event = (size(timewarp_indices, 2) == 3) && ...
                  isequal(sort(soi_evt(:))', [1 3]);

% --- handle 2-D vs 3-D Y_epochs ---
Y_is_3d = (ndims(Y_epochs) == 3);
if Y_is_3d
    [nMuscles, nSamples_Y, nEpochs_Y] = size(Y_epochs);
else
    nMuscles  = 1;
    [nSamples_Y, nEpochs_Y] = size(Y_epochs);
end

[nSamples_X, nEpochs_X] = size(X_epochs);
assert(nEpochs_X == nEpochs_Y, ...
    'X_epochs and Y_epochs must have the same number of epochs.');
assert(nSamples_X == nSamples_Y, ...
    'X_epochs and Y_epochs must have the same number of samples per epoch.');
nEpochs  = nEpochs_X;
nSamples = nSamples_X;
assert(size(timewarp_indices, 1) == nEpochs, ...
    'timewarp_indices must have one row per epoch.');
assert(size(timewarp_indices, 2) >= max(soi_evt), ...
    'timewarp_indices must have at least %d columns to use soi_event_indices = [%d %d].', ...
    max(soi_evt), soi_evt(1), soi_evt(2));

% =========================================================================
% STEP 1 — Identify the longest-SOI epoch and build the reference % grid
% =========================================================================
SOI_lengths = timewarp_indices(:, soi_evt(2)) - timewarp_indices(:, soi_evt(1));   % samples
[max_SOI, max_idx] = max(SOI_lengths);

if params.verbose
    fprintf('compute_TF_timenorm:\n');
    fprintf('  %d epochs, %d samples each at %.0f Hz\n', nEpochs, nSamples, srate);
    fprintf('  SOI defined by events [%d, %d] of timewarp_indices\n', ...
        soi_evt(1), soi_evt(2));
    fprintf('  SOI length: min = %d, max = %d (epoch %d), mean = %.1f samples\n', ...
        min(SOI_lengths), max_SOI, max_idx, mean(SOI_lengths));
end

% Reference epoch event indices (in samples)
e_on_ref  = timewarp_indices(max_idx, soi_evt(1));   % SOI onset
e_off_ref = timewarp_indices(max_idx, soi_evt(2));   % SOI offset

% Build the reference % axis spanning the full reference epoch:
% inside SOI:   linspace(0, 100, e_off_ref - e_on_ref + 1)
% before SOI:   negative values with same delta
% after  SOI:   values > 100 with same delta
soi_pct       = linspace(0, 100, e_off_ref - e_on_ref + 1);
delta_pct     = soi_pct(2) - soi_pct(1);
before_pct    = -delta_pct*(e_on_ref - 1) : delta_pct : -delta_pct;
after_pct     = 100 + delta_pct : delta_pct : 100 + (nSamples - e_off_ref) * delta_pct;
ref_pct_full  = [before_pct, soi_pct, after_pct];

if params.verbose
    fprintf('  Reference %% axis: %.2f%% to %.2f%% (delta = %.4f%%)\n', ...
        ref_pct_full(1), ref_pct_full(end), delta_pct);
end

% =========================================================================
% STEP 2 — Compute upsampling ratios and allocate output arrays.
% =========================================================================
fs_eff_ratio = max_SOI ./ SOI_lengths;   % >= 1 (= 1 for max_idx)
fs_eff       = srate * fs_eff_ratio;

% =========================================================================
% STEP 3 — Compute TF of the reference epoch.
%          The reference epoch (longest SOI) is not upsampled (ratio = 1).
%          Its timefreq() call defines the % output grid (percent_grid)
%          that all other epochs are forced to share via 'timesout'.
%          We compute TF on [X | Y_muscles] in a single call and store
%          the result directly into TF_X / TF_Y at index max_idx so it is
%          not recomputed in the main loop.
% =========================================================================
if Y_is_3d
    Y_ref = squeeze(Y_epochs(:, :, max_idx))';   % [nSamples x nMuscles]
else
    Y_ref = Y_epochs(:, max_idx);                % [nSamples x 1]
end
all_ref = [double(X_epochs(:, max_idx)), double(Y_ref)];   % [nSamples x (1+nMuscles)]

[tf_ref, freqs, percent_grid] = timefreq( ...
    all_ref, srate, ...
    'tlimits',   [ref_pct_full(1), ref_pct_full(end)], ...
    'freqs',     freqs_in, ...
    'nfreqs',    nfreqs, ...
    'ntimesout', nTout, ...
    'freqscale', params.freqscale, ...
    'cycles',    cycles);

nTimesOut = length(percent_grid);

if params.verbose
    fprintf('  Output grid: %d frequencies, %d %% time points\n', ...
        length(freqs), nTimesOut);
end

% Allocate output arrays now that we know nTimesOut and length(freqs)
TF_X = nan(length(freqs), nTimesOut, nEpochs);
if Y_is_3d
    TF_Y = nan(length(freqs), nTimesOut, nEpochs, nMuscles);
else
    TF_Y = nan(length(freqs), nTimesOut, nEpochs);
end

% Allocate middle-event percent output (only used if track_mid_event is true)
if track_mid_event
    mid_event_pct = nan(nEpochs, 1);
else
    mid_event_pct = [];
end

% Store reference epoch TF result
TF_X(:, :, max_idx) = tf_ref(:, :, 1);
if Y_is_3d
    TF_Y(:, :, max_idx, :) = permute(tf_ref(:, :, 2:end), [1 2 4 3]);
else
    TF_Y(:, :, max_idx) = tf_ref(:, :, 2);
end

% Compute middle-event percent for the reference epoch (no upsampling)
if track_mid_event
    e_mid_ref = timewarp_indices(max_idx, 2);
    if ~isnan(e_mid_ref) && e_mid_ref > 0
        mid_pct_ref = (e_mid_ref - e_on_ref) * delta_pct;
        [~, k] = min(abs(percent_grid - mid_pct_ref));
        mid_event_pct(max_idx) = percent_grid(k);
    end
end

% =========================================================================
% STEP 4 — Per-epoch loop: upsample, build per-epoch % axis, run timefreq
%          with 'timesout' = percent_grid so every epoch's TF matches.
%          The reference epoch (max_idx) is skipped — already computed above.
%
% Speed tricks:
%  - Cache (p, q) per unique upsampling ratio so we do not call rat() every
%    epoch. Also cache resample's filter design via persistent storage in
%    resample itself (MATLAB does this internally for repeated p, q pairs).
%  - Stack X (1 row) and all Y muscles (nMuscles rows) and call timefreq
%    on the multi-row matrix in a single call per epoch. timefreq accepts
%    [time x channels] input and returns [freq x time x channels].
%  - Wrap each epoch in try/catch so a single bad epoch does not abort the
%    whole run. Also pre-check that the returned TF has the expected size.
% =========================================================================

% Pre-compute (p, q) per epoch (vectorised and cheap)
P_all = zeros(nEpochs, 1);
Q_all = zeros(nEpochs, 1);
for ep = 1:nEpochs
    if SOI_lengths(ep) <= 0 || isnan(SOI_lengths(ep))
        continue;
    end
    [P_all(ep), Q_all(ep)] = rat(fs_eff_ratio(ep), params.resample_tol);
end

skipped_epochs = false(nEpochs, 1);
skip_reasons   = strings(nEpochs, 1);

for ep = 1:nEpochs

    % Skip the reference epoch — already computed in STEP 3
    if ep == max_idx
        continue;
    end

    if params.verbose && rem(ep, 100) == 0
        fprintf('  Processing epoch %d / %d...\n', ep, nEpochs);
    end

    try
        % --- Validate event indices ---
        e_on  = timewarp_indices(ep, soi_evt(1));
        e_off = timewarp_indices(ep, soi_evt(2));
        if isnan(e_on) || isnan(e_off) || e_off <= e_on || P_all(ep) == 0
            skipped_epochs(ep) = true;
            skip_reasons(ep)   = "invalid events";
            continue;
        end

        p      = P_all(ep);
        q      = Q_all(ep);
        fs_new = srate * p / q;

        % --- Build a multi-channel matrix [nSamples x (1 + nMuscles)] ---
        if Y_is_3d
            Y_ep = squeeze(Y_epochs(:, :, ep))';   % [nSamples x nMuscles]
        else
            Y_ep = Y_epochs(:, ep);                % [nSamples x 1]
        end
        all_chans = [double(X_epochs(:, ep)), double(Y_ep)];   % [nSamples x (1+nMuscles)]

        % --- Resample all channels at once (resample handles columns) ---
        all_up = resample(all_chans, p, q);
        L_up   = size(all_up, 1);

        % --- Build per-epoch % axis ---
        e_on_up  = round(e_on  * p / q);
        e_off_up = round(e_off * p / q);
        if e_off_up <= e_on_up || e_on_up < 1 || e_off_up > L_up
            skipped_epochs(ep) = true;
            skip_reasons(ep)   = "event indices out of range after upsampling";
            continue;
        end

        soi_pct_e    = linspace(0, 100, e_off_up - e_on_up + 1);
        delta_pct_e  = soi_pct_e(2) - soi_pct_e(1);
        before_pct_e = -delta_pct_e*(e_on_up - 1) : delta_pct_e : -delta_pct_e;
        after_pct_e  = 100 + delta_pct_e : delta_pct_e : ...
                       100 + (L_up - e_off_up) * delta_pct_e;
        full_pct_e   = [before_pct_e, soi_pct_e, after_pct_e];

        % --- Compute middle-event percent (snapped to percent_grid) ---
        if track_mid_event
            e_mid = timewarp_indices(ep, 2);
            if ~isnan(e_mid) && e_mid > 0
                e_mid_up    = round(e_mid * p / q);
                mid_pct_raw = (e_mid_up - e_on_up) * delta_pct_e;
                [~, k]      = min(abs(percent_grid - mid_pct_raw));
                mid_event_pct(ep) = percent_grid(k);
            end
        end

        % --- timefreq on all channels at once ---
        % timefreq accepts [time x channels] and returns [freq x time x channels]
        [tf_all, ~, ~] = timefreq( ...
            all_up, fs_new, ...
            'tlimits',   [full_pct_e(1), full_pct_e(end)], ...
            'timesout',  percent_grid, ...
            'freqs',     freqs_in, ...
            'nfreqs',    nfreqs, ...
            'freqscale', params.freqscale, ...
            'cycles',    cycles);

        % --- Size check: timefreq sometimes drops 1 time point at edges ---
        if size(tf_all, 2) ~= nTimesOut
            skipped_epochs(ep) = true;
            skip_reasons(ep)   = sprintf( ...
                "timefreq returned %d time points (expected %d)", ...
                size(tf_all, 2), nTimesOut);
            continue;
        end

        % --- Unpack: first channel is X, rest are Y muscles ---
        TF_X(:, :, ep) = tf_all(:, :, 1);
        if Y_is_3d
            % tf_all(:,:,2:end) is [freq x time x nMuscles]
            % TF_Y target slice for this epoch is [freq x time x nMuscles]
            TF_Y(:, :, ep, :) = permute(tf_all(:, :, 2:end), [1 2 4 3]);
        else
            TF_Y(:, :, ep) = tf_all(:, :, 2);
        end

    catch ME
        skipped_epochs(ep) = true;
        skip_reasons(ep)   = string(ME.message);
        if params.verbose
            fprintf('  Epoch %d skipped: %s\n', ep, ME.message);
        end
    end
end

% --- Skip summary ---
nSkipped = sum(skipped_epochs);
if nSkipped > 0 && params.verbose
    fprintf('\n  %d epochs skipped:\n', nSkipped);
    [u_reasons, ~, ic] = unique(skip_reasons(skipped_epochs));
    counts = accumarray(ic, 1);
    for r = 1:length(u_reasons)
        fprintf('    %d epoch(s): %s\n', counts(r), u_reasons(r));
    end
end

if params.verbose
    nValid = sum(~isnan(squeeze(TF_X(1, 1, :))));
    fprintf('  Done. %d / %d epochs processed successfully.\n', nValid, nEpochs);
    fprintf('  TF_X size: [%s]\n', num2str(size(TF_X)));
    if Y_is_3d
        fprintf('  TF_Y size: [%s]  (last dim = muscles)\n', num2str(size(TF_Y)));
    else
        fprintf('  TF_Y size: [%s]\n', num2str(size(TF_Y)));
    end
end

end