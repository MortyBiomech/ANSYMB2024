function [TF_X, TF_Y, time_grid_ms, freqs] = ...
    compute_TF_locked(X_epochs, Y_epochs, params, tlimits_ms)
% COMPUTE_TF_LOCKED
%   Time-frequency decomposition with no time normalisation.
%   All epochs share the original epoch window and the same time axis.
%   EEGLAB's timefreq() is called in batch mode on all epochs at once
%   ([nSamples x nEpochs] input), which is more efficient than a per-epoch loop.
%
%   This implements the "No timewarp, no upsampling, locked to the 1st event"
%   condition from the design diagram.
%
% Inputs:
%   X_epochs   : [nSamples x nEpochs]  — EEG IC activations
%   Y_epochs   : [nMuscles x nSamples x nEpochs]  or  [nSamples x nEpochs]
%                EMG signal(s)
%   params     : struct
%                  .srate       (Hz)
%                  .freqs       ([fmin fmax] Hz)
%                  .nfreqs      (number of frequency bins)
%                  .cycles      (Morlet cycles, e.g. [2 0.8])
%                  .freqscale   ('log' or 'linear')
%                  .ntimesout   (number of output time points)
%   tlimits_ms : [tmin tmax] epoch window in milliseconds
%
% Outputs:
%   TF_X        : [nFreqs x nTimes x nEpochs]           complex EEG TF
%   TF_Y        : [nFreqs x nTimes x nEpochs x nMuscles] complex EMG TF
%                 (or [nFreqs x nTimes x nEpochs] if Y_epochs is 2-D)
%   time_grid_ms: [1 x nTimes] output time axis in milliseconds
%   freqs       : [1 x nFreqs] frequency axis in Hz
%
% Note:
%   The returned time_grid_ms is the EEGLAB timefreq() output time vector.
%   Pass it as the "percent_grid" argument to compute_and_plot_CMC_cluster_perm()
%   (the x-axis will then show milliseconds rather than cycle %).

% --- handle 2-D vs 3-D Y_epochs ---
Y_is_3d = (ndims(Y_epochs) == 3);
if Y_is_3d
    [nMuscles, nSamples_Y, nEpochs_Y] = size(Y_epochs);
else
    [nSamples_Y, nEpochs_Y] = size(Y_epochs);
    nMuscles  = 1;
    Y_epochs  = reshape(Y_epochs, 1, nSamples_Y, nEpochs_Y);
end
[nSamples_X, nEpochs_X] = size(X_epochs);
assert(nEpochs_X == nEpochs_Y,   'X and Y must have the same number of epochs.');
assert(nSamples_X == nSamples_Y, 'X and Y must have the same number of samples.');
nEpochs = nEpochs_X;

% ── Compute TF for X (all epochs in one batch call) ───────────────────────
% timefreq accepts [nSamples x nEpochs] and returns [nFreqs x nTimes x nEpochs].
[TF_X, freqs, time_grid_ms] = timefreq( ...
    double(X_epochs), params.srate, ...
    'tlimits',   tlimits_ms, ...
    'freqs',     params.freqs, ...
    'nfreqs',    params.nfreqs, ...
    'ntimesout', params.ntimesout, ...
    'freqscale', params.freqscale, ...
    'cycles',    params.cycles);
% TF_X : [nFreqs x nTimes x nEpochs]

nF = length(freqs);
nT = length(time_grid_ms);

% ── Compute TF for Y (one muscle at a time to avoid memory blowout) ───────
if Y_is_3d
    TF_Y = nan(nF, nT, nEpochs, nMuscles);
else
    TF_Y = nan(nF, nT, nEpochs);
end

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
    % tf_m : [nFreqs x nTimes x nEpochs]

    if Y_is_3d
        TF_Y(:, :, :, m) = tf_m;
    else
        TF_Y = tf_m;   % single muscle, 3-D output
    end
end

% Ensure time_grid_ms is a row vector
time_grid_ms = time_grid_ms(:)';

end
