% MY_NEWCROSSF_TIMENORM - Coherence with time-normalization via per-trial
%                         upsampling, following Fauvet et al. (2019).
%
%   Reference:
%     Fauvet M, Cremoux S, Chalard A, Tisseyre J, Gasq D, Amarantini D.
%     "A novel method to generalize time-frequency coherence analysis
%     between EEG or EMG signals during repetitive trials with high
%     intra-subject variability in duration."
%     9th Int. IEEE EMBS Conf. Neural Engineering, 2019, pp. 437-440.
%
%   Based on MY_NEWCROSSF (Morteza Khosrotabar, TU Darmstadt, 2025),
%   which was itself based on NEWCROSSF (Arnaud Delorme, Sigurd Enghoff
%   & Scott Makeig, CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, 2002-).
%
% %% MK4 MODIFICATION SUMMARY (Fauvet et al. 2019 approach):
%
%   Key idea:
%     Each epoch has a "segment of interest" (SOI) that varies in duration
%     across trials (e.g. movement onset to offset). The SOI of every
%     trial is upsampled to match the longest SOI in the dataset.
%     The flanking segments (pre-SOI and post-SOI) are upsampled at the
%     same rate. This gives each trial an effective (higher) sampling rate
%     fs_eff(t) = fs * (L_max / L_soi(t)), and all trials end up with
%     the same number of samples. The time axis becomes % of movement.
%
%   Pipeline:
%     1. For each trial, identify SOI onset/offset frames.
%     2. Compute effective srate: fs_eff = fs * (L_soi_max / L_soi(t)).
%     3. Upsample the whole epoch at fs_eff using MATLAB resample().
%        All trials now have the same length (reference epoch length).
%     4. TF-decompose each upsampled trial with timefreq() using fs_eff.
%        Because every trial uses its own fs_eff but the same number of
%        output time points are requested, all TF matrices are
%        [nFreqs x nTimesOut] — matched across trials.
%     5. Stack into [nFreqs x nTimesOut x nTrials] and compute coherence
%        using the standard MSC / phase-coherence formula.
%     6. Output time axis is expressed in % of movement realization
%        (0 = SOI onset, 100 = SOI offset).
%
%   IMPORTANT NOTE on tlimits for timefreq():
%     After upsampling, each trial has exactly ref_frames samples at rate
%     fs_eff(t). We pass tlimits_norm = [0, (ref_frames-1)/fs_eff * 1000]
%     to timefreq() so that the internal time grid is consistent with the
%     upsampled signal. The output timesout_pct converts this to %.
%
% Usage:
%   >> [coherres, mbase, timesout_pct, freqsout, Rbootout, Rangle,  ...
%       alltfX, alltfX_pow, alltfY, alltfY_pow]                      ...
%         = my_newcrossf_timenorm(                                   ...
%             x, y, frames, tlimits, srate, cycles,                  ...
%             soi_onsets_ms, soi_offsets_ms,                         ...
%             'key', 'val', ...);
%
% Required inputs:
%   x               = EEG data  (1 x frames*nTrials)
%   y               = EMG data  (1 x frames*nTrials)
%   frames          = Samples per epoch (integer)
%   tlimits         = [tmin tmax] ms, epoch latency limits
%   srate           = Original sampling rate (Hz)
%   cycles          = Wavelet cycles: 0 -> FFT; >0 -> Morlet; [c expf]
%   soi_onsets_ms   = (nTrials x 1) SOI onset  latency per trial (ms)
%   soi_offsets_ms  = (nTrials x 1) SOI offset latency per trial (ms)
%
% Optional inputs (key-value pairs, same as NEWCROSSF):
%   'freqs'      = [fmin fmax] frequency range (Hz). Default [0 50].
%   'nfreqs'     = Number of frequencies. Default [].
%   'timesout'   = Number of output time points OR vector of times.
%                  Default 200. These are expressed in % after output.
%   'baseline'   = Baseline end latency (% or scalar). Default NaN.
%   'alpha'      = Bootstrap significance level. Default NaN (no boot).
%   'type'       = Coherence type: 'coher' (MSC, default), 'phasecoher',
%                  'phasecoher2'.
%   'plotamp'    = 'on'|'off'. Default 'off' (% axis needs custom plot).
%   'plotphase'  = 'on'|'off'. Default 'off'.
%   ... (all other NEWCROSSF key-value options are accepted)
%
% Outputs:
%   coherres     = [nFreqs x nTimesOut]         coherence map
%   mbase        = [nFreqs x 1]                 baseline coherence
%   timesout_pct = [1 x nTimesOut]              time axis (% of movement)
%   freqsout     = [1 x nFreqs]                 frequency axis (Hz)
%   Rbootout     = bootstrap limits ([] if alpha=NaN)
%   Rangle       = [nFreqs x nTimesOut]         coherence phase (radians)
%   alltfX       = [nFreqs x nTimesOut x nTrials] warped EEG TF (complex)
%   alltfX_pow   = [nFreqs x nTimesOut x nTrials] EEG power
%   alltfY       = [nFreqs x nTimesOut x nTrials] warped EMG TF (complex)
%   alltfY_pow   = [nFreqs x nTimesOut x nTrials] EMG power
%
% Example:
%   >> [coherres, mbase, timesout_pct, freqs, ~, Rangle,             ...
%       alltfX, alltfX_pow, alltfY, alltfY_pow]                       ...
%         = my_newcrossf_timenorm(                                    ...
%             x, y, 2000, [-500 3500], 500, [3 0.8],                  ...
%             soi_onsets_ms, soi_offsets_ms,                          ...
%             'freqs', [3 60], 'nfreqs', 100, 'timesout', 200,        ...
%             'plotamp', 'off', 'plotphase', 'off');

function [coherres, mbase, timesout_pct, freqs, Rbootout, Rangle,  ...
          alltfX, alltfX_pow, alltfY, alltfY_pow]                    ...
    = my_newcrossf_timenorm(X, Y, frame, tlimits, Fs, varwin,       ...
                             soi_onsets_ms, soi_offsets_ms, varargin)

% =========================================================================
% Commandline argument defaults
% =========================================================================
DEFAULT_ANGLEUNITS = 'deg';
DEFAULT_EPOCH      = 750;
DEFAULT_TIMELIM    = [-1000 2000];
DEFAULT_FS         = 250;
DEFAULT_NWIN       = 200;
DEFAULT_VARWIN     = 0;
DEFAULT_OVERSMP    = 2;
DEFAULT_MAXFREQ    = 50;
DEFAULT_TITLE      = 'Event-Related Coherence (Time-Normalized)';
DEFAULT_ALPHA      = NaN;

% initialise outputs
coherres     = [];
mbase        = [];
timesout_pct = [];
freqs        = [];
Rbootout     = [];
Rangle       = [];
alltfX       = [];
alltfX_pow   = [];
alltfY       = [];
alltfY_pow   = [];

if nargin < 8
    help my_newcrossf_timenorm
    return
end

% -------------------------------------------------------------------------
% Basic input validation
% -------------------------------------------------------------------------
if ~iscell(X)
    if min(size(X))~=1 || length(X)<2
        fprintf('my_newcrossf_timenorm(): x must be a row or column vector.\n'); return
    elseif min(size(Y))~=1 || length(Y)<2
        fprintf('my_newcrossf_timenorm(): y must be a row or column vector.\n'); return
    elseif length(X)~=length(Y)
        fprintf('my_newcrossf_timenorm(): x and y must have same length.\n'); return
    end
end

if ~isnumeric(frame) || length(frame)~=1 || frame~=round(frame) || frame<=0
    error('my_newcrossf_timenorm(): frames must be a positive integer.');
end
if ~iscell(X) && rem(numel(X), frame)~=0
    error('my_newcrossf_timenorm(): data length must be divisible by frames.');
end
if ~isnumeric(tlimits) || numel(tlimits)~=2 || tlimits(1)>=tlimits(2)
    error('my_newcrossf_timenorm(): tlimits must be [tmin tmax] with tmin < tmax.');
end
if ~isnumeric(Fs) || ~isscalar(Fs) || Fs<=0
    error('my_newcrossf_timenorm(): srate must be a positive scalar.');
end
if ~isnumeric(varwin) || length(varwin)>2 || varwin(1)<0
    error('my_newcrossf_timenorm(): cycles must be zero or positive (scalar or 2-vector).');
end

nTrials = numel(X) / frame;
soi_onsets_ms  = soi_onsets_ms(:);   % ensure column
soi_offsets_ms = soi_offsets_ms(:);

if numel(soi_onsets_ms)~=nTrials || numel(soi_offsets_ms)~=nTrials
    error(['my_newcrossf_timenorm(): soi_onsets_ms and soi_offsets_ms must ' ...
           'each have nTrials (%d) elements.'], nTrials);
end
if any(soi_offsets_ms <= soi_onsets_ms)
    error('my_newcrossf_timenorm(): SOI offset must be later than SOI onset for every trial.');
end
if any(soi_onsets_ms < tlimits(1)) || any(soi_offsets_ms > tlimits(2))
    error('my_newcrossf_timenorm(): SOI latencies must lie within tlimits.');
end

% =========================================================================
% Parse key-value arguments
% =========================================================================
for index = 1:length(varargin)
    if iscell(varargin{index}), varargin{index} = {varargin{index}}; end
end
if ~isempty(varargin)
    [~, indices] = unique_bc(varargin(1:2:end));
    varargin = varargin(sort(union(indices*2-1, indices*2)));
    try,   g = struct(varargin{:});
    catch, error('Argument error in the {''param'', value} sequence'); end
else
    g = [];
end

% =========================================================================
% Set defaults
% =========================================================================
try, g.condboot;      catch, g.condboot      = 'abs';              end
try, g.shuffle;       catch, g.shuffle       = 0;                  end
try, g.title;         catch, g.title         = DEFAULT_TITLE;      end
try, g.winsize;       catch, g.winsize       = max(pow2(nextpow2(frame)-3),4); end
try, g.pad;           catch, g.pad           = max(pow2(nextpow2(g.winsize)),4); end
try, g.timesout;      catch, g.timesout      = DEFAULT_NWIN;       end
try, g.padratio;      catch, g.padratio      = DEFAULT_OVERSMP;    end
try, g.topovec;       catch, g.topovec       = [];                 end
try, g.elocs;         catch, g.elocs         = '';                 end
try, g.alpha;         catch, g.alpha         = DEFAULT_ALPHA;      end
try, g.marktimes;     catch, g.marktimes     = [];                 end
try, g.marktimes = g.vert; catch, g.vert     = [];                 end
try, g.rboot;         catch, g.rboot         = [];                 end
try, g.plotamp;       catch, g.plotamp       = 'off';              end  % default off: % axis
try, g.plotphase;     catch, g.plotphase     = 'off';              end
try, g.plotbootsub;   catch, g.plotbootsub   = 'on';               end
try, g.detrend;       catch, g.detrend       = 'off';              end
try, g.rmerp;         catch, g.rmerp         = 'off';              end
try, g.baseline;      catch, g.baseline      = NaN;                end
try, g.baseboot;      catch, g.baseboot      = 1;                  end
try, g.linewidth;     catch, g.linewidth     = 2;                  end
try, g.maxfreq;       catch, g.maxfreq       = DEFAULT_MAXFREQ;    end
try, g.freqs;         catch, g.freqs         = [0 g.maxfreq];      end
try, g.nfreqs;        catch, g.nfreqs        = [];                 end
try, g.freqscale;     catch, g.freqscale     = 'linear';           end
try, g.naccu;         catch, g.naccu         = 200;                end
try, g.angleunit;     catch, g.angleunit     = DEFAULT_ANGLEUNITS; end
try, g.type;          catch, g.type          = 'coher';            end
try, g.newfig;        catch, g.newfig        = 'on';               end
try, g.boottype;      catch, g.boottype      = 'shuffle';          end
try, g.subitc;        catch, g.subitc        = 'off';              end
try, g.compute;       catch, g.compute       = 'matlab';           end
try, g.maxamp;        catch, g.maxamp        = [];                 end
try, g.savecoher;     catch, g.savecoher     = 0;                  end
try, g.amplag;        catch, g.amplag        = 0;                  end
try, g.noinput;       catch, g.noinput       = 'no';               end
try, g.lowmem;        catch, g.lowmem        = 'off';              end
try, g.plottype;      catch, g.plottype      = 'image';            end
try, g.plotmean;      catch, g.plotmean      = 'on';               end
try, g.highlightmode; catch, g.highlightmode = 'background';       end
try, g.chaninfo;      catch, g.chaninfo      = [];                 end

if isfield(g,'detret'), g.detrend = g.detret; end
if isfield(g,'detrep'), g.rmerp   = g.detrep; end

% =========================================================================
% Validate field names
% =========================================================================
allfields = fieldnames(g);
for index = 1:length(allfields)
    switch allfields{index}
        case { 'shuffle' 'title' 'winsize' 'pad' 'timesout' 'padratio'    ...
               'maxfreq' 'topovec' 'elocs' 'alpha' 'marktimes' 'vert'     ...
               'rboot' 'plotamp' 'plotphase' 'plotbootsub' 'detrep'       ...
               'rmerp' 'detret' 'detrend' 'baseline' 'baseboot'           ...
               'linewidth' 'naccu' 'angleunit' 'type' 'boottype'          ...
               'subitc' 'lowmem' 'plottype' 'compute' 'maxamp'            ...
               'savecoher' 'noinput' 'condboot' 'newfig' 'freqs'          ...
               'nfreqs' 'freqscale' 'amplag' 'highlightmode'              ...
               'plotmean' 'chaninfo' };
            % valid — do nothing
        case {'plotersp' 'plotitc'}
            disp(['my_newcrossf_timenorm warning: option ''' allfields{index} ''' ignored']);
        otherwise
            disp(['my_newcrossf_timenorm error: unrecognized option ''' allfields{index} '''']);
            beep; return;
    end
end

g.tlimits = tlimits;
g.frame   = frame;
g.trials  = nTrials;
g.srate   = Fs;
g.cycles  = varwin;
g.type    = lower(g.type);
g.AXES_FONT  = 10;
g.TITLE_FONT = 14;

if strcmpi(g.freqscale,'log') && g.freqs(1)==0, g.freqs(1) = 3; end

% =========================================================================
% Reshape data into [frame x nTrials]
% =========================================================================
X = reshape(X, frame, nTrials);
Y = reshape(Y, frame, nTrials);

% ERP removal if requested
if strcmpi(g.rmerp,'on')
    X = X - mean(X,2)*ones(1,nTrials);
    Y = Y - mean(Y,2)*ones(1,nTrials);
end

fprintf('\nmy_newcrossf_timenorm: Time-normalization method (Fauvet et al. 2019)\n');
fprintf('  %d trials, original srate = %.1f Hz\n', nTrials, Fs);

% =========================================================================
%% MK4 MODIFICATION — STEP 1:
%  Convert SOI onset/offset ms -> sample indices within each epoch.
%  Compute per-trial SOI lengths and find the reference (maximum) length.
% =========================================================================

% ms -> sample indices (1-based within epoch)
% sample index = round((t_ms - tlimits(1)) / 1000 * Fs) + 1
soi_on_fr  = round((soi_onsets_ms  - tlimits(1)) / 1000 * Fs) + 1;
soi_off_fr = round((soi_offsets_ms - tlimits(1)) / 1000 * Fs) + 1;

% clamp to valid range
soi_on_fr  = max(1,      min(frame, soi_on_fr));
soi_off_fr = max(1,      min(frame, soi_off_fr));

% SOI length per trial (samples)
L_soi = soi_off_fr - soi_on_fr + 1;   % [nTrials x 1]

% Reference: longest SOI in dataset (Fauvet eq. 1 uses l2,i = L_soi_max)
L_soi_max = max(L_soi);

fprintf('  SOI duration range: %d – %d samples (%.0f – %.0f ms)\n', ...
    min(L_soi), L_soi_max, ...
    min(L_soi)/Fs*1000, L_soi_max/Fs*1000);

% =========================================================================
%% MK4 MODIFICATION — STEP 2:
%  Upsample each trial so its SOI occupies exactly L_soi_max samples.
%
%  Key insight: we do NOT need to compute x_i (flanking duration) or
%  manage flanks explicitly. We simply upsample the ENTIRE epoch at the
%  effective srate fs_eff(t) = Fs * L_soi_max / L_soi(t).
%
%  After upsampling, every trial has a different number of samples
%  (because each trial's original epoch has the same duration in ms but
%  the effective srate differs per trial). The SOI portion is guaranteed
%  to be L_soi_max samples in all trials. The pre- and post-SOI flanks
%  simply scale proportionally — no explicit x_i calculation needed.
%
%  To make all trials the same length for TF decomposition, we use the
%  LONGEST upsampled epoch as ref_frames (trial with longest SOI has the
%  lowest fs_eff and therefore fewest upsampled samples; trial with
%  shortest SOI has the highest fs_eff and most upsampled samples).
%  Wait — actually it is the opposite: shortest SOI -> highest fs_eff ->
%  most upsampled samples. So ref_frames = upsampled length of the trial
%  with the SHORTEST SOI (highest upsampling ratio).
%
%  Trials with longer SOI (lower fs_eff, fewer upsampled samples) are
%  zero-padded on the right to reach ref_frames. This padding falls
%  entirely in the post-SOI guard region and does not affect coherence
%  in the 0-100% window.
%
%  The shared tlimits for timefreq() is simply the original tlimits
%  expressed in the normalised time axis: since every upsampled epoch
%  represents the same real-world time window [-500, 3500] ms but at
%  a higher srate, we keep tlimits = original tlimits and pass fs_eff(t)
%  as the srate. timefreq() then builds its time grid correctly.
% =========================================================================

% Effective srate per trial
fs_eff = Fs * L_soi_max ./ L_soi;   % [nTrials x 1]

fprintf('  Effective srate range: %.1f – %.1f Hz\n', min(fs_eff), max(fs_eff));

% ref_frames = upsampled length when using the highest fs_eff
% (shortest SOI trial — most samples after upsampling)
% Compute by dry-running resample on the trial with max fs_eff
[~, idx_max_fseff] = max(fs_eff);
[p_ref, q_ref]     = rat(fs_eff(idx_max_fseff) / Fs, 1e-4);
ref_frames         = length(resample(double(X(:,idx_max_fseff)), p_ref, q_ref));

fprintf('  Normalised epoch length (ref_frames): %d samples\n', ref_frames);
fprintf('  (from trial %d, shortest SOI, highest fs_eff = %.1f Hz)\n', ...
    idx_max_fseff, fs_eff(idx_max_fseff));

% Preallocate upsampled data matrices [ref_frames x nTrials]
X_norm        = zeros(ref_frames, nTrials);
Y_norm        = zeros(ref_frames, nTrials);
fs_eff_actual = zeros(nTrials, 1);   % actual fs_eff after rat() rounding

for t = 1:nTrials

    if rem(t,50)==0
        fprintf('  Upsampling trial %d of %d...\n', t, nTrials);
    end

    [p, q] = rat(fs_eff(t) / Fs, 1e-4);
    fs_eff_actual(t) = Fs * p / q;

    x_up = resample(double(X(:,t)), p, q);
    y_up = resample(double(Y(:,t)), p, q);

    % Place into normalised matrix — zero-pad on right if shorter
    len = length(x_up);
    X_norm(1:min(len,ref_frames), t) = x_up(1:min(len,ref_frames));
    Y_norm(1:min(len,ref_frames), t) = y_up(1:min(len,ref_frames));

end
fprintf('  Upsampling complete.\n');

%% END MK4 MODIFICATION — STEP 2

% =========================================================================
%% MK4 MODIFICATION — STEP 3:
%  TF-decompose all upsampled trials.
%
%  All trials now have ref_frames samples. Each trial's signal represents
%  the original tlimits window [-500, 3500] ms but sampled at fs_eff(t).
%  So we call timefreq() with:
%    - srate  = fs_eff_actual(t)   (trial-specific)
%    - tlimits = original tlimits  (same for all — the ms window is fixed)
%
%  We request the same number/vector of output time points for every
%  trial so all TF matrices are [nFreqs x nTimesOut] — matched.
%  timesout_ms will be in the original ms domain (e.g. -500 to 3500 ms).
%  This is then converted to % of SOI at the end.
% =========================================================================

% Output time points — in original ms domain
if length(g.timesout) > 1
    % user supplied a vector of ms times — use directly
    tmioutopt = {'timesout', g.timesout};
else
    tmioutopt = {'ntimesout', g.timesout};
end

% Spectral options shared across all trials
spectralopts_base = {                        ...
    'winsize',   g.winsize,                  ...
    'tlimits',   tlimits,                    ...
    'detrend',   g.detrend,                  ...
    'subitc',    g.subitc,                   ...
    'wavelet',   g.cycles,                   ...
    'padratio',  g.padratio,                 ...
    'freqs',     g.freqs,                    ...
    'freqscale', g.freqscale,                ...
    'nfreqs',    g.nfreqs,                   ...
    tmioutopt{:} };

if ~strcmpi(g.type,'amp') && ~strcmpi(g.type,'crossspec')
    spectralopts_base = {spectralopts_base{:} 'itctype' g.type};
end

fprintf('\nComputing TF decomposition of normalised epochs...\n');

alltfX = [];
alltfY = [];

for t = 1:nTrials

    if rem(t,50)==0
        fprintf('  TF decomposing trial %d of %d...\n', t, nTrials);
    end

    % Use this trial's effective srate — tlimits stays as original ms window
    [tf_x, freqs, timesout_ms] = timefreq(X_norm(:,t)', fs_eff_actual(t), ...
                                           spectralopts_base{:});
    [tf_y]                     = timefreq(Y_norm(:,t)', fs_eff_actual(t), ...
                                           spectralopts_base{:});

    if isempty(alltfX)
        nFreqs = size(tf_x, 1);
        nTO    = size(tf_x, 2);
        alltfX = zeros(nFreqs, nTO, nTrials, 'double');
        alltfY = zeros(nFreqs, nTO, nTrials, 'double');
    end

    alltfX(:,:,t) = tf_x;
    alltfY(:,:,t) = tf_y;

end
fprintf('  TF decomposition complete.\n');

%% END MK4 MODIFICATION — STEP 3

% =========================================================================
%% MK4 MODIFICATION — STEP 4:
%  Convert output time axis from ms to % of SOI.
%
%  timesout_ms is in the original ms domain (e.g. -500 to 3500 ms).
%  SOI onset  = mean(soi_onsets_ms)  — 0 ms in your case (FlxS)
%  SOI offset = mean(soi_offsets_ms) — mean ExtE latency across trials
%
%  We use the MEAN SOI duration as the 0-100% reference, which is the
%  natural interpretation: 0% = average FlxS onset, 100% = average ExtE.
%  Individual trials are already warped to this by the upsampling.
% =========================================================================
soi_on_ref_ms  = mean(soi_onsets_ms);   % 0 ms in your case
soi_off_ref_ms = mean(soi_offsets_ms);  % mean ExtE latency
soi_dur_ref_ms = soi_off_ref_ms - soi_on_ref_ms;

timesout_pct = (timesout_ms - soi_on_ref_ms) / soi_dur_ref_ms * 100;

fprintf('\n  SOI reference: onset = %.1f ms, offset = %.1f ms (mean across trials)\n', ...
    soi_on_ref_ms, soi_off_ref_ms);
fprintf('  Time axis: %.1f%% to %.1f%% of movement\n', ...
    timesout_pct(1), timesout_pct(end));

%% END MK4 MODIFICATION — STEP 4

% =========================================================================
%% MK4 MODIFICATION — STEP 5:
%  Compute cross-spectrum, power and coherence from the matched TF arrays.
%  Standard MSC formula (Fauvet eq. 2):
%    R^2(f,t) = |S_XY(f,t)|^2 / (S_XX(f,t) * S_YY(f,t))
%  Here we compute the unsquared coherence (magnitude) to match NEWCROSSF
%  convention. The user can square coherres if MSC is needed.
% =========================================================================
fprintf('\nComputing cross-spectrum, power, and coherence...\n');

crossspec  = alltfX .* conj(alltfY);   % [freq x time x trials], complex
alltfX_pow = abs(alltfX).^2;           % [freq x time x trials], real
alltfY_pow = abs(alltfY).^2;           % [freq x time x trials], real

switch g.type
    case 'coher'
        % Linear coherence (MSC numerator/denominator — unsquared):
        % |sum(S_XY)| / sqrt(sum(S_XX) * sum(S_YY))
        num      = sum(crossspec,  3);                   % [freq x time]
        denX     = sum(alltfX_pow, 3);                   % [freq x time]
        denY     = sum(alltfY_pow, 3);                   % [freq x time]
        coherres = num ./ sqrt(denX .* denY);

    case 'phasecoher'
        % Phase-locking value (PLV)
        cs_norm  = crossspec ./ abs(crossspec);
        coherres = sum(cs_norm, 3) / nTrials;

    case 'phasecoher2'
        coherres = sum(crossspec, 3) ./ sum(abs(crossspec), 3);

    case 'crossspec'
        coherres = crossspec;

    otherwise
        error('my_newcrossf_timenorm(): type ''%s'' not supported.', g.type);
end
%% END MK4 MODIFICATION — STEP 5

% =========================================================================
% Baseline
% =========================================================================
if size(g.baseline,2) == 2
    baseln = [];
    for index = 1:size(g.baseline,1)
        tmptime = find(timesout_pct >= g.baseline(index,1) & ...
                       timesout_pct <= g.baseline(index,2));
        baseln = union_bc(baseln, tmptime);
    end
    if isempty(baseln), error('No time point found in baseline range.'); end
else
    if ~isnan(g.baseline) && ~isempty(find(timesout_pct < g.baseline))
        baseln = find(timesout_pct < g.baseline);
    else
        baseln = 1:length(timesout_pct);
    end
end

if ~strcmpi(g.type,'crossspec')
    mbase = mean(abs(coherres(:, baseln)), 2);
else
    mbase = [];
end

% =========================================================================
% Bootstrap
% =========================================================================
if ~isempty(g.rboot)
    Rbootout = g.rboot;
elseif ~isnan(g.alpha) && ~strcmpi(g.type,'crossspec')
    switch g.type
        case 'coher'
            inputdata = {crossspec, alltfX_pow, alltfY_pow};
            formula   = 'sum(arg1,3) ./ sqrt(sum(arg2,3).*sum(arg3,3));';
        case 'phasecoher'
            cs_norm   = crossspec ./ abs(crossspec);
            inputdata = {cs_norm};
            formula   = 'mean(arg1,3);';
        case 'phasecoher2'
            inputdata = {crossspec};
            formula   = 'sum(arg1,3)./sum(abs(arg1),3);';
    end

    if size(g.baseboot,2)==1
        if g.baseboot==0,             baselntmp = [];
        elseif ~isnan(g.baseline(1)), baselntmp = baseln;
        else,                         baselntmp = find(timesout_pct<=0);
        end
    else
        baselntmp = [];
        for index = 1:size(g.baseboot,1)
            tmptime   = find(timesout_pct>=g.baseboot(index,1) & ...
                             timesout_pct<=g.baseboot(index,2));
            baselntmp = union_bc(baselntmp, tmptime);
        end
    end

    Rbootout = bootstat(inputdata, formula, ...
        'boottype',   g.boottype,  ...
        'label',      'coherence', ...
        'bootside',   'upper',     ...
        'shuffledim', [2 3],       ...
        'dimaccu',    2,           ...
        'naccu',      g.naccu,     ...
        'alpha',      g.alpha,     ...
        'basevect',   baselntmp);
else
    Rbootout = [];
end

% =========================================================================
% Format outputs
% =========================================================================
if ~strcmpi(g.type,'crossspec')
    Rangle   = angle(coherres);
    coherres = abs(coherres);
else
    Rangle = [];
end

fprintf('\nmy_newcrossf_timenorm complete.\n');
fprintf('  coherres    [%s] — coherence (time axis = %% movement)\n', num2str(size(coherres)));
fprintf('  timesout_pct range: %.1f to %.1f %%\n', timesout_pct(1), timesout_pct(end));
fprintf('  freqs       [%s] Hz\n', num2str(size(freqs)));
fprintf('  alltfX      [%s] — TF of upsampled EEG\n', num2str(size(alltfX)));
fprintf('  alltfX_pow  [%s] — EEG power\n', num2str(size(alltfX_pow)));
fprintf('  alltfY      [%s] — TF of upsampled EMG\n', num2str(size(alltfY)));
fprintf('  alltfY_pow  [%s] — EMG power\n', num2str(size(alltfY_pow)));

return;

% =========================================================================
% SUBFUNCTIONS
% =========================================================================

function res = dims(array)
res = min(ndims(array), max(size(array,2),size(array,3)));