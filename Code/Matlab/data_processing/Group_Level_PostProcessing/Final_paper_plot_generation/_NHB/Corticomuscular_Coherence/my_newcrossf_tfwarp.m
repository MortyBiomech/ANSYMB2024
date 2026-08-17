% MY_NEWCROSSF_TFWARP - Coherence with timewarp applied to per-trial TF
%                       coefficients BEFORE cross-spectrum computation.
%
%   Based on MY_NEWCROSSF (Morteza Khosrotabar, TU Darmstadt, 2026),
%   which was itself based on NEWCROSSF (Arnaud Delorme, Sigurd Enghoff
%   & Scott Makeig, CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, 2002-).
%
% %% MK3 MODIFICATION SUMMARY:
%
%   Pipeline comparison:
%
%   MY_NEWCROSSF (spectral-domain warping of derived quantities):
%     X, Y  ->  TF decompose (all trials)
%           ->  cross-spectrum  = X_tf .* conj(Y_tf)   [freq x time x trials]
%           ->  power_X         = |X_tf|^2              [freq x time x trials]
%           ->  power_Y         = |Y_tf|^2              [freq x time x trials]
%           ->  WARP cross-spectrum, power_X, power_Y trial by trial
%           ->  coherence from warped quantities
%
%   THIS FUNCTION (MY_NEWCROSSF_TFWARP) — warp TF coefficients first:
%     X, Y  ->  TF decompose each trial individually  ->  X_tf, Y_tf
%           ->  WARP X_tf trial by trial (amplitude via timewarp matrix,
%               phase via angtimewarp, per frequency)
%           ->  WARP Y_tf trial by trial (same warp per trial)
%           ->  cross-spectrum  = X_tf_warped .* conj(Y_tf_warped)
%           ->  power_X         = |X_tf_warped|^2
%           ->  power_Y         = |Y_tf_warped|^2
%           ->  coherence from these quantities
%
%   All modifications marked with %% MK3 MODIFICATION comments.
%
% Usage:
%   >> [coherres, mbase, timesout, freqsout, Rbootout, Rangle,   ...
%       alltfX_warped, alltfX_pow,                                ...
%       alltfY_warped, alltfY_pow,                                ...
%       crossspec]                                                ...
%         = my_newcrossf_tfwarp(x, y, frames, tlimits, srate,    ...
%                               cycles, 'key', 'val', ...);
%
% Required inputs: (identical to NEWCROSSF)
%   x            = First single-channel data (1, frames*nepochs) — EEG
%   y            = Second single-channel data (1, frames*nepochs) — EMG
%   frames       = Frames per epoch
%   tlimits      = [mintime maxtime] (ms) epoch latency limits
%   srate        = Data sampling rate (Hz)
%   cycles       = 0 -> FFTs; >0 -> wavelet cycles; [cycles expfactor]
%
% Timewarp optional inputs:
%   'timewarp'   = [nTrials x nEvents] per-trial event latencies in ms.
%                  Typically EEG.timewarp.latencies(epochs_to_keep,:).
%   'timewarpms' = [1 x nEvents] target latencies in ms. Typically
%                  EEG.timewarp.warpto (group median). If empty, the
%                  column-wise median of 'timewarp' is used.
%   'timewarpidx'= [vector] indices of events to mark with vertical lines.
%
% All other optional inputs identical to NEWCROSSF / MY_NEWCROSSF.
%
% Outputs:
%   coherres      = [freq x time] coherence (MSC / phase coh. etc.)
%   mbase         = [freq x 1]   mean baseline coherence per frequency
%   timesout      = [1 x time]   output time vector (ms)
%   freqsout      = [1 x freq]   output frequency vector (Hz)
%   Rbootout      = bootstrap significance limits ([] if alpha=NaN)
%   Rangle        = [freq x time] coherence phase angle (radians)
%   alltfX_warped = [freq x time x trials] warped complex TF of EEG
%   alltfX_pow    = [freq x time x trials] |alltfX_warped|^2
%   alltfY_warped = [freq x time x trials] warped complex TF of EMG
%   alltfY_pow    = [freq x time x trials] |alltfY_warped|^2
%   crossspec     = [freq x time x trials] per-trial cross-spectrum
%                   (formed from warped TF coefficients)
%
% Example:
%   >> [coherres, mbase, timesout, freqs, ~, Rangle,              ...
%       alltfX_warped, alltfX_pow,                                 ...
%       alltfY_warped, alltfY_pow, crossspec]                      ...
%         = my_newcrossf_tfwarp(                                   ...
%           x, y, 2000, [-500 3500], 500, [3 0.8],                ...
%           'freqs', [3 60], 'nfreqs', 100,                        ...
%           'timesout', ersp_times,                                ...
%           'timewarp',   EEG.timewarp.latencies(epochs_to_keep,:),...
%           'timewarpms', EEG.timewarp.warpto,                     ...
%           'plotamp', 'off', 'plotphase', 'off');

function [coherres, mbase, timesout, freqs, Rbootout, Rangle,    ...
          alltfX_warped, alltfX_pow, alltfY_warped, alltfY_pow,  ...
          crossspec]                                             ...
          = my_newcrossf_tfwarp(X, Y, frame, tlimits, Fs, varwin, varargin)

% =========================================================================
% Commandline argument defaults (identical to newcrossf)
% =========================================================================
DEFAULT_ANGLEUNITS = 'deg';
DEFAULT_EPOCH      = 750;
DEFAULT_TIMELIM    = [-1000 2000];
DEFAULT_FS         = 250;
DEFAULT_NWIN       = 200;
DEFAULT_VARWIN     = 0;
DEFAULT_OVERSMP    = 2;
DEFAULT_MAXFREQ    = 50;
DEFAULT_TITLE      = 'Event-Related Coherence';
DEFAULT_ALPHA      = NaN;

% initialise outputs
coherres      = [];
mbase         = [];
timesout      = [];
freqs         = [];
Rbootout      = [];
Rangle        = [];
alltfX_warped = [];
alltfX_pow    = [];
alltfY_warped = [];
alltfY_pow    = [];
crossspec     = [];

if (nargin < 2)
    help my_newcrossf_tfwarp
    return
end

if ~iscell(X)
    if (min(size(X)) ~= 1 || length(X) < 2)
        fprintf('my_newcrossf_tfwarp(): x must be a row or column vector.\n'); return
    elseif (min(size(Y)) ~= 1 || length(Y) < 2)
        fprintf('my_newcrossf_tfwarp(): y must be a row or column vector.\n'); return
    elseif (length(X) ~= length(Y))
        fprintf('my_newcrossf_tfwarp(): x and y must have same length.\n'); return
    end
end

if (nargin < 3),    frame = DEFAULT_EPOCH;
elseif (~isnumeric(frame) || length(frame)~=1 || frame~=round(frame))
    fprintf('my_newcrossf_tfwarp(): frames must be an integer.\n'); return
elseif (frame <= 0)
    fprintf('my_newcrossf_tfwarp(): frames must be positive.\n'); return
elseif ~iscell(X) && (rem(size(X,2),frame)~=0) && (rem(size(X,1),frame)~=0)
    fprintf('my_newcrossf_tfwarp(): data length must be divisible by frames.\n'); return
end

if (nargin < 4),    tlimits = DEFAULT_TIMELIM;
elseif (~isnumeric(tlimits) || sum(size(tlimits))~=3)
    error('my_newcrossf_tfwarp(): tlimits must be a 2-element vector.');
elseif (tlimits(1) >= tlimits(2))
    error('my_newcrossf_tfwarp(): tlimits must be [min max].');
end

if (nargin < 5),    Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) || length(Fs)~=1)
    error('my_newcrossf_tfwarp(): srate must be a number.');
elseif (Fs <= 0)
    error('my_newcrossf_tfwarp(): srate must be positive.');
end

if (nargin < 6),    varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) || length(varwin)>2)
    error('my_newcrossf_tfwarp(): cycles must be a number or (1,2) vector.');
elseif (varwin < 0)
    error('my_newcrossf_tfwarp(): cycles must be zero or positive.');
end

% =========================================================================
% Parse key-value arguments
% =========================================================================
vararginori = varargin;
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
try, g.plotamp;       catch, g.plotamp       = 'on';               end
try, g.plotphase;     catch, g.plotphase     = 'on';               end
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

% timewarp parameters
try, g.timewarp;      catch, g.timewarp      = [];                 end
try, g.timewarpms;    catch, g.timewarpms    = [];                 end
try, g.timewarpidx;   catch, g.timewarpidx   = [];                 end
g.timewarpfr       = {};
g.timeStretchMarks = [];
g.timeStretchRefs  = [];

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
               'plotmean' 'chaninfo'                                       ...
               'timewarp' 'timewarpms' 'timewarpidx'                      ...
               'timewarpfr' 'timeStretchMarks' 'timeStretchRefs' };
            % valid — do nothing
        case {'plotersp' 'plotitc'}
            disp(['my_newcrossf_tfwarp warning: timef option ''' allfields{index} ''' ignored']);
        otherwise
            disp(['my_newcrossf_tfwarp error: unrecognized option ''' allfields{index} '''']);
            beep; return;
    end
end

g.tlimits    = tlimits;
g.frame      = frame;
if ~iscell(X), g.trials = prod(size(X))/g.frame;
else,          g.trials = prod(size(X{1}))/g.frame; end
g.srate      = Fs;
g.cycles     = varwin;
g.type       = lower(g.type);
g.boottype   = lower(g.boottype);
g.rmerp      = lower(g.rmerp);
g.detrend    = lower(g.detrend);
g.plotphase  = lower(g.plotphase);
g.plotbootsub= lower(g.plotbootsub);
g.subitc     = lower(g.subitc);
g.plotamp    = lower(g.plotamp);
g.compute    = lower(g.compute);
g.AXES_FONT  = 10;
g.TITLE_FONT = 14;

if g.trials == 1 && ~strcmpi(g.type,'crossspec')
    disp('Continuous data: switching to crossspectrum'); g.type = 'crossspec';
end
if strcmpi(g.freqscale,'log') && g.freqs(1)==0, g.freqs(1) = 3; end

if ndims(X)==3
    X = reshape(X, size(X,1), size(X,2)*size(X,3));
    Y = reshape(Y, size(Y,1), size(Y,2)*size(Y,3));
end

% =========================================================================
% Convert timewarp ms -> signal frames -> timeStretch fields
% (Identical logic to newtimef / my_newcrossf)
% =========================================================================
if ~isempty(g.timewarp)
    fprintf('\nmy_newcrossf_tfwarp: Timewarp parameters detected.\n');
    fprintf('  Warping will be applied to X_tf and Y_tf before cross-spectrum.\n');

    evntms = g.timewarp;   % [nTrials x nEvents] in ms
    warpfr = round((evntms - g.tlimits(1)) / 1000 * g.srate) + 1;
    g.timewarpfr{1} = warpfr';   % [nEvents x nTrials]

    if ~isempty(g.timewarpms)
        refms  = g.timewarpms;
        reffr  = round((refms - g.tlimits(1)) / 1000 * g.srate) + 1;
        g.timewarpfr{2} = reffr';
        fprintf('  Using supplied timewarpms as warp targets.\n');
    end

    if ~isempty(g.timewarpidx)
        g.timewarpfr{3} = g.timewarpidx;
    end

    g.timeStretchMarks = g.timewarpfr{1};   % [nEvents x nTrials]
    if length(g.timewarpfr) > 1
        g.timeStretchRefs = g.timewarpfr{2};
    else
        g.timeStretchRefs = median(g.timeStretchMarks, 2);
        fprintf('  No timewarpms supplied — using median event latencies.\n');
    end

    % validate bounds
    if max(g.timeStretchMarks(:)) > frame-2 || min(g.timeStretchMarks(:)) < 3
        error('my_newcrossf_tfwarp: Warp event latencies must be inside the epoch.');
    end
    if ~isempty(g.timeStretchRefs)
        if max(g.timeStretchRefs(:)) > frame-2 || min(g.timeStretchRefs(:)) < 3
            error('my_newcrossf_tfwarp: Warp reference latencies must be inside the epoch.');
        end
    end

    % set marker times for plotting (ms)
    if ~isempty(g.timewarpidx), plotidx = g.timewarpidx;
    else, plotidx = 1:size(g.timeStretchMarks,1); end
    g.marktimes = ((g.timeStretchRefs(plotidx)-1)/g.srate + g.tlimits(1)/1000)*1000;
    fprintf('  Warp marker times (ms):'); fprintf(' %g', g.marktimes); fprintf('\n');
else
    fprintf('\nmy_newcrossf_tfwarp: No timewarp provided — running without warping.\n');
end

% =========================================================================
% Shuffle / detrend (identical to newcrossf)
% =========================================================================
if g.shuffle ~= 0
    fprintf('Shuffling x and y trials %d times\n', g.shuffle);
    XX = reshape(X, 1, frame, length(X)/g.frame);
    YY = Y; X = []; Y = [];
    for index = 1:g.shuffle
        XX = shuffle(XX,3);
        X  = [X XX(:,:)];
        Y  = [Y YY];
    end
end

switch g.rmerp
    case 'on'
        X = reshape(X, g.frame, length(X)/g.frame);
        X = X - mean(X,2)*ones(1,length(X(:))/g.frame);
        Y = reshape(Y, g.frame, length(Y)/g.frame);
        Y = Y - mean(Y,2)*ones(1,length(Y(:))/g.frame);
end

fprintf('\nComputing Event-Related Coherence (type: %s)\n', g.type);
fprintf('Based on %d trials\n', g.trials);

% =========================================================================
% Build spectral options for timefreq
% =========================================================================
if length(g.timesout) > 1
    tmioutopt = {'timesout',  g.timesout};
else
    tmioutopt = {'ntimesout', g.timesout};
end

spectraloptions = { tmioutopt{:},                   ...
    'winsize',   g.winsize,                          ...
    'tlimits',   g.tlimits,                          ...
    'detrend',   g.detrend,                          ...
    'subitc',    g.subitc,                           ...
    'wavelet',   g.cycles,                           ...
    'padratio',  g.padratio,                         ...
    'freqs',     g.freqs,                            ...
    'freqscale', g.freqscale,                        ...
    'nfreqs',    g.nfreqs };

if ~strcmpi(g.type,'amp') && ~strcmpi(g.type,'crossspec')
    spectraloptions = {spectraloptions{:} 'itctype' g.type};
end

% =========================================================================
%% MK3 MODIFICATION — STEP 1:
%  TF decompose X and Y across all trials (no warping at this stage).
%  timefreq() returns [freq x time x trials] complex coefficients.
% =========================================================================
fprintf('\nComputing TF decomposition of EEG...\n');
X = reshape(X, g.frame, g.trials);
[alltfX, freqs, timesout] = timefreq(X, g.srate, spectraloptions{:});
% alltfX : [nFreqs x nTimes x nTrials], complex

fprintf('Computing TF decomposition of EMG...\n');
Y = reshape(Y, g.frame, g.trials);
[alltfY] = timefreq(Y, g.srate, spectraloptions{:});
% alltfY : [nFreqs x nTimes x nTrials], complex

nFreqs = size(alltfX, 1);
nTimes = size(alltfX, 2);
nTrials = g.trials;

% =========================================================================
%% MK3 MODIFICATION — STEP 2:
%  Warp the complex TF coefficients of X and Y trial by trial.
%
%  For each trial t:
%    - Build the warp matrix M (piecewise linear, [nTimes x nTimes])
%      from this trial's event frames -> reference frames.
%    - Warp amplitude of X_tf : TSr_X = M * |X_tf|'  (matrix multiply)
%    - Warp phase   of X_tf : TSphi_X via angtimewarp (circular)
%    - Reconstruct X_tf_warped = TSr_X .* exp(i * TSphi_X)
%    - Same for Y_tf.
%
%  The warp matrix and reference positions are IDENTICAL for X and Y
%  within each trial, so both channels are time-aligned to the same grid.
% =========================================================================
if ~isempty(g.timeStretchMarks)

    fprintf('\nWarping TF coefficients of X and Y trial by trial...\n');

    % Map reference event frames (signal samples) -> TF time-grid indices
    indexout = round(eeg_lat2point(timesout, 1, g.srate, g.tlimits, 1E-3));

    timemarks = g.timeStretchMarks;   % [nEvents x nTrials]
    timerefs  = g.timeStretchRefs;    % [nEvents x 1]

    % Convert timerefs (signal frames) to TF-grid positions
    [~, refsPos] = min(abs( ...
        repmat(timerefs,  [1 nTimes]) - ...
        repmat(indexout,  [length(timerefs) 1]) ...
        )');                              % [1 x nEvents]
    refsPos(end+1) = 1;
    refsPos(end+1) = nTimes;
    refsPos = sort(refsPos);

    % Preallocate warped TF arrays
    alltfX_warped = zeros(nFreqs, nTimes, nTrials, 'like', alltfX);
    alltfY_warped = zeros(nFreqs, nTimes, nTrials, 'like', alltfY);

    for t = 1:nTrials

        if rem(t, 50) == 0
            fprintf('  Warping trial %d of %d...\n', t, nTrials);
        end

        % Convert this trial's event frames to TF-grid positions
        [~, marksPos] = min(abs( ...
            repmat(timemarks(:,t), [1 nTimes]) - ...
            repmat(indexout, [size(timemarks,1) 1]) ...
            )');                           % [1 x nEvents]
        marksPos(end+1) = 1;
        marksPos(end+1) = nTimes;
        marksPos = sort(marksPos);

        % Build piecewise-linear warp matrix [nTimes x nTimes]
        % (shared for X and Y — same trial, same warp)
        M = timewarp(marksPos, refsPos);

        % --- Warp X_tf for this trial ---
        Xt    = alltfX(:,:,t);          % [nFreqs x nTimes], complex
        Xt_r   = abs(Xt);               % amplitude
        Xt_phi = angle(Xt);             % phase

        % Warp amplitude via matrix multiply (linear interp in time)
        Xt_r_w = transpose(M * Xt_r'); % [nFreqs x nTimes]

        % Warp phase via circular interpolation, per frequency
        Xt_phi_w = zeros(nFreqs, nTimes);
        for fi = 1:nFreqs
            Xt_phi_w(fi,:) = angtimewarp(marksPos, refsPos, Xt_phi(fi,:));
        end

        % Reconstruct warped complex TF for X
        alltfX_warped(:,:,t) = Xt_r_w .* exp(1i * Xt_phi_w);

        % --- Warp Y_tf for this trial (same M, same marksPos/refsPos) ---
        Yt    = alltfY(:,:,t);
        Yt_r   = abs(Yt);
        Yt_phi = angle(Yt);

        Yt_r_w = transpose(M * Yt_r');

        Yt_phi_w = zeros(nFreqs, nTimes);
        for fi = 1:nFreqs
            Yt_phi_w(fi,:) = angtimewarp(marksPos, refsPos, Yt_phi(fi,:));
        end

        alltfY_warped(:,:,t) = Yt_r_w .* exp(1i * Yt_phi_w);

    end
    fprintf('Warping complete.\n');

else
    % No warping: warped TF = unwarped TF
    alltfX_warped = alltfX;
    alltfY_warped = alltfY;
    fprintf('No warping applied.\n');
end
%% END MK3 MODIFICATION — STEP 2

% =========================================================================
%% MK3 MODIFICATION — STEP 3:
%  Compute cross-spectrum and power from the warped TF coefficients.
%  This is the standard NEWCROSSF formula, but applied to warped data.
% =========================================================================
fprintf('\nComputing cross-spectrum and power from warped TF coefficients...\n');

crossspec  = alltfX_warped .* conj(alltfY_warped);  % [freq x time x trials]
alltfX_pow = abs(alltfX_warped).^2;                 % [freq x time x trials]
alltfY_pow = abs(alltfY_warped).^2;                 % [freq x time x trials]

% Compute coherence
fprintf('Computing coherence...\n');
switch g.type
    case 'coher'
        % Linear coherence (MSC):
        % C(f,t) = sum_t[ X_w(f,t) * conj(Y_w(f,t)) ]
        %          / sqrt( sum_t|X_w|^2 * sum_t|Y_w|^2 )
        num      = sum(crossspec,  3);                      % [freq x time]
        denX     = sum(alltfX_pow, 3);                      % [freq x time]
        denY     = sum(alltfY_pow, 3);                      % [freq x time]
        coherres = num ./ sqrt(denX .* denY);               % [freq x time]

    case 'phasecoher'
        % Phase-locking value (PLV):
        % normalize each trial's cross-spectrum to unit amplitude first
        cs_norm  = crossspec ./ abs(crossspec);             % [freq x time x trials]
        coherres = sum(cs_norm, 3) / nTrials;               % [freq x time]

    case 'phasecoher2'
        coherres = sum(crossspec, 3) ./ ...
                   sum(abs(crossspec), 3);                  % [freq x time]

    case 'crossspec'
        coherres = crossspec;                               % [freq x time x trials]

    case 'amp'
        % Amplitude correlation uses unwarped amplitudes (no cross-spectrum)
        alltfX_amp = abs(alltfX_warped);
        alltfY_amp = abs(alltfY_warped);
        coherres   = ampcorr(alltfX_amp, alltfY_amp, freqs, timesout, g);
        g.alpha    = NaN;
end
%% END MK3 MODIFICATION — STEP 3

% =========================================================================
% Baseline (identical to newcrossf)
% =========================================================================
if size(g.baseline,2) == 2
    baseln = [];
    for index = 1:size(g.baseline,1)
        tmptime = find(timesout >= g.baseline(index,1) & ...
                       timesout <= g.baseline(index,2));
        baseln = union_bc(baseln, tmptime);
    end
    if isempty(baseln), error('No point found in baseline'); end
else
    if ~isempty(find(timesout < g.baseline))
        baseln = find(timesout < g.baseline);
    else
        baseln = 1:length(timesout);
    end
end

if ~strcmpi(g.type,'crossspec') && ~strcmpi(g.type,'amp')
    mbase = mean(abs(coherres(:, baseln)'));
else
    mbase = [];
end

% =========================================================================
% Bootstrap (identical to newcrossf logic, operating on warped quantities)
% =========================================================================
if ~isempty(g.rboot)
    Rbootout = g.rboot;
elseif ~isnan(g.alpha) && ~strcmpi(g.type,'crossspec') && ~strcmpi(g.type,'amp')
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
        else,                         baselntmp = find(timesout<=0);
        end
    else
        baselntmp = [];
        for index = 1:size(g.baseboot,1)
            tmptime   = find(timesout>=g.baseboot(index,1) & ...
                             timesout<=g.baseboot(index,2));
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
% Plot (identical to newcrossf)
% =========================================================================
if strcmpi(g.plotamp,'on') || strcmpi(g.plotphase,'on')
    if ~strcmpi(g.type,'crossspec') && ~strcmpi(g.type,'amp')
        if strcmpi(g.plottype,'image')
            plotall(coherres, Rbootout, timesout, freqs, mbase, g);
        else
            plotallcurves(coherres, Rbootout, timesout, freqs, mbase, g);
        end
    end
end

% =========================================================================
% Format outputs
% =========================================================================
if ~strcmpi(g.type,'crossspec') && ~strcmpi(g.type,'amp')
    Rangle   = angle(coherres);
    coherres = abs(coherres);
else
    Rangle   = [];
end

fprintf('\nmy_newcrossf_tfwarp complete.\n');
fprintf('  coherres      [%s] — coherence from warped TF coefficients\n', num2str(size(coherres)));
fprintf('  crossspec     [%s] — per-trial cross-spectrum (warped)\n',      num2str(size(crossspec)));
fprintf('  alltfX_warped [%s] — warped complex TF of EEG\n',              num2str(size(alltfX_warped)));
fprintf('  alltfX_pow    [%s] — |alltfX_warped|^2\n',                     num2str(size(alltfX_pow)));
fprintf('  alltfY_warped [%s] — warped complex TF of EMG\n',              num2str(size(alltfY_warped)));
fprintf('  alltfY_pow    [%s] — |alltfY_warped|^2\n',                     num2str(size(alltfY_pow)));

return;

% =========================================================================
% SUBFUNCTIONS (identical to newcrossf / my_newcrossf)
% =========================================================================

function [coherres, lagmap] = ampcorr(alltfX, alltfY, freqs, timesout, g)
coherres = zeros(length(freqs), length(timesout), length(g.amplag));
alpha    = zeros(length(freqs), length(timesout), length(g.amplag));
countlag = 1;
for lag = g.amplag
    fprintf('Computing %d point lag amplitude correlation...\n', lag);
    for i1 = 1:length(freqs)
        for i2 = max(1,1-lag):min(length(timesout)-lag, length(timesout))
            if ~isnan(g.alpha)
                [tmp1,tmp2] = corrcoef(squeeze(alltfX(i1,i2,:)), ...
                                       squeeze(alltfY(i1,i2+lag,:)));
                coherres(i1,i2,countlag) = tmp1(1,2);
                alpha(i1,i2,countlag)    = tmp2(1,2);
            else
                tmp1 = corrcoef(squeeze(alltfX(i1,i2,:)), ...
                                squeeze(alltfY(i1,i2+lag,:)));
                coherres(i1,i2,countlag) = tmp1(1,2);
            end
        end
    end
    countlag = countlag+1;
end
if length(g.amplag)>1
    [coherres,lagmap] = max(coherres,[],3);
    dimsize = length(freqs)*length(timesout);
    alpha   = reshape(alpha((lagmap(:)-1)*dimsize+[1:dimsize]'), ...
                      length(freqs), length(timesout));
    lagmap  = g.amplag(lagmap);
    coherres = coherres.*exp(1i*lagmap/max(abs(g.amplag)));
else
    lagmap = [];
end
if ~isnan(g.alpha)
    tmpind = find(alpha(:)>g.alpha);
    coherres(tmpind) = 0;
end

function plotall(R, Rboot, times, freqs, mbase, g)
switch lower(g.plotphase)
    case 'on'
        switch lower(g.plotamp)
            case 'on',  ordinate1=0.67; ordinate2=0.1; height=0.33; g.plot=1;
            case 'off', ordinate2=0.1;  height=0.9;   g.plot=1;
        end
    case 'off'
        ordinate1=0.1; height=0.9;
        switch lower(g.plotamp)
            case 'on',  ordinate1=0.1; height=0.9; g.plot=1;
            case 'off', g.plot=0;
        end
end
Rangle = angle(R);
if ~isreal(R)
    R=abs(R); Rraw=R; setylim=1;
    if ~isnan(g.baseline), R=R-repmat(mbase',[1 g.timesout]); end
else
    Rraw=R; setylim=0;
end
if g.plot
    set(gcf,'DefaultAxesFontSize',g.AXES_FONT)
    try, icadefs; colormap(feval(DEFAULT_COLORMAP,256)); catch, end
    pos=get(gca,'position'); q=[pos(1) pos(2) 0 0]; s=[pos(3) pos(4) pos(3) pos(4)];
    axis('off')
end
switch lower(g.plotamp)
    case 'on'
        RR=R;
        if ~isnan(g.alpha)
            switch dims(Rboot)
                case 3, RR(find(RR>Rboot(:,:,1)&RR<Rboot(:,:,2)))=0;
                case 2, RR(find(RR<Rboot))=0;
                case 1, RR(find(RR<repmat(Rboot(:),[1 size(RR,2)])))=0;
            end
        end
        h(6)=axes('Units','Normalized','Position',[.1 ordinate1 .8 height].*s+q);
        map=hsv(300); map=flipud([map(251:end,:);map(1:250,:)]); map(151,:)=map(151,:)*0.9;
        colormap(map);
        if ~strcmpi(g.freqscale,'log')
            try, imagesc(times,freqs,RR,max(max(RR))*[-1 1]);
            catch, imagesc(times,freqs,RR,[-1 1]); end
        else
            try, imagesclogy(times,freqs,RR,max(max(RR))*[-1 1]);
            catch, imagesclogy(times,freqs,RR,[-1 1]); end
        end
        set(gca,'ydir','norm');
        if ~isempty(g.maxamp), caxis([-g.maxamp g.maxamp]); end
        tmpscale=caxis;
        hold on
        plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth)
        for i=1:length(g.marktimes)
            plot([g.marktimes(i) g.marktimes(i)],[0 freqs(end)],'--m','LineWidth',g.linewidth);
        end
        hold off
        set(h(6),'YTickLabel',[],'YTick',[])
        set(h(6),'XTickLabel',[],'XTick',[])
        h(8)=axes('Position',[.95 ordinate1 .05 height].*s+q);
        if setylim, cbar(h(8),151:300,[0 tmpscale(2)]);
        else,       cbar(h(8),1:300,[-tmpscale(2) tmpscale(2)]); end
        h(10)=axes('Units','Normalized','Position',[.1 ordinate1-0.1 .8 .1].*s+q);
        Emax=max(R); Emin=min(R);
        plot(times,Emin,times,Emax,'LineWidth',g.linewidth); hold on;
        plot([times(1) times(end)],[0 0],'LineWidth',0.7);
        plot([0 0],[-500 500],'--m','LineWidth',g.linewidth);
        for i=1:length(g.marktimes)
            plot([g.marktimes(i) g.marktimes(i)],[-500 500],'--m','LineWidth',g.linewidth);
        end
        if ~isnan(g.alpha) && dims(Rboot)>1
            switch dims(Rboot)
                case 2, plot(times,mean(Rboot(:,:),1),'g','LineWidth',g.linewidth);
                        plot(times,mean(Rboot(:,:),1),'k:','LineWidth',g.linewidth);
                case 3, plot(times,mean(Rboot(:,:,1),1),'g','LineWidth',g.linewidth);
                        plot(times,mean(Rboot(:,:,1),1),'k:','LineWidth',g.linewidth);
                        plot(times,mean(Rboot(:,:,2),1),'g','LineWidth',g.linewidth);
                        plot(times,mean(Rboot(:,:,2),1),'k:','LineWidth',g.linewidth);
            end
            axis([min(times) max(times) 0 max([Emax(:)' Rboot(:)'])*1.2])
        else
            axis([min(times) max(times) 0 max(Emax)*1.2])
        end
        tick=get(h(10),'YTick'); set(h(10),'YTick',[tick(1);tick(end)])
        set(h(10),'YAxisLocation','right')
        xlabel('Time (ms)'); ylabel('coh.')
        h(11)=axes('Units','Normalized','Position',[0 ordinate1 .1 height].*s+q);
        E=abs(mbase);
        if ~strcmpi(g.freqscale,'log'), plot(freqs,E,'b','LineWidth',g.linewidth);
        else, semilogx(freqs,E,'b','LineWidth',g.linewidth); end
        if ~isnan(g.alpha)
            hold on
            if ~strcmpi(g.freqscale,'log')
                switch dims(Rboot)
                    case 1, plot(freqs,Rboot(:),'g','LineWidth',g.linewidth);
                            plot(freqs,Rboot(:),'k:','LineWidth',g.linewidth);
                    case 2, plot(freqs,mean(Rboot(:,:),2),'g','LineWidth',g.linewidth);
                            plot(freqs,mean(Rboot(:,:),2),'k:','LineWidth',g.linewidth);
                end
            end
            if ~isnan(max(E)), axis([freqs(1) freqs(end) 0 max([E Rboot(:)'])*1.2]); end
        else
            if ~isnan(max(E)), axis([freqs(1) freqs(end) 0 max(E)*1.2]); end
        end
        set(gca,'xdir','rev');
        tick=get(h(11),'YTick'); set(h(11),'YTick',[tick(1);tick(end)]);
        set(h(11),'View',[90 90])
        xlabel('Freq. (Hz)'); ylabel('coh.')
end
switch lower(g.plotphase)
    case 'on'
        h(13)=axes('Units','Normalized','Position',[.1 ordinate2 .8 height].*s+q);
        if setylim
            if strcmp(g.angleunit,'ms')
                Rangle=(Rangle/(2*pi)).*repmat(1000./freqs(:)',1,length(times));
                maxangle=max(max(abs(Rangle)));
            elseif strcmpi(g.angleunit,'deg')
                Rangle=Rangle*180/pi; maxangle=180;
            else, maxangle=pi;
            end
            Rangle(find(Rraw==0))=0;
            if ~strcmpi(g.freqscale,'log'), imagesc(times,freqs,Rangle,[-maxangle maxangle]);
            else, imagesclogy(times,freqs,Rangle,[-maxangle maxangle]); end
            hold on
            plot([0 0],[0 freqs(end)],'--m','LineWidth',g.linewidth);
            for i=1:length(g.marktimes)
                plot([g.marktimes(i) g.marktimes(i)],[0 freqs(end)],'--m','LineWidth',g.linewidth);
            end
            set(gca,'ydir','norm');
            ylabel('Freq. (Hz)'); xlabel('Time (ms)')
            h(14)=axes('Position',[.95 ordinate2 .05 height].*s+q);
            cbar(h(14),0,[-maxangle maxangle]);
        else
            axis off; text(0,0.5,'Real values, no angles');
        end
end
if g.plot
    try, icadefs; set(gcf,'color',BACKCOLOR); catch, end
    if length(g.title)>0
        axes('Position',pos,'Visible','Off');
        h(13)=text(-.05,1.01,g.title);
        set(h(13),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',g.TITLE_FONT)
    end
    try, axcopy(gcf); catch, end
end

function plotallcurves(R, Rboot, times, freqs, mbase, g)
Rangle=angle(R);
if ~isreal(R)
    R=abs(R); Rraw=R;
    if ~isnan(g.baseline), R=R-repmat(mbase',[1 g.timesout]); end
else
    Rraw=R; setylim=0;
end
if times(end)>10000, times=times/1000; timeunit='s'; else, timeunit='ms'; end
alllegend={};
if strcmpi(g.plotmean,'on') && freqs(1)~=freqs(end)
    alllegend={[num2str(freqs(1)) '-' num2str(freqs(end)) 'Hz']};
else
    for index=1:length(freqs), alllegend{index}=[num2str(freqs(index)) 'Hz']; end
end
if strcmpi(g.plotamp,'on')
    if strcmpi(g.plotphase,'on'), subplot(2,1,1); end
    if isempty(g.maxamp), g.maxamp=0; end
    plotcurve(times,R,'maskarray',Rboot,'title','Coherence amplitude',...
        'xlabel',['Time (' timeunit ')'],'ylabel','0-1','ylim',g.maxamp,...
        'vert',g.vert,'marktimes',g.marktimes,'legend',alllegend,...
        'linewidth',g.linewidth,'highlightmode',g.highlightmode,'plotmean',g.plotmean);
end
if strcmpi(g.plotphase,'on')
    if strcmpi(g.plotamp,'on'), subplot(2,1,2); end
    plotcurve(times,Rangle/pi*180,'maskarray',Rboot,'val2mask',R,...
        'title','Coherence phase','xlabel',['Time (' timeunit ')'],...
        'ylabel','Angle (deg.)','ylim',[-180 180],...
        'vert',g.vert,'marktimes',g.marktimes,'legend',alllegend,...
        'linewidth',g.linewidth,'highlightmode',g.highlightmode,'plotmean',g.plotmean);
end
if strcmpi(g.plotamp,'on')||strcmpi(g.plotphase,'on')
    try, icadefs; set(gcf,'color',BACKCOLOR); catch, end
    if length(g.title)>0, h(13)=textsc(g.title,'title'); end
    try, axcopy(gcf); catch, end
end

function res = dims(array)
res = min(ndims(array), max(size(array,2),size(array,3)));