% MY_NEWCROSSF - Modified version of NEWCROSSF with timewarp support.
%
%   This function is based on NEWCROSSF (Arnaud Delorme, Sigurd Enghoff &
%   Scott Makeig, CNL/Salk Institute 1998-2001; SCCN/INC/UCSD, 2002-)
%   and adds time-warping capability by passing 'timestretch' parameters
%   to the internal TIMEFREQ calls, mirroring the approach used in
%   NEWTIMEF.
%
% %% MK MODIFICATION SUMMARY (Morteza Khosrotabar, TU Darmstadt, 2026):
%   Three targeted changes were made to the original NEWCROSSF:
%   (1) Added 'timewarp', 'timewarpms', 'timewarpidx' to accepted params.
%   (2) Added ms-to-frames conversion block (same logic as NEWTIMEF).
%   (3) Passed 'timestretch' to both TIMEFREQ calls for X and Y.
%   All modifications are marked with %% MK MODIFICATION comments.
%
% Usage:
%   >> [coh, mcoh, timesout, freqsout, cohboot, cohangles, ...
%          allcoher, alltfX, alltfY] = my_newcrossf(x, y, frames, ...
%          tlimits, srate, cycles, 'key', 'val', ...);
%
% Required inputs: (identical to NEWCROSSF)
%   x       = First single-channel data set (1, frames*nepochs)
%   y       = Second single-channel data set (1, frames*nepochs)
%   frames  = Frames per epoch
%   tlimits = [mintime maxtime] (ms) epoch latency limits
%   srate   = Data sampling rate (Hz)
%   cycles  = 0 -> FFTs; >0 -> wavelet cycles; [cycles expfactor]
%
% %% MK MODIFICATION (1): New optional timewarp inputs
%   'timewarp'   = [matrix] (nTrials x nEvents) matrix of event latencies
%                  (in ms) at which events occur in each trial. These are
%                  the subject-specific event latencies from
%                  EEG.timewarp.latencies. A value of 0 (placeholder) means
%                  no warping is applied.
%   'timewarpms' = [vector] (1 x nEvents) target latencies (in ms) to warp
%                  events to. These are the group median latencies from
%                  EEG.timewarp.warpto. If empty, median of 'timewarp'
%                  columns is used.
%   'timewarpidx'= [vector] indices of events to mark with vertical lines.
%
% All other optional inputs are identical to NEWCROSSF.
% See >> help newcrossf for full parameter documentation.
%
% Outputs: (identical to NEWCROSSF)
%   coh         = (nfreqs x timesout) coherence magnitude matrix
%   mcoh        = vector of mean baseline coherence at each frequency
%   timesout    = vector of output latencies (ms)
%   freqsout    = vector of frequency bin centers (Hz)
%   cohboot     = bootstrap coherence significance limits
%   cohangles   = (nfreqs x timesout) coherence phase angles (radians)
%   allcoher    = single-trial coherence
%   alltfX      = single-trial TF decomposition of X
%   alltfY      = single-trial TF decomposition of Y
%
% Example call with timewarping:
%   >> [coh, mcoh, timesout, freqs] = my_newcrossf(x, y, 2000, ...
%        [-500 3500], 500, [3 0.8], ...
%        'freqscale', 'log', 'nfreqs', 100, ...
%        'timewarp',   EEG.timewarp.latencies, ...
%        'timewarpms', EEG.timewarp.warpto);
%
% Authors: Based on NEWCROSSF by Arnaud Delorme, Sigurd Enghoff &
%          Scott Makeig, CNL/Salk Institute 1998-2001;
%          SCCN/INC/UCSD, La Jolla, 2002-
%          Modified by Morteza Khosrotabar, TU Darmstadt, 2026.

function [R, mbase, timesout, freqs, Rbootout, Rangle, coherresout, alltfX, alltfY] = ...
    my_newcrossf_v1(X, Y, frame, tlimits, Fs, varwin, varargin)

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

if (nargin < 2)
    help my_newcrossf
    return
end
coherresout = [];

if ~iscell(X)
    if (min(size(X)) ~= 1 || length(X) < 2)
        fprintf('my_newcrossf(): x must be a row or column vector.\n');
        return
    elseif (min(size(Y)) ~= 1 || length(Y) < 2)
        fprintf('my_newcrossf(): y must be a row or column vector.\n');
        return
    elseif (length(X) ~= length(Y))
        fprintf('my_newcrossf(): x and y must have same length.\n');
        return
    end
end

if (nargin < 3)
    frame = DEFAULT_EPOCH;
elseif (~isnumeric(frame) || length(frame) ~= 1 || frame ~= round(frame))
    fprintf('my_newcrossf(): Value of frames must be an integer.\n');
    return
elseif (frame <= 0)
    fprintf('my_newcrossf(): Value of frames must be positive.\n');
    return
elseif ~iscell(X) && (rem(size(X,2), frame) ~= 0) && (rem(size(X,1), frame) ~= 0)
    fprintf('my_newcrossf(): Length of data must be divisible by frames.\n');
    return
end

if (nargin < 4)
    tlimits = DEFAULT_TIMELIM;
elseif (~isnumeric(tlimits) || sum(size(tlimits)) ~= 3)
    error('my_newcrossf(): tlimits must be a 2-element vector.');
elseif (tlimits(1) >= tlimits(2))
    error('my_newcrossf(): tlimits interval must be [min, max].');
end

if (nargin < 5)
    Fs = DEFAULT_FS;
elseif (~isnumeric(Fs) || length(Fs) ~= 1)
    error('my_newcrossf(): srate must be a number.');
elseif (Fs <= 0)
    error('my_newcrossf(): srate must be positive.');
end

if (nargin < 6)
    varwin = DEFAULT_VARWIN;
elseif (~isnumeric(varwin) || length(varwin) > 2)
    error('my_newcrossf(): cycles must be a number or (1,2) vector.');
elseif (varwin < 0)
    error('my_newcrossf(): cycles must be zero or positive.');
end

% =========================================================================
% Parse key-value arguments
% =========================================================================
vararginori = varargin;
for index = 1:length(varargin)
    if iscell(varargin{index}), varargin{index} = {varargin{index}}; end
end
if ~isempty(varargin)
    [tmp, indices] = unique_bc(varargin(1:2:end));
    varargin = varargin(sort(union(indices*2-1, indices*2)));
    try,   g = struct(varargin{:});
    catch, error('Argument error in the {''param'', value} sequence'); end
else
    g = [];
end

% =========================================================================
% Set defaults for all parameters
% =========================================================================
try, g.condboot;        catch, g.condboot       = 'abs';          end
try, g.shuffle;         catch, g.shuffle        = 0;              end
try, g.title;           catch, g.title          = DEFAULT_TITLE;  end
try, g.winsize;         catch, g.winsize        = max(pow2(nextpow2(frame)-3), 4); end
try, g.pad;             catch, g.pad            = max(pow2(nextpow2(g.winsize)), 4); end
try, g.timesout;        catch, g.timesout       = DEFAULT_NWIN;   end
try, g.padratio;        catch, g.padratio       = DEFAULT_OVERSMP; end
try, g.topovec;         catch, g.topovec        = [];             end
try, g.elocs;           catch, g.elocs          = '';             end
try, g.alpha;           catch, g.alpha          = DEFAULT_ALPHA;  end
try, g.marktimes;       catch, g.marktimes      = [];             end
try, g.marktimes = g.vert; catch, g.vert        = [];             end
try, g.rboot;           catch, g.rboot          = [];             end
try, g.plotamp;         catch, g.plotamp        = 'on';           end
try, g.plotphase;       catch, g.plotphase      = 'on';           end
try, g.plotbootsub;     catch, g.plotbootsub    = 'on';           end
try, g.detrend;         catch, g.detrend        = 'off';          end
try, g.rmerp;           catch, g.rmerp          = 'off';          end
try, g.baseline;        catch, g.baseline       = NaN;            end
try, g.baseboot;        catch, g.baseboot       = 1;              end
try, g.linewidth;       catch, g.linewidth      = 2;              end
try, g.maxfreq;         catch, g.maxfreq        = DEFAULT_MAXFREQ; end
try, g.freqs;           catch, g.freqs          = [0 g.maxfreq];  end
try, g.nfreqs;          catch, g.nfreqs         = [];             end
try, g.freqscale;       catch, g.freqscale      = 'linear';       end
try, g.naccu;           catch, g.naccu          = 200;            end
try, g.angleunit;       catch, g.angleunit      = DEFAULT_ANGLEUNITS; end
try, g.type;            catch, g.type           = 'phasecoher';   end
try, g.newfig;          catch, g.newfig         = 'on';           end
try, g.boottype;        catch, g.boottype       = 'shuffle';      end
try, g.subitc;          catch, g.subitc         = 'off';          end
try, g.compute;         catch, g.compute        = 'matlab';       end
try, g.maxamp;          catch, g.maxamp         = [];             end
try, g.savecoher;       catch, g.savecoher      = 0;              end
try, g.amplag;          catch, g.amplag         = 0;              end
try, g.noinput;         catch, g.noinput        = 'no';           end
try, g.lowmem;          catch, g.lowmem         = 'off';          end
try, g.plottype;        catch, g.plottype       = 'image';        end
try, g.plotmean;        catch, g.plotmean       = 'on';           end
try, g.highlightmode;   catch, g.highlightmode  = 'background';   end
try, g.chaninfo;        catch, g.chaninfo       = [];             end

%% MK MODIFICATION (1): defaults for new timewarp parameters
try, g.timewarp;        catch, g.timewarp       = [];             end
try, g.timewarpms;      catch, g.timewarpms     = [];             end
try, g.timewarpidx;     catch, g.timewarpidx    = [];             end
% Internal fields (will be populated in conversion block below)
g.timewarpfr        = {};
g.timeStretchMarks  = [];
g.timeStretchRefs   = [];
%% END MK MODIFICATION (1)

if isfield(g, 'detret'), g.detrend = g.detret; end
if isfield(g, 'detrep'), g.rmerp   = g.detrep; end
if ~isnan(g.alpha) && frame == length(X)
    error('Cannot compute significance for continuous data.');
end

% =========================================================================
% Validate all field names — MK MODIFICATION (1): added timewarp fields
% =========================================================================
allfields = fieldnames(g);
for index = 1:length(allfields)
    switch allfields{index}
        case { 'shuffle' 'title' 'winsize' 'pad' 'timesout' 'padratio' ...
               'maxfreq' 'topovec' 'elocs' 'alpha' 'marktimes' 'vert' ...
               'rboot' 'plotamp' 'plotphase' 'plotbootsub' 'detrep' ...
               'rmerp' 'detret' 'detrend' 'baseline' 'baseboot' ...
               'linewidth' 'naccu' 'angleunit' 'type' 'boottype' ...
               'subitc' 'lowmem' 'plottype' 'compute' 'maxamp' ...
               'savecoher' 'noinput' 'condboot' 'newfig' 'freqs' ...
               'nfreqs' 'freqscale' 'amplag' 'highlightmode' ...
               'plotmean' 'chaninfo' ...
               ... %% MK MODIFICATION (1): timewarp fields added to accepted list
               'timewarp' 'timewarpms' 'timewarpidx' ...
               'timewarpfr' 'timeStretchMarks' 'timeStretchRefs' };
            % valid — do nothing
        case {'plotersp' 'plotitc'}
            disp(['my_newcrossf warning: timef option ''' ...
                allfields{index} ''' ignored']);
        otherwise
            disp(['my_newcrossf error: unrecognized option ''' ...
                allfields{index} '''']);
            beep; return;
    end
end

g.tlimits = tlimits;
g.frame   = frame;
if ~iscell(X)
    g.trials = prod(size(X)) / g.frame;
else
    g.trials = prod(size(X{1})) / g.frame;
end
g.srate      = Fs;
g.cycles     = varwin;
g.type       = lower(g.type);
g.boottype   = lower(g.boottype);
g.rmerp      = lower(g.rmerp);
g.detrend    = lower(g.detrend);
g.plotphase  = lower(g.plotphase);
g.plotbootsub = lower(g.plotbootsub);
g.subitc     = lower(g.subitc);
g.plotamp    = lower(g.plotamp);
g.compute    = lower(g.compute);
g.AXES_FONT  = 10;
g.TITLE_FONT = 14;

if g.trials == 1 && ~strcmpi(g.type, 'crossspec')
    disp('Continuous data: switching to crossspectrum');
    g.type = 'crossspec';
end
if strcmpi(g.freqscale, 'log') && g.freqs(1) == 0, g.freqs(1) = 3; end

% reshape 3D inputs
if ndims(X) == 3
    X = reshape(X, size(X,1), size(X,2)*size(X,3));
    Y = reshape(Y, size(Y,1), size(Y,2)*size(Y,3));
end

% update title based on type
if strcmpi(g.title, DEFAULT_TITLE)
    switch g.type
        case 'coher',       g.title = 'Event-Related Coherence';
        case 'phasecoher',  g.title = 'Event-Related Phase Coherence';
        case 'phasecoher2', g.title = 'Event-Related Phase Coherence 2';
        case 'amp',         g.title = 'Event-Related Amplitude Correlation';
        case 'crossspec',   g.title = 'Event-Related Cross-Spectrum';
    end
end

% argument consistency checks
if ~ischar(g.title) && ~iscell(g.title)
    error('Title must be a string or cell array.');
end
if (~isnumeric(g.alpha) || length(g.alpha) ~= 1)
    error('my_newcrossf(): alpha must be a number.');
elseif (round(g.naccu*g.alpha) < 2)
    fprintf('Alpha out of normal range [%g,0.5]\n', 2/g.naccu);
    g.naccu = round(2/g.alpha);
    fprintf('Increasing bootstrap iterations to %d\n', g.naccu);
end
if g.alpha > 0.5 || g.alpha <= 0
    error('Alpha out of allowed range (0.00, 0.5).');
end
switch lower(g.newfig)
    case {'on','off'}, ;
    otherwise, error('newfig must be on or off');
end
switch g.angleunit
    case {'ms','deg','rad'}, ;
    otherwise, error('Angleunit must be deg, rad, or ms');
end
switch g.type
    case {'coher','phasecoher','phasecoher2','amp','crossspec'}, ;
    otherwise, error('Invalid coherence type.');
end
switch g.boottype
    case {'shuffle','shufftrials','rand','randall'}, ;
    otherwise, error('Invalid boottype.');
end
if ~isnumeric(g.shuffle) && ~iscell(g.shuffle)
    error('Shuffle argument must be numeric.');
end
switch g.compute
    case {'matlab','c'}, ;
    otherwise, error('compute must be matlab or c');
end
if ~strcmpi(g.condboot,'abs') && ~strcmpi(g.condboot,'angle') && ...
        ~strcmpi(g.condboot,'complex')
    error('condboot must be abs, angle, or complex.');
end
if g.tlimits(2) - g.tlimits(1) < 30
    disp('my_newcrossf WARNING: time range very small (<30ms). Check units.');
end

% =========================================================================
%% MK MODIFICATION (2): Convert timewarp ms -> frames -> timeStretch params
%  (Mirrors the exact same logic used in newtimef.m)
% =========================================================================
if ~isempty(g.timewarp)
    fprintf('\nmy_newcrossf: Applying time warping to X and Y...\n');

    % Convert event latencies from ms to frame indices
    evntms = g.timewarp;   % [nTrials x nEvents] matrix in ms
    warpfr = round((evntms - g.tlimits(1)) / 1000 * g.srate) + 1;
    g.timewarpfr{1} = warpfr';   % [nEvents x nTrials]

    % Convert target latencies (warpto) from ms to frame indices
    if ~isempty(g.timewarpms)
        refms  = g.timewarpms;   % [1 x nEvents] vector in ms
        reffr  = round((refms - g.tlimits(1)) / 1000 * g.srate) + 1;
        g.timewarpfr{2} = reffr';
    end

    % Store optional plot indices
    if ~isempty(g.timewarpidx)
        g.timewarpfr{3} = g.timewarpidx;
    end

    % Populate timeStretch fields used by timefreq()
    g.timeStretchMarks = g.timewarpfr{1};   % [nEvents x nTrials]
    if length(g.timewarpfr) > 1
        g.timeStretchRefs = g.timewarpfr{2};
        fprintf('my_newcrossf: Using supplied warp target latencies.\n');
    else
        % Fallback: use median of event latencies as reference
        g.timeStretchRefs = median(g.timeStretchMarks, 2);
        fprintf(['my_newcrossf: No timewarpms supplied. ' ...
            'Using median event latencies as warp targets.\n']);
    end

    % Validate frame bounds
    if max(g.timeStretchMarks(:)) > frame-2 || ...
            min(g.timeStretchMarks(:)) < 3
        error(['my_newcrossf: Time warping events must be inside ' ...
            'the epoch (not at boundaries).']);
    end
    if ~isempty(g.timeStretchRefs)
        if max(g.timeStretchRefs(:)) > frame-2 || ...
                min(g.timeStretchRefs(:)) < 3
            error(['my_newcrossf: Warp reference latencies must be ' ...
                'within the epoch boundaries.']);
        end
    end

    % Determine vertical marker times for plotting (in ms)
    if ~isempty(g.timewarpidx)
        plotidx = g.timewarpidx;
    else
        plotidx = 1:size(g.timeStretchMarks, 1);
    end
    g.marktimes = ((g.timeStretchRefs(plotidx) - 1) / g.srate + ...
        g.tlimits(1)/1000) * 1000;
    fprintf('my_newcrossf: Warp marker times (ms): ');
    fprintf('%g ', g.marktimes);
    fprintf('\n');
else
    % No timewarping — clear timeStretch fields
    g.timeStretchMarks = [];
    g.timeStretchRefs  = [];
    fprintf('my_newcrossf: No timewarp provided. Running without warping.\n');
end
%% END MK MODIFICATION (2)

% =========================================================================
% Shuffle trials if requested (identical to newcrossf)
% =========================================================================
if g.shuffle ~= 0
    fprintf('x and y data trials being shuffled %d times\n', g.shuffle);
    XX = reshape(X, 1, frame, length(X)/g.frame);
    YY = Y;
    X  = [];
    Y  = [];
    for index = 1:g.shuffle
        XX = shuffle(XX, 3);
        X  = [X XX(:,:)];
        Y  = [Y YY];
    end
end

% detrend over epochs if requested
switch g.rmerp
    case 'on'
        X = reshape(X, g.frame, length(X)/g.frame);
        X = X - mean(X, 2) * ones(1, length(X(:))/g.frame);
        Y = reshape(Y, g.frame, length(Y)/g.frame);
        Y = Y - mean(Y, 2) * ones(1, length(Y(:))/g.frame);
end

% =========================================================================
% Display info
% =========================================================================
fprintf('\nComputing the Event-Related\n');
switch g.type
    case 'phasecoher',  fprintf('Phase Coherence (ITC) based on %d trials\n', g.trials);
    case 'phasecoher2', fprintf('Phase Coherence 2 (ITC) based on %d trials\n', g.trials);
    case 'coher',       fprintf('Linear Coherence based on %d trials\n', g.trials);
    case 'amp',         fprintf('Amplitude correlation based on %d trials\n', g.trials);
    case 'crossspec',   fprintf('Cross-spectral images based on %d trials\n', g.trials);
end
if ~isnan(g.alpha)
    fprintf('Bootstrap confidence limits: alpha = %g\n', g.alpha);
else
    fprintf('Bootstrap confidence limits will NOT be computed.\n');
end
switch g.plotphase
    case 'on', fprintf(['Coherence angles will be imaged in ', g.angleunit, '\n']);
end

% =========================================================================
% Build spectral options for timefreq calls
% =========================================================================
if length(g.timesout) > 1
    tmioutopt = {'timesout',  g.timesout};
else
    tmioutopt = {'ntimesout', g.timesout};
end

spectraloptions = { tmioutopt{:}, ...
    'winsize',   g.winsize, ...
    'tlimits',   g.tlimits, ...
    'detrend',   g.detrend, ...
    'subitc',    g.subitc, ...
    'wavelet',   g.cycles, ...
    'padratio',  g.padratio, ...
    'freqs',     g.freqs, ...
    'freqscale', g.freqscale, ...
    'nfreqs',    g.nfreqs };

if ~strcmpi(g.type, 'amp') && ~strcmpi(g.type, 'crossspec')
    spectraloptions = {spectraloptions{:} 'itctype' g.type};
end

%% MK MODIFICATION (3): Build timestretch option and pass to BOTH timefreq calls
% This is the key change: timefreq receives the timeStretchMarks and
% timeStretchRefs so that both X and Y are warped identically before
% TF decomposition — coherence is then computed on the warped signals.
if ~isempty(g.timeStretchMarks)
    timestretchopt = {'timestretch', ...
        {g.timeStretchMarks', g.timeStretchRefs}};
    fprintf('\nProcessing first input (EEG) with timewarping\n');
    X = reshape(X, g.frame, g.trials);
    [alltfX, freqs, timesout] = timefreq(X, g.srate, ...
        spectraloptions{:}, timestretchopt{:});

    fprintf('\nProcessing second input (EMG) with timewarping\n');
    Y = reshape(Y, g.frame, g.trials);
    [alltfY] = timefreq(Y, g.srate, ...
        spectraloptions{:}, timestretchopt{:});
else
    % Original newcrossf behaviour — no timewarping
    fprintf('\nProcessing first input (EEG)\n');
    X = reshape(X, g.frame, g.trials);
    [alltfX, freqs, timesout] = timefreq(X, g.srate, spectraloptions{:});

    fprintf('\nProcessing second input (EMG)\n');
    Y = reshape(Y, g.frame, g.trials);
    [alltfY] = timefreq(Y, g.srate, spectraloptions{:});
end
%% END MK MODIFICATION (3)

% =========================================================================
% Compute coherence (identical to newcrossf)
% =========================================================================
tmpprod = alltfX .* conj(alltfY);
if nargout > 6 || strcmpi(g.type, 'phasecoher2') || strcmpi(g.type, 'phasecoher')
    coherresout = alltfX .* conj(alltfY);
end

switch g.type
    case 'crossspec'
        coherres = alltfX .* conj(alltfY);

    case 'coher'
        coherres = sum(alltfX .* conj(alltfY), 3) ./ ...
            sqrt(sum(abs(alltfX).^2, 3) .* sum(abs(alltfY).^2, 3));

    case 'amp'
        alltfX   = abs(alltfX);
        alltfY   = abs(alltfY);
        coherres = ampcorr(alltfX, alltfY, freqs, timesout, g);
        g.alpha  = NaN;
        coherresout = [];

    case 'phasecoher2'
        coherres = sum(coherresout, 3) ./ sum(abs(coherresout), 3);

    case 'phasecoher'
        coherres = sum(coherresout ./ abs(coherresout), 3) / g.trials;
end

% =========================================================================
% Baseline (identical to newcrossf)
% =========================================================================
if size(g.baseline, 2) == 2
    baseln = [];
    for index = 1:size(g.baseline, 1)
        tmptime = find(timesout >= g.baseline(index,1) & ...
            timesout <= g.baseline(index,2));
        baseln = union_bc(baseln, tmptime);
    end
    if isempty(baseln)
        error('No point found in baseline');
    end
else
    if ~isempty(find(timesout < g.baseline))
        baseln = find(timesout < g.baseline);
    else
        baseln = 1:length(timesout);
    end
end
if ~isnan(g.alpha) && isempty(baseln)
    fprintf('my_newcrossf(): no windows in baseline (times<%g).\n', g.baseline);
    return
end
mbase = mean(abs(coherres(:, baseln)'));

% =========================================================================
% Bootstrap (identical to newcrossf)
% =========================================================================
if ~isempty(g.rboot)
    Rbootout = g.rboot;
else
    if ~isnan(g.alpha)
        switch g.type
            case 'coher'
                inputdata = {alltfX alltfY};
                formula   = ['sum(arg1 .* conj(arg2), 3) ./ sqrt( ' ...
                    'sum(abs(arg1).^2,3) .* sum(abs(arg2).^2,3) );'];
            case 'amp'
                inputdata = {abs(alltfX) abs(alltfY)};
            case 'phasecoher2'
                inputdata = {alltfX alltfY};
                formula   = ['tmp = arg1 .* conj(arg2);' ...
                    'res = sum(tmp, 3) ./ sum(abs(tmp),3);'];
            case 'phasecoher'
                inputdata = {alltfX./abs(alltfX) alltfY./abs(alltfY)};
                formula   = 'mean(arg1 .* conj(arg2),3);';
            case 'crossspec'
                inputdata = {alltfX./abs(alltfX) alltfY./abs(alltfY)};
                formula   = 'arg1 .* conj(arg2);';
        end

        % baseline for bootstrap
        if size(g.baseboot, 2) == 1
            if g.baseboot == 0,             baselntmp = [];
            elseif ~isnan(g.baseline(1)),   baselntmp = baseln;
            else,                           baselntmp = find(timesout <= 0);
            end
        else
            baselntmp = [];
            for index = 1:size(g.baseboot, 1)
                tmptime   = find(timesout >= g.baseboot(index,1) & ...
                    timesout <= g.baseboot(index,2));
                baselntmp = union_bc(baselntmp, tmptime);
            end
        end

        if strcmpi(g.boottype, 'shuffle') || strcmpi(g.boottype, 'rand')
            Rbootout = bootstat(inputdata, formula, ...
                'boottype',  g.boottype, ...
                'label',     'coherence', ...
                'bootside',  'upper', ...
                'shuffledim', [2 3], ...
                'dimaccu',   2, ...
                'naccu',     g.naccu, ...
                'alpha',     g.alpha, ...
                'basevect',  baselntmp);
        elseif strcmpi(g.boottype, 'randall')
            Rbootout = bootstat(inputdata, formula, ...
                'boottype',  'rand', ...
                'bootside',  'upper', ...
                'shuffledim', 3, ...
                'naccu',     g.naccu, ...
                'alpha',     g.alpha, ...
                'basevect',  baselntmp);
        else
            Rbootout = bootstat(inputdata, formula, ...
                'boottype',  'shuffle', ...
                'bootside',  'upper', ...
                'shuffledim', 3, ...
                'naccu',     g.naccu, ...
                'alpha',     g.alpha, ...
                'basevect',  baselntmp);
        end
    else
        Rbootout = [];
    end
end

% =========================================================================
% Plot (identical to newcrossf)
% =========================================================================
if strcmpi(g.plotamp, 'on') || strcmpi(g.plotphase, 'on')
    if strcmpi(g.plottype, 'image')
        plotall(coherres, Rbootout, timesout, freqs, mbase, g);
    else
        plotallcurves(coherres, Rbootout, timesout, freqs, mbase, g);
    end
end

% =========================================================================
% Format outputs (identical to newcrossf)
% =========================================================================
Rangle = angle(coherres);
R      = abs(coherres);

return;

% =========================================================================
% SUBFUNCTIONS (identical to newcrossf — copied verbatim)
% =========================================================================

% ------------------
% amplitude correlation
% ------------------
function [coherres, lagmap] = ampcorr(alltfX, alltfY, freqs, timesout, g)

coherres = zeros(length(freqs), length(timesout), length(g.amplag));
alpha    = zeros(length(freqs), length(timesout), length(g.amplag));
countlag = 1;
for lag = g.amplag
    fprintf('Computing %d point lag amplitude correlation...\n', lag);
    for i1 = 1:length(freqs)
        for i2 = max(1, 1-lag):min(length(timesout)-lag, length(timesout))
            if ~isnan(g.alpha)
                [tmp1, tmp2] = corrcoef(...
                    squeeze(alltfX(i1,i2,:)), ...
                    squeeze(alltfY(i1,i2+lag,:)));
                coherres(i1,i2,countlag) = tmp1(1,2);
                alpha(i1,i2,countlag)    = tmp2(1,2);
            else
                tmp1 = corrcoef(...
                    squeeze(alltfX(i1,i2,:)), ...
                    squeeze(alltfY(i1,i2+lag,:)));
                coherres(i1,i2,countlag) = tmp1(1,2);
            end
        end
    end
    countlag = countlag + 1;
end

if length(g.amplag) > 1
    [coherres, lagmap] = max(coherres, [], 3);
    dimsize = length(freqs)*length(timesout);
    alpha   = reshape(alpha((lagmap(:)-1)*dimsize + [1:dimsize]'), ...
        length(freqs), length(timesout));
    lagmap  = g.amplag(lagmap);
    coherres = coherres .* exp(j*lagmap/max(abs(g.amplag)));
else
    lagmap = [];
end

if ~isnan(g.alpha)
    tmpind = find(alpha(:) > g.alpha);
    coherres(tmpind) = 0;
end

% ------------------
% plotting functions
% ------------------
function plotall(R, Rboot, times, freqs, mbase, g)

switch lower(g.plotphase)
    case 'on'
        switch lower(g.plotamp)
            case 'on',  ordinate1 = 0.67; ordinate2 = 0.1; height = 0.33; g.plot = 1;
            case 'off', ordinate2 = 0.1;  height = 0.9;  g.plot = 1;
        end
    case 'off'
        ordinate1 = 0.1; height = 0.9;
        switch lower(g.plotamp)
            case 'on',  ordinate1 = 0.1; height = 0.9; g.plot = 1;
            case 'off', g.plot = 0;
        end
end

Rangle = angle(R);
if ~isreal(R)
    R      = abs(R);
    Rraw   = R;
    setylim = 1;
    if ~isnan(g.baseline)
        R = R - repmat(mbase', [1 g.timesout]);
    end
else
    Rraw    = R;
    setylim = 0;
end

if g.plot
    fprintf('\nNow plotting...\n');
    set(gcf, 'DefaultAxesFontSize', g.AXES_FONT)
    try, icadefs; colormap(feval(DEFAULT_COLORMAP, 256)); catch, end
    pos = get(gca, 'position');
    q   = [pos(1) pos(2) 0 0];
    s   = [pos(3) pos(4) pos(3) pos(4)];
    axis('off')
end

switch lower(g.plotamp)
    case 'on'
        RR = R;
        if ~isnan(g.alpha)
            switch dims(Rboot)
                case 3, RR(find(RR > Rboot(:,:,1) & RR < Rboot(:,:,2))) = 0;
                    Rraw(find(RR > Rboot(:,:,1) & RR < Rboot(:,:,2))) = 0;
                case 2, RR(find(RR < Rboot)) = 0;
                    Rraw(find(RR < Rboot)) = 0;
                case 1, RR(find(RR < repmat(Rboot(:), [1 size(RR,2)]))) = 0;
                    Rraw(find(RR < repmat(Rboot(:), [1 size(Rraw,2)]))) = 0;
            end
        end

        h(6) = axes('Units','Normalized','Position', ...
            [.1 ordinate1 .8 height].*s+q);
        map = hsv(300);
        map = flipud([map(251:end,:); map(1:250,:)]);
        map(151,:) = map(151,:)*0.9;
        colormap(map);

        if ~strcmpi(g.freqscale, 'log')
            try, imagesc(times, freqs, RR, max(max(RR))*[-1 1]);
            catch, imagesc(times, freqs, RR, [-1 1]); end
        else
            try, imagesclogy(times, freqs, RR, max(max(RR))*[-1 1]);
            catch, imagesclogy(times, freqs, RR, [-1 1]); end
        end
        set(gca, 'ydir', 'norm');
        if ~isempty(g.maxamp), caxis([-g.maxamp g.maxamp]); end
        tmpscale = caxis;
        hold on
        plot([0 0], [0 freqs(end)], '--m', 'LineWidth', g.linewidth)
        for i = 1:length(g.marktimes)
            plot([g.marktimes(i) g.marktimes(i)], [0 freqs(end)], ...
                '--m', 'LineWidth', g.linewidth);
        end
        hold off
        set(h(6), 'YTickLabel', [], 'YTick', [])
        set(h(6), 'XTickLabel', [], 'XTick', [])

        h(8) = axes('Position', [.95 ordinate1 .05 height].*s+q);
        if setylim
            cbar(h(8), 151:300, [0 tmpscale(2)]);
        else
            cbar(h(8), 1:300, [-tmpscale(2) tmpscale(2)]);
        end

        h(10) = axes('Units','Normalized', ...
            'Position', [.1 ordinate1-0.1 .8 .1].*s+q);
        Emax = max(R);
        Emin = min(R);
        plot(times, Emin, times, Emax, 'LineWidth', g.linewidth); hold on;
        plot([times(1) times(end)], [0 0], 'LineWidth', 0.7);
        plot([0 0], [-500 500], '--m', 'LineWidth', g.linewidth);
        for i = 1:length(g.marktimes)
            plot([g.marktimes(i) g.marktimes(i)], [-500 500], ...
                '--m', 'LineWidth', g.linewidth);
        end
        if ~isnan(g.alpha) && dims(Rboot) > 1
            switch dims(Rboot)
                case 2, plot(times, mean(Rboot(:,:),1), 'g', 'LineWidth', g.linewidth);
                    plot(times, mean(Rboot(:,:),1), 'k:', 'LineWidth', g.linewidth);
                case 3, plot(times, mean(Rboot(:,:,1),1), 'g', 'LineWidth', g.linewidth);
                    plot(times, mean(Rboot(:,:,1),1), 'k:', 'LineWidth', g.linewidth);
                    plot(times, mean(Rboot(:,:,2),1), 'g', 'LineWidth', g.linewidth);
                    plot(times, mean(Rboot(:,:,2),1), 'k:', 'LineWidth', g.linewidth);
            end
            axis([min(times) max(times) 0 max([Emax(:)' Rboot(:)'])*1.2])
        else
            axis([min(times) max(times) 0 max(Emax)*1.2])
        end
        tick = get(h(10), 'YTick');
        set(h(10), 'YTick', [tick(1); tick(end)])
        set(h(10), 'YAxisLocation', 'right')
        xlabel('Time (ms)')
        ylabel('coh.')

        h(11) = axes('Units','Normalized', ...
            'Position', [0 ordinate1 .1 height].*s+q);
        E = abs(mbase);
        if ~strcmpi(g.freqscale, 'log')
            plot(freqs, E, 'b', 'LineWidth', g.linewidth);
        else
            semilogx(freqs, E, 'b', 'LineWidth', g.linewidth);
        end
        if ~isnan(g.alpha)
            hold on
            if ~strcmpi(g.freqscale, 'log')
                switch dims(Rboot)
                    case 1, plot(freqs, Rboot(:), 'g', 'LineWidth', g.linewidth);
                        plot(freqs, Rboot(:), 'k:', 'LineWidth', g.linewidth);
                    case 2, plot(freqs, mean(Rboot(:,:),2), 'g', 'LineWidth', g.linewidth);
                        plot(freqs, mean(Rboot(:,:),2), 'k:', 'LineWidth', g.linewidth);
                end
            end
            if ~isnan(max(E))
                axis([freqs(1) freqs(end) 0 max([E Rboot(:)'])*1.2]);
            end
        else
            if ~isnan(max(E))
                axis([freqs(1) freqs(end) 0 max(E)*1.2]);
            end
        end
        set(gca, 'xdir', 'rev');
        tick = get(h(11), 'YTick');
        set(h(11), 'YTick', [tick(1); tick(end)]);
        set(h(11), 'View', [90 90])
        xlabel('Freq. (Hz)')
        ylabel('coh.')
end

switch lower(g.plotphase)
    case 'on'
        h(13) = axes('Units','Normalized', ...
            'Position', [.1 ordinate2 .8 height].*s+q);
        if setylim
            if strcmp(g.angleunit, 'ms')
                Rangle = (Rangle/(2*pi)) .* repmat(1000./freqs(:)', 1, length(times));
                maxangle = max(max(abs(Rangle)));
            elseif strcmpi(g.angleunit, 'deg')
                Rangle   = Rangle*180/pi;
                maxangle = 180;
            else
                maxangle = pi;
            end
            Rangle(find(Rraw==0)) = 0;
            if ~strcmpi(g.freqscale, 'log')
                imagesc(times, freqs, Rangle, [-maxangle maxangle]);
            else
                imagesclogy(times, freqs, Rangle, [-maxangle maxangle]);
            end
            hold on
            plot([0 0], [0 freqs(end)], '--m', 'LineWidth', g.linewidth);
            for i = 1:length(g.marktimes)
                plot([g.marktimes(i) g.marktimes(i)], [0 freqs(end)], ...
                    '--m', 'LineWidth', g.linewidth);
            end
            set(gca, 'ydir', 'norm');
            ylabel('Freq. (Hz)')
            xlabel('Time (ms)')
            h(14) = axes('Position', [.95 ordinate2 .05 height].*s+q);
            cbar(h(14), 0, [-maxangle maxangle]);
        else
            axis off;
            text(0, 0.5, 'Real values, no angles');
        end
end

if g.plot
    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
    if length(g.title) > 0
        axes('Position', pos, 'Visible', 'Off');
        h(13) = text(-.05, 1.01, g.title);
        set(h(13), 'VerticalAlignment', 'bottom')
        set(h(13), 'HorizontalAlignment', 'left')
        set(h(13), 'FontSize', g.TITLE_FONT)
    end
    if ~isempty(g.topovec) && strcmpi(g.plotamp,'on') && strcmpi(g.plotphase,'on')
        h(15) = subplot('Position', [-.1 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(1), g.elocs, 'electrodes', 'off', ...
                'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
        else
            topoplot(g.topovec(1,:), g.elocs, 'electrodes', 'off', 'chaninfo', g.chaninfo);
        end
        axis('square')
        h(16) = subplot('Position', [.9 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(2), g.elocs, 'electrodes', 'off', ...
                'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
        else
            topoplot(g.topovec(2,:), g.elocs, 'electrodes', 'off', 'chaninfo', g.chaninfo);
        end
        axis('square')
    end
    try, axcopy(gcf); catch, end
end

% ---------------
function plotallcurves(R, Rboot, times, freqs, mbase, g)
Rangle = angle(R);
pos    = get(gca, 'position');
q      = [pos(1) pos(2) 0 0];
s      = [pos(3) pos(4) pos(3) pos(4)];
if ~isreal(R)
    R    = abs(R);
    Rraw = R;
    if ~isnan(g.baseline)
        R = R - repmat(mbase', [1 g.timesout]);
    end
else
    Rraw    = R;
    setylim = 0;
end
if times(end) > 10000, times = times/1000; timeunit = 's';
else, timeunit = 'ms'; end

alllegend = {};
if strcmpi(g.plotmean, 'on') && freqs(1) ~= freqs(end)
    alllegend = {[num2str(freqs(1)) '-' num2str(freqs(end)) 'Hz']};
else
    for index = 1:length(freqs)
        alllegend{index} = [num2str(freqs(index)) 'Hz'];
    end
end

fprintf('\nNow plotting...\n');
if strcmpi(g.plotamp, 'on')
    if strcmpi(g.plotphase, 'on'), subplot(2,1,1); end
    if isempty(g.maxamp), g.maxamp = 0; end
    plotcurve(times, R, 'maskarray', Rboot, 'title', 'Coherence amplitude', ...
        'xlabel', ['Time (' timeunit ')'], 'ylabel', '0-1', 'ylim', g.maxamp, ...
        'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
        'linewidth', g.linewidth, 'highlightmode', g.highlightmode, ...
        'plotmean', g.plotmean);
end
if strcmpi(g.plotphase, 'on')
    if strcmpi(g.plotamp, 'on'), subplot(2,1,2); end
    plotcurve(times, Rangle/pi*180, 'maskarray', Rboot, 'val2mask', R, ...
        'title', 'Coherence phase', 'xlabel', ['Time (' timeunit ')'], ...
        'ylabel', 'Angle (deg.)', 'ylim', [-180 180], ...
        'vert', g.vert, 'marktimes', g.marktimes, 'legend', alllegend, ...
        'linewidth', g.linewidth, 'highlightmode', g.highlightmode, ...
        'plotmean', g.plotmean);
end
if strcmpi(g.plotamp,'on') || strcmpi(g.plotphase,'on')
    try, icadefs; set(gcf,'color',BACKCOLOR); catch, end
    if length(g.title) > 0
        h(13) = textsc(g.title, 'title');
    end
    if ~isempty(g.topovec)
        h(15) = subplot('Position', [-.1 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(1), g.elocs, 'electrodes', 'off', ...
                'style', 'blank', 'emarkersize1chan', 10);
        else
            topoplot(g.topovec(1,:), g.elocs, 'electrodes', 'off');
        end
        axis('square')
        h(16) = subplot('Position', [.9 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(2), g.elocs, 'electrodes', 'off', ...
                'style', 'blank', 'emarkersize1chan', 10);
        else
            topoplot(g.topovec(2,:), g.elocs, 'electrodes', 'off');
        end
        axis('square')
    end
    try, axcopy(gcf); catch, end
end

% ---------------
function res = dims(array)
res = min(ndims(array), max(size(array,2), size(array,3)));