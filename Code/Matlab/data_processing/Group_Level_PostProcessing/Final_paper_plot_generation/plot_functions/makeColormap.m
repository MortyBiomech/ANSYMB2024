function cmap = makeColormap(rgb_anchors, n, method, param)
% MAKECOLORMAP  Interpolate a custom colormap with non-linear spacing.
%   rgb_anchors : Kx3 in 0–255 or 0–1 (auto-detected)
%   n           : number of output colors (default 256)
%   method      : 'linear' (default) | 'log' | 'gamma' | 'sigmoid'
%   param       : curve parameter
%                 - 'log'    -> epsilon (default 1e-3)
%                 - 'gamma'  -> gamma (>0, default 2)
%                 - 'sigmoid'-> steepness k (default 10)
%
% Example:
%   rgb = [ 0 132 142; 27 189 201; 35 207 220; 40 217 231; ...
%           70 229 241; 90 235 246; 110 239 249; 125 242 251 ];
%   cmap = makeColormap(rgb, 256, 'gamma', 1.8); colormap(cmap)

    if nargin < 2 || isempty(n), n = 256; end
    if nargin < 3 || isempty(method), method = 'linear'; end

    % Normalize to [0,1] if needed
    if max(rgb_anchors(:)) > 1, rgb = rgb_anchors./255; else, rgb = rgb_anchors; end

    % Anchor parameter (uniform along the anchor list)
    K  = size(rgb,1);
    tA = linspace(0,1,K);

    % Target parameter (nonlinear warping)
    u = linspace(0,1,n);   % base uniform grid
    switch lower(method)
        case 'linear'
            t = u;
        case 'log'
            if nargin < 4 || isempty(param), param = 1e-3; end  % epsilon to avoid log(0)
            epsl = max(param, 1e-9);
            t = (log10(u*(1-epsl) + epsl) - log10(epsl)) / (log10(1) - log10(epsl));
        case 'gamma'
            if nargin < 4 || isempty(param), param = 2; end     % >1 packs colors near start
            t = u.^param;
        case 'sigmoid'
            if nargin < 4 || isempty(param), param = 10; end    % larger = sharper mid focus
            k = param;
            t = (1 ./ (1 + exp(-k*(u-0.5))) - 1/(1+exp(k/2))) / ...
                (1/(1+exp(-k/2)) - 1/(1+exp(k/2)));
        otherwise
            error('Unknown method: %s', method);
    end
    t = min(max(t,0),1);   % clamp

    % Interpolate (pchip avoids flat spots; 'linear' also OK)
    cmap = interp1(tA, rgb, t, 'pchip');

    % Clamp to [0,1]
    cmap = min(max(cmap,0),1);
end
