function cmap = mortezaMap(n)
% MORTEZAMAP  Perceptually-smooth, high-contrast colormap
%   cmap = mortezaMap(n) returns an n-by-3 colormap (default n = 256)
%   built from your 9 anchor colors. It interpolates in CIELAB and
%   enforces a monotonic lightness ramp for clear value ordering.
%
%   Requires Image Processing Toolbox (rgb2lab/lab2rgb). If you don’t
%   have it, let me know and I’ll give you a toolbox-free version.

if nargin < 1, n = 256; end

% --- 1) Anchor colors in natural hue order (blue → red)
hex = [ ...
 "3849E1"  % Palatinate blue
 "8BC9FF"  % Maya blue
 "9FE1FF"  % Uranian Blue
 "D1FF88"  % Mindaro
 "AFFA37"  % Green Yellow
 "FFE491"  % Jasmine
 "FFB99D"  % Apricot
 "FF6E35"  % Orange (Crayola)
 "A80600"];% Turkey red

% --- 2) Hex → RGB in [0,1]
rgb = reshape(sscanf(join(hex,'')','%2x'),3,[])';
rgb = rgb ./ 255;

% --- 3) sRGB → CIELAB
Lab = rgb2lab(rgb);

% --- 4) Reparameterize along the path for smooth interpolation
%     (arc length in a*b* plane to keep hue transitions smooth)
ab  = Lab(:,2:3);
d   = [0; sqrt(sum(diff(ab,1,1).^2,2))];
s   = cumsum(d);            % cumulative distance
s   = s ./ s(end);          % normalize to [0,1]

% --- 5) Make lightness (L*) strictly monotonic for contrast
L_target = linspace(30, 95, n)';   % tweak endpoints if you like

% --- 6) Interpolate a* and b* with splines; combine with L_target
si = linspace(0,1,n)';
ai = interp1(s, Lab(:,2), si, 'pchip');
bi = interp1(s, Lab(:,3), si, 'pchip');
L  = L_target;

Lab_i = [L, ai, bi];

% --- 7) Back to RGB, clip to gamut
cmap = lab2rgb(Lab_i);
cmap = min(max(cmap,0),1);

end
