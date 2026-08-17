figure()
Z = 120;
set(gcf, 'Position', [100, 100, 8.5*Z, 5*Z])
trial = 143;
N_epochs = length(Epochs_FlextoFlex_based{trial}.EXP_stream.Ref_angle);
ref = [];
angle = [];
time = [];
for i = 1:N_epochs
    time = cat(2, time, Epochs_FlextoFlex_based{trial}.EXP_stream.Times{i});
    ref = cat(2, ref, Epochs_FlextoFlex_based{trial}.EXP_stream.Ref_angle{i});
    angle = cat(2, angle, Epochs_FlextoFlex_based{trial}.EXP_stream.Encoder_angle{i});
end

% t : time vector
% x : measured signal (column or row vector)
x = ref;
x = x(:);                      % make sure column
t = time;
dt = mean(diff(t));           % sampling interval
fs = 1/dt;                    % sampling frequency
N  = length(x);

% 1) Remove DC and take FFT
x0 = x - mean(x);
X  = fft(x0);

% 2) Find dominant frequency (fundamental)
f  = (0:N-1)*(fs/N);
half = 1:floor(N/2);
[~,k] = max(abs(X(half)));    % index of largest peak
f0 = f(half(k));              % fundamental frequency

% 3) Fit sine at that frequency using linear regression
w0 = 2*pi*f0;
M  = [sin(w0*t(:)) cos(w0*t(:)) ones(N,1)];  % [sin cos offset]
theta = M\x;                  % least-squares solution

A_sin  = theta(1);
A_cos  = theta(2);
offset = theta(3);

A     = hypot(A_sin, A_cos);  % amplitude
phi   = atan2(A_cos, A_sin);  % phase

% 4) Reconstruct clean sinusoid
t_new = linspace(t(1), t(end), length(t));
x_clean = offset + A*sin(w0*t_new + phi);

ref_filt = x_clean;

color_ref_filt = [0, 0, 0, 0.5];
plot(t_new, ref_filt, 'LineWidth', 5, 'Color', color_ref_filt)
% plot3(t_new, zeros(size(t_new)), ref_filt, 'LineWidth', 5, 'Color', color_ref_filt)

hold on

v = 135;
color_angle = [0, 0, 1, 0.5]
plot(t_new(1:v), angle(1:v), 'LineWidth', 5, 'LineStyle', '-', ...
    'Color', color_angle)
% plot3(t_new(1:v), zeros(size(t_new(1:v))), angle(1:v), 'LineWidth', 5, 'LineStyle', '-', ...
%     'Color', color_angle)

m = max(ref_filt);
n = min(ref_filt);
ylim([n-15 m+15])
xlim([t(1) t(end-15)])



x0 = t_new(v);                   
y0 = angle(v);                   

% %%%% Upper arrow showing the knee extension
% hTot = 25;                 % total "arrow" height
% yMid = y0 + hTot/2;        % where arrow starts
% yTop = y0 + hTot;          % arrow tip
% colArrow = [0.30 0.70 0.30];%[0.4660 0.6740 0.1880];
% 
% % 1) lower shaft (line only, under the marker)
% plot([x0 x0],[y0 yMid+5], ...
%      'Color',colArrow,'LineWidth', 5);
% % 2) upper part as pretty annotation arrow
% ax  = gca;
% pos = ax.Position;
% 
% xf = @(x) (x-ax.XLim(1))/diff(ax.XLim)*pos(3)+pos(1);
% yf = @(y) (y-ax.YLim(1))/diff(ax.YLim)*pos(4)+pos(2);
% 
% xN  = xf(x0);
% yN0 = yf(yMid);
% yN1 = yf(yTop);
% 
% annotation('arrow',[xN xN],[yN0 yN1], ...
%     'Color',colArrow, ...
%     'LineWidth',5, ...
%     'HeadStyle','cback2', ...
%     'HeadLength',16,'HeadWidth',14);


% %%%% lower arrow showing the flexion
% hTot = -25;                 % total "arrow" height
% yMid = y0 + hTot/2;        % where arrow starts
% yTop = y0 + hTot;          % arrow tip
% colArrow = [0.85 0.33 0.10];
% 
% % 1) lower shaft (line only, under the marker)
% plot([x0 x0],[y0 yMid-5], ...
%      'Color',colArrow,'LineWidth', 5);
% % 2) upper part as pretty annotation arrow
% ax  = gca;
% pos = ax.Position;
% 
% xf = @(x) (x-ax.XLim(1))/diff(ax.XLim)*pos(3)+pos(1);
% yf = @(y) (y-ax.YLim(1))/diff(ax.YLim)*pos(4)+pos(2);
% 
% xN  = xf(x0);
% yN0 = yf(yMid);
% yN1 = yf(yTop);
% 
% annotation('arrow',[xN xN],[yN0 yN1], ...
%     'Color',colArrow, ...
%     'LineWidth',5, ...
%     'HeadStyle','cback2', ...
%     'HeadLength',16,'HeadWidth',14);



% glowColor = [0 0.4470 0.7410];   % bluish
glowColor = [0, 0.1, 0.8];   % bluish


% --- draw the glow (big, very transparent circles behind) ---
for k = 6:-1:1
    scatter(x0, y0, 200 + 80*k, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', glowColor, ...
        'MarkerFaceAlpha', 0.03*k);  % more alpha near center
    % scatter3(x0, 0, y0, 200 + 80*k, ...
    %     'Marker', 'o', ...
    %     'MarkerEdgeColor', 'none', ...
    %     'MarkerFaceColor', glowColor, ...
    %     'MarkerFaceAlpha', 0.03*k);  % more alpha near center
end

% --- solid center point ---
scatter(x0, y0, 120, ...
    'Marker', 'o', ...
    'MarkerFaceColor', glowColor, ...
    'MarkerEdgeColor', 'w', ...   % white edge if you like
    'LineWidth', 1.5);
% scatter3(x0, 0, y0, 120, ...
%     'Marker', 'o', ...
%     'MarkerFaceColor', glowColor, ...
%     'MarkerEdgeColor', 'w', ...   % white edge if you like
%     'LineWidth', 1.5);




set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
    'XTickLabel', [], 'YTickLabel', [], ...
    'TickLength', [0 0])


ax = gca;
ax.LineWidth = 0.1;          % keep thin grid lines
ax.Box = 'off';              % we'll draw our own box
ax.XAxis.Visible = 'off';   % hide x-axis (ticks + line + labels)
ax.YAxis.Visible = 'off';   % hide y-axis

ax.Units = 'normalized';
pos = ax.Position;           % [x y w h] in figure normalized units

% annotation('rectangle', pos, ...
%     'LineWidth', 5, ...      % bold outer box
%     'Color', 'k');



      