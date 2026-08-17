% PLOT_CMC_PAIR - Generate and save all informative CMC plots for one
%   cluster-muscle pair. X-axis is expressed as percentage of movement
%   cycle (0% = first valid TF time point ~ FlxS,
%                100% = last warped event = ExtE).
%
%   Time axis is cropped to [timesout(1), warpingvalues(end)] — i.e.
%   from the first valid TF output point to the last warped event.
%   This excludes pre-movement baseline and post-movement tail.
%
% Usage (standalone from saved struct):
%   >> load('sub-18_CMC.mat');
%   >> s = cmc(1,1);
%   >> plot_CMC_pair(s.coherres, s.crossspec_warped,        ...
%                   s.alltfX, s.alltfX_pow_warped,          ...
%                   s.alltfY, s.alltfY_pow_warped,          ...
%                   s.timesout, s.freqsout,                 ...
%                   s.warpingvalues, s.condition_vector,    ...
%                   s.cluster_name, s.muscle_name,          ...
%                   s.subject_id, params, figSavePath);
%
% Inputs:
%   coherres          - [freq x time] warped MSC coherence
%   crossspec_warped  - [freq x time x trials] warped cross-spectrum
%   alltfX            - [freq x time x trials] EEG TF unwarped
%   alltfX_pow_warped - [freq x time x trials] warped EEG power
%   alltfY            - [freq x time x trials] EMG TF unwarped
%   alltfY_pow_warped - [freq x time x trials] warped EMG power
%   timesout          - [1 x nTimes] ms — full TF time vector
%   freqsout          - [1 x nFreqs] Hz
%   warpingvalues     - [1 x nEvents] group median event latencies ms
%                       warpingvalues(1) = FlxS (should be ~0 or timesout(1))
%                       warpingvalues(end) = ExtE = 100% of cycle
%   condition_vector  - [nTrials x 1] condition label per trial
%   clusterName       - string e.g. 'Left_Prim_Motor'
%   muscleName        - string e.g. 'VL'
%   subject_id        - integer
%   params            - struct with optional fields:
%       .event_labels - cell array e.g. {'FlxS','FlxE/ExtS','ExtE'}
%       .freqs        - [min max] Hz for display limits
%   figSavePath       - folder to save figures
%
% Figures produced:
%   Fig 1 - TF coherence image (% cycle x-axis)
%   Fig 2 - Condition-specific TF coherence
%   Fig 3 - Mean coherence spectrum
%   Fig 4 - Beta-band CMC time course (% cycle x-axis)
%   Fig 5 - Coherence phase image (% cycle x-axis)
%   Fig 6 - EEG and EMG power spectra
%
% Author: Morteza Khosrotabar, TU Darmstadt, 2025

function plot_CMC_pair(coherres, crossspec_warped,                     ...
    alltfX, alltfX_pow_warped,                                          ...
    alltfY, alltfY_pow_warped,                                          ...
    timesout, freqsout,                                                 ...
    warpingvalues, condition_vector,                                    ...
    clusterName, muscleName, subject_id,                               ...
    params, figSavePath)

% =========================================================================
% Setup
% =========================================================================
if ~exist(figSavePath,'dir'), mkdir(figSavePath); end

pairTitle = sprintf('sub-%d | %s — %s', subject_id, ...
    strrep(clusterName,'_',' '), muscleName);
freqLim   = [freqsout(1), freqsout(end)];

% -------------------------------------------------------------------------
% Define time crop window:
%   start: timesout(1)         — first valid TF output (~58ms after epoch start)
%   end:   warpingvalues(end)  — last warped event (ExtE) = 100% of cycle
% -------------------------------------------------------------------------
t_start = timesout(1);                % ms — first valid TF point
t_end   = warpingvalues(end);         % ms — last warped event

timeIdx      = timesout >= t_start & timesout <= t_end;
timesout_c   = timesout(timeIdx);     % cropped ms time vector
coherres_c   = coherres(:, timeIdx);
crossspec_c  = crossspec_warped(:, timeIdx, :);
alltfXp_c    = alltfX_pow_warped(:, timeIdx, :);
alltfYp_c    = alltfY_pow_warped(:, timeIdx, :);

% -------------------------------------------------------------------------
% Convert ms time axis to percent of movement cycle
%   0%   = t_start (first valid TF point ≈ FlxS)
%   100% = t_end   (last warped event = ExtE)
% -------------------------------------------------------------------------
ms_range   = t_end - t_start;                        % total ms span
pct_axis   = (timesout_c - t_start) / ms_range * 100; % [1 x nTimes_c] %

% convert warpingvalues to percent for event lines
evPct = (warpingvalues - t_start) / ms_range * 100;  % [1 x nEvents] %
% clip to [0 100] in case any event is slightly outside
evPct = max(0, min(100, evPct));

pctLim = [0, 100];   % x-axis limits in %

% event labels
if isfield(params,'event_labels') && ...
        length(params.event_labels) == length(warpingvalues)
    evLabels = params.event_labels;
else
    evLabels = arrayfun(@(x) sprintf('E%d',x), ...
        1:length(warpingvalues), 'UniformOutput', false);
end

% conditions
conditions = unique(condition_vector);
nConds     = length(conditions);

% beta band
betaIdx    = freqsout >= 15 & freqsout <= 30;

% coherence color limits
clim_coh = [0, max(coherres_c(:))];

fprintf('    Time crop  : %.1f to %.1f ms\n', t_start, t_end);
fprintf('    Pct axis   : %.1f to %.1f %%\n', pct_axis(1), pct_axis(end));
fprintf('    Event %%   :'); fprintf(' %.1f%%', evPct); fprintf('\n');

% =========================================================================
% Figure 1 — TF coherence image (all trials, % cycle x-axis)
% =========================================================================
fig1 = figure('Name',['Fig1_CMC_TF_' clusterName '_' muscleName], ...
    'Color','w','Units','inches','Position',[0 0 9 4.5]);

contourf(pct_axis, freqsout, coherres_c, 200, 'linecolor','none');
set(gca,'yscale','log','ydir','norm',                                   ...
    'xlim',pctLim,'ylim',freqLim,                                       ...
    'ytick',[4 8 14 30 60],                                             ...
    'clim',clim_coh,                                                    ...
    'FontName','Arial','FontSize',12,'box','on');
colormap(hot);
c = colorbar;
ylabel(c,'Coherence (MSC)','FontSize',11,'FontName','Arial');
xlabel('Cycle (%)','FontSize',13,'FontWeight','bold','FontName','Arial');
ylabel('Frequency (Hz)','FontSize',13,'FontWeight','bold','FontName','Arial');
title(sprintf('TF Coherence — %s', pairTitle), ...
    'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
hold on;
add_event_lines(evPct, evLabels, freqLim, 'k');
hold off;
savethisfig(fig1, figSavePath, 'Fig1_CMC_TF');

% =========================================================================
% Figure 2 — Condition-specific TF coherence (% cycle x-axis)
% =========================================================================
if nConds > 1
    fig2 = figure('Name',['Fig2_CMC_byCond_' clusterName '_' muscleName], ...
        'Color','w','Units','inches','Position',[0 0 4*nConds 4.5]);

    coh_by_cond = zeros(length(freqsout), length(pct_axis), nConds);
    for k = 1:nConds
        condIdx = condition_vector == conditions(k);
        cs_k    = crossspec_c(:,:,condIdx);
        xp_k    = alltfXp_c(:,:,condIdx);
        yp_k    = alltfYp_c(:,:,condIdx);
        num_k   = sum(cs_k, 3);
        den_k   = sqrt(sum(xp_k,3) .* sum(yp_k,3));
        coh_by_cond(:,:,k) = abs(num_k ./ den_k);
    end
    clim_c = [0, max(coh_by_cond(:))];

    for k = 1:nConds
        condIdx = condition_vector == conditions(k);
        subplot(1, nConds, k);
        contourf(pct_axis, freqsout, coh_by_cond(:,:,k), ...
            200, 'linecolor','none');
        set(gca,'yscale','log','ydir','norm',                           ...
            'xlim',pctLim,'ylim',freqLim,                               ...
            'ytick',[4 8 14 30 60],'clim',clim_c,                      ...
            'FontName','Arial','FontSize',11,'box','on');
        colormap(hot);
        if k==nConds
            c2=colorbar;
            ylabel(c2,'Coherence','FontSize',10,'FontName','Arial');
        end
        title(sprintf('Cond %g  (n=%d)', conditions(k), sum(condIdx)), ...
            'FontSize',12,'FontName','Arial','FontWeight','bold');
        if k==1
            ylabel('Frequency (Hz)','FontSize',12,'FontWeight','bold','FontName','Arial');
        end
        xlabel('Cycle (%)','FontSize',12,'FontWeight','bold','FontName','Arial');
        hold on; add_event_lines(evPct, evLabels, freqLim, 'k'); hold off;
    end
    sgtitle(sprintf('Per-condition CMC — %s', pairTitle), ...
        'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
    savethisfig(fig2, figSavePath, 'Fig2_CMC_byCond');
end

% =========================================================================
% Figure 3 — Mean coherence spectrum (frequency profile, averaged over cycle)
% =========================================================================
fig3 = figure('Name',['Fig3_CMC_spectrum_' clusterName '_' muscleName], ...
    'Color','w','Units','inches','Position',[0 0 5.5 4.5]);

meanCohFreq = mean(coherres_c, 2);   % [nFreqs x 1]

semilogx(freqsout, meanCohFreq, 'k-', 'LineWidth',2.5);
set(gca,'xlim',freqLim,'xtick',[4 8 14 30 60],                        ...
    'FontName','Arial','FontSize',12,'box','off');
hold on;
yl = ylim;
fill([15 30 30 15],[yl(1) yl(1) yl(2) yl(2)],                         ...
    [0.85 0.85 1],'EdgeColor','none','FaceAlpha',0.4);
semilogx(freqsout, meanCohFreq, 'k-', 'LineWidth',2.5);
text(17, yl(2)*0.93,'Beta','FontSize',9,'Color',[0.4 0.4 0.8],        ...
    'FontName','Arial','FontWeight','bold');
hold off;
xlabel('Frequency (Hz)','FontSize',13,'FontWeight','bold','FontName','Arial');
ylabel('Mean Coherence','FontSize',13,'FontWeight','bold','FontName','Arial');
title(sprintf('Mean CMC Spectrum — %s', pairTitle), ...
    'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
grid on;
savethisfig(fig3, figSavePath, 'Fig3_CMC_spectrum');

% =========================================================================
% Figure 4 — Beta-band CMC time course (% cycle x-axis)
% =========================================================================
fig4 = figure('Name',['Fig4_CMC_beta_' clusterName '_' muscleName], ...
    'Color','w','Units','inches','Position',[0 0 8 3.5]);

meanBeta = mean(coherres_c(betaIdx,:), 1);   % [1 x nTimes_c]

plot(pct_axis, meanBeta, 'k-', 'LineWidth',2);
hold on;
fill([pct_axis, fliplr(pct_axis)],                                     ...
     [meanBeta, zeros(1,length(meanBeta))],                             ...
     [0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none');
add_event_lines(evPct, evLabels, [0, max(meanBeta)*1.2], 'k');
hold off;
set(gca,'xlim',pctLim,'ylim',[0, max(meanBeta)*1.2],                  ...
    'xtick',[0 25 50 75 100],                                           ...
    'FontName','Arial','FontSize',12,'box','off');
xlabel('Cycle (%)','FontSize',13,'FontWeight','bold','FontName','Arial');
ylabel('Mean Beta Coherence','FontSize',13,'FontWeight','bold','FontName','Arial');
title(sprintf('Beta CMC Time Course — %s', pairTitle), ...
    'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
grid on;
savethisfig(fig4, figSavePath, 'Fig4_CMC_beta');

% =========================================================================
% Figure 5 — Coherence phase image (% cycle x-axis)
%   Mean phase angle of cross-spectrum across trials.
%   Positive = EEG leads EMG; Negative = EMG leads EEG.
% =========================================================================
fig5 = figure('Name',['Fig5_CMC_phase_' clusterName '_' muscleName], ...
    'Color','w','Units','inches','Position',[0 0 9 4.5]);

meanPhase_deg = angle(mean(crossspec_c, 3)) * 180 / pi;  % [freq x time]

contourf(pct_axis, freqsout, meanPhase_deg, 200, 'linecolor','none');
set(gca,'yscale','log','ydir','norm',                                   ...
    'xlim',pctLim,'ylim',freqLim,                                       ...
    'ytick',[4 8 14 30 60],                                             ...
    'clim',[-180 180],                                                  ...
    'FontName','Arial','FontSize',12,'box','on');
colormap(hsv(256));
c3=colorbar;
ylabel(c3,'Phase angle (deg)','FontSize',11,'FontName','Arial');
xlabel('Cycle (%)','FontSize',13,'FontWeight','bold','FontName','Arial');
ylabel('Frequency (Hz)','FontSize',13,'FontWeight','bold','FontName','Arial');
title(sprintf('CMC Phase — %s', pairTitle), ...
    'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
annotation('textbox',[0.01 0.01 0.5 0.06],'String', ...
    'Positive = EEG leads EMG | Negative = EMG leads EEG', ...
    'FontSize',8,'EdgeColor','none','FontName','Arial','Color',[0.4 0.4 0.4]);
hold on; add_event_lines(evPct, evLabels, freqLim, 'w'); hold off;
savethisfig(fig5, figSavePath, 'Fig5_CMC_phase');

% =========================================================================
% Figure 6 — EEG and EMG power spectra (normalised for overlay)
% =========================================================================
fig6 = figure('Name',['Fig6_power_spectra_' clusterName '_' muscleName], ...
    'Color','w','Units','inches','Position',[0 0 6 4.5]);

meanEEGpow   = mean(mean(abs(alltfX).^2, 3), 2);     % [nFreqs x 1]
meanEMGpow   = mean(mean(alltfY_pow_warped, 3), 2);  % [nFreqs x 1]
meanEEGpow_n = meanEEGpow / max(meanEEGpow);
meanEMGpow_n = meanEMGpow / max(meanEMGpow);

semilogx(freqsout, meanEEGpow_n, 'b-', 'LineWidth',2.5); hold on;
semilogx(freqsout, meanEMGpow_n, 'r-', 'LineWidth',2.5);
yl2=ylim;
fill([15 30 30 15],[yl2(1) yl2(1) yl2(2) yl2(2)],                     ...
    [0.85 0.85 1],'EdgeColor','none','FaceAlpha',0.35);
semilogx(freqsout, meanEEGpow_n, 'b-', 'LineWidth',2.5);
semilogx(freqsout, meanEMGpow_n, 'r-', 'LineWidth',2.5);
text(17, yl2(2)*0.93,'Beta','FontSize',9,'Color',[0.4 0.4 0.8],       ...
    'FontName','Arial','FontWeight','bold');
hold off;
set(gca,'xlim',freqLim,'xtick',[4 8 14 30 60],                        ...
    'FontName','Arial','FontSize',12,'box','off');
legend({'EEG IC (norm.)','EMG (norm.)'},'Location','northeast',        ...
    'FontSize',10,'FontName','Arial');
xlabel('Frequency (Hz)','FontSize',13,'FontWeight','bold','FontName','Arial');
ylabel('Normalised Power','FontSize',13,'FontWeight','bold','FontName','Arial');
title(sprintf('EEG & EMG Power — %s', pairTitle), ...
    'FontSize',13,'FontName','Arial','FontWeight','bold','Interpreter','none');
grid on;
savethisfig(fig6, figSavePath, 'Fig6_power_spectra');

fprintf('    Figures saved to:\n    %s\n', figSavePath);

end

% =========================================================================
% Helper: add vertical event lines with labels (x-axis in %)
% =========================================================================
function add_event_lines(evPct, evLabels, yLim, lineColor)
yPos = yLim(1) + (yLim(2)-yLim(1)) * 0.78;  % label height
for L = 1:length(evPct)
    if L==1 || L==length(evPct)
        xline(evPct(L), '-', 'Color',lineColor, 'LineWidth',1.5);
    else
        xline(evPct(L), '--', 'Color',lineColor, 'LineWidth',1.5);
    end
    text(evPct(L)+0.8, yPos, evLabels{L},                             ...
        'FontSize',8,'FontName','Arial','FontWeight','bold',           ...
        'Rotation',90,'Color',lineColor);
end
end

% =========================================================================
% Helper: save figure as png and fig
% =========================================================================
function savethisfig(fig, folder, fname)
if ~exist(folder,'dir'), mkdir(folder); end
saveas(fig, fullfile(folder, [fname '.png']));
saveas(fig, fullfile(folder, [fname '.fig']));
end