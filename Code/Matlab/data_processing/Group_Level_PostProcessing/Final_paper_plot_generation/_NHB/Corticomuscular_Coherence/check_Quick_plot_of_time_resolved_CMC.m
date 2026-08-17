%% Quick plot of time-resolved corticomuscular coherence
%  Run this after my_newcrossf has been called and outputs are in workspace.
%
%  Expected variables in workspace:
%   coh           - [nfreqs x ntimes] coherence magnitude
%   timesout      - [1 x ntimes] time vector (ms)
%   freqsout      - [1 x nfreqs] frequency vector (Hz)
%   warpingvalues - [1 x nEvents] group median event latencies (ms)
%                   e.g. [0, middleEventMs, lastEventMs]
%
%  Author: Morteza Khosrotabar, TU Darmstadt, 2026

%% ---- settings ----
clusterName = 'Left Prim Motor';    % change as needed
muscleName  = 'VL';                 % change as needed
freqLim     = [4 200];               % Hz — covers delta to lower gamma (CMC relevant)
eventLabels = {'FlxS', 'FlxE / ExtS', 'ExtE'};

%% ---- crop time axis to [FlxS, ExtE] i.e. [warpingvalues(1), warpingvalues(end)] ----
% Before FlxS and after ExtE are outside the gait cycle of interest.
% warpingvalues(1) should be 0 (FlxS), warpingvalues(end) is ExtE latency.
warpingvalues = EEG_epoched_main.timewarp.warpto;
warpingvalues(1) = timesout(1);
% timeLim      = [warpingvalues(1), warpingvalues(end)];  % ms
timeLim      = [timesout(1), timesout(end)];  % ms
timeIdx      = timesout >= timeLim(1) & timesout <= timeLim(2);
timesout_crop = timesout(timeIdx);
% coh_crop      = coh(:, timeIdx); 
% coh_crop      = coherres(:, timeIdx); 
coh_crop      = Rangle(:, timeIdx);

%% ---- crop frequency axis ----
freqIdx       = freqsout >= freqLim(1) & freqsout <= freqLim(2);
freqsout_crop = freqsout(freqIdx);
coh_crop      = coh_crop(freqIdx, :);

%% ---- color limits ----
clim = [0, max(coh_crop(:))];
% clim = [min(coh_crop(:)), max(coh_crop(:))];

%% ================================================================
%  Figure 1: Time-frequency coherence image
% ================================================================
figure('Name', ['CMC: ' clusterName ' — ' muscleName], ...
    'Color', 'w', 'Units', 'inches', 'Position', [1 1 8 4]);

contourf(timesout_crop, freqsout_crop, coh_crop, 200, 'linecolor', 'none');
set(gca, ...
    'yscale',   'log', ...
    'ydir',     'norm', ...
    'xlim',     timeLim, ...
    'ylim',     freqLim, ...
    'ytick',    [4 8 14 30 60 150], ...
    'clim',     clim, ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'box',      'on');

colormap("parula");
c = colorbar;
ylabel(c, 'Coherence', 'FontSize', 12, 'FontName', 'Arial');
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
title(['Figure 2: ', clusterName ' — ' muscleName ' CMC'], ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

% event lines
hold on;
for L = 1:length(warpingvalues)
    if L == 1 || L == length(warpingvalues)
        xline(warpingvalues(L), '-', 'LineWidth', 1.5, 'Color', [1, 1, 1]);
    else
        xline(warpingvalues(L), '--', 'LineWidth', 1.5, 'Color', [1, 1, 1]);
    end
    text(warpingvalues(L) + 70, freqLim(2)*0.2, eventLabels{L}, ...
        'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold', ...
        'Rotation', 90, 'Color', [1, 1, 1]);
end
hold off;


%% ================================================================
%  Figure 2: Mean coherence spectrum (averaged across gait cycle)
% ================================================================
figure('Name', ['Mean CMC spectrum: ' clusterName ' — ' muscleName], ...
    'Color', 'w', 'Units', 'inches', 'Position', [1 1 5 4]);

meanCohByFreq = mean(coh_crop, 2);  % [nfreqs x 1]

semilogx(freqsout_crop, meanCohByFreq, 'k-', 'LineWidth', 2);
set(gca, ...
    'xlim',     freqLim, ...
    'xtick',    [4 8 14 30 60], ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'box',      'off');

hold on;
xline(15, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
xline(30, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
text(15.5, max(meanCohByFreq)*0.95, 'Beta', ...
    'FontSize', 9, 'Color', [0.5 0.5 0.5], 'FontName', 'Arial');
hold off;

xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Mean Coherence',  'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
title(['Mean CMC spectrum: ' clusterName ' — ' muscleName], ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
grid on;

%% ================================================================
%  Figure 3: Beta-band CMC time course across the gait cycle
% ================================================================
figure('Name', ['Beta CMC time course: ' clusterName ' — ' muscleName], ...
    'Color', 'w', 'Units', 'inches', 'Position', [1 1 7 3]);

betaIdx      = freqsout_crop >= 15 & freqsout_crop <= 30;
meanCohBeta  = mean(coh_crop(betaIdx, :), 1);  % [1 x ntimes_crop]

plot(timesout_crop, meanCohBeta, 'k-', 'LineWidth', 2);
hold on;
fill([timesout_crop, fliplr(timesout_crop)], ...
     [meanCohBeta,   zeros(1, length(meanCohBeta))], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

for L = 1:length(warpingvalues)
    if L == 1 || L == length(warpingvalues)
        xline(warpingvalues(L), '-k', 'LineWidth', 1.2);
    else
        xline(warpingvalues(L), '--k', 'LineWidth', 1.2);
    end
    text(warpingvalues(L) + 10, max(meanCohBeta)*0.9, eventLabels{L}, ...
        'FontSize', 8, 'FontName', 'Arial', 'FontWeight', 'bold', ...
        'Rotation', 90, 'Color', 'k');
end
hold off;

set(gca, ...
    'xlim',     timeLim, ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'box',      'off');

xlabel('Time (ms)',           'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Mean Beta Coherence', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
title(['Beta CMC time course: ' clusterName ' — ' muscleName], ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
grid on;