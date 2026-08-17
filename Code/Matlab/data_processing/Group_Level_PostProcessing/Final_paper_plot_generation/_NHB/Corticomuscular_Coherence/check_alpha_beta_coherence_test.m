%% =========================================================
%  EEG/IC PSD + EMG PSD + EEG/EMG Coherence
%  Inputs:
%     x  -> EEG or IC epochs      [nSamples x nTrials]
%     y  -> EMG epochs (aligned)  [nSamples x nTrials]
%     fs -> sampling rate (Hz), e.g. 500
%
%  This script:
%     1) computes trialwise EEG/IC PSD
%     2) computes trialwise EMG PSD
%     3) computes trialwise coherence
%     4) plots the averages
%     5) reports alpha and beta band summaries
%% =========================================================

%% ---------------- USER INPUTS ----------------
% Example:
% x  = IC_epochs;        % [samples x trials]
% y  = EMG_interpEpochs; % [samples x trials]
% fs = 500;

ic = 2; muscle = 1;
ch = 24; % Cz
time_indx = 1:2000;
x = squeeze(EEG_epoched_main.icaact(ic, time_indx, indx_to_keep));
% x = squeeze(EEG_epoched_main.data(ch, time_indx, indx_to_keep));
y = squeeze(EMG_epoched_data_intrp(muscle, time_indx, indx_to_keep));
fs = 500;

alphaBand = [8 14];
betaBand  = [15 35];

% Welch parameters
win_sec      = 0.25;              % 500 ms
overlap_frac = 0.75;              % 50% overlap

%% ---------------- CHECK INPUTS ----------------
if size(x,1) ~= size(y,1) || size(x,2) ~= size(y,2)
    error('x and y must have the same size: [samples x trials].');
end

if any(isnan(x(:)))
    error('x contains NaNs.');
end

if any(isnan(y(:)))
    error('y contains NaNs.');
end

[nSamples, nTrials] = size(x);

fprintf('Number of samples per trial: %d\n', nSamples);
fprintf('Number of trials: %d\n', nTrials);
fprintf('Sampling rate: %.2f Hz\n', fs);

%% ---------------- WELCH SETTINGS ----------------
win      = round(win_sec * fs);
noverlap = round(overlap_frac * win);
nfft     = max(512, 2^nextpow2(win));

fprintf('Welch window: %d samples (%.3f s)\n', win, win/fs);
fprintf('Overlap: %d samples\n', noverlap);
fprintf('NFFT: %d\n', nfft);

%% ---------------- PREALLOCATE ----------------
Pxx_all = [];
Pyy_all = [];
Cxy_all = [];
Sxy_all = [];

%% ---------------- LOOP OVER TRIALS ----------------
for tr = 1:nTrials
    
    sig_x = double(x(:,tr));
    sig_y = double(y(:,tr));
    
    % Remove DC offset per trial
    sig_x = sig_x - mean(sig_x, 'omitnan');
    sig_y = sig_y - mean(sig_y, 'omitnan');
    
    % Skip bad trials if needed
    if any(~isfinite(sig_x)) || any(~isfinite(sig_y))
        warning('Trial %d contains non-finite values. Skipping.', tr);
        continue
    end
    
    % EEG/IC PSD
    [Pxx, F] = pwelch(sig_x, hamming(win), noverlap, nfft, fs);
    
    % EMG PSD
    [Pyy, ~] = pwelch(sig_y, hamming(win), noverlap, nfft, fs);
    
    % Coherence
    [Cxy, ~] = mscohere(sig_x, sig_y, hamming(win), noverlap, nfft, fs);
    
    % Cross-spectrum
    [Sxy, ~] = cpsd(sig_x, sig_y, hamming(win), noverlap, nfft, fs);
    
    Pxx_all(:,end+1) = Pxx;
    Pyy_all(:,end+1) = Pyy;
    Cxy_all(:,end+1) = Cxy;
    Sxy_all(:,end+1) = Sxy;
end

%% ---------------- CHECK VALID TRIALS ----------------
nValid = size(Pxx_all,2);

if nValid == 0
    error('No valid trials remained after checks.');
end

fprintf('Valid trials used: %d / %d\n', nValid, nTrials);

%% ---------------- AVERAGES ----------------
meanPxx = mean(Pxx_all, 2);
meanPyy = mean(Pyy_all, 2);
meanCxy = mean(Cxy_all, 2);
meanSxy = mean(Sxy_all, 2);

stdPxx = std(Pxx_all, 0, 2);
stdPyy = std(Pyy_all, 0, 2);
stdCxy = std(Cxy_all, 0, 2);

%% ---------------- BAND INDICES ----------------
alphaIdx = F >= alphaBand(1) & F <= alphaBand(2);
betaIdx  = F >= betaBand(1)  & F <= betaBand(2);

%% ---------------- BAND SUMMARIES ----------------
% Mean spectra/coherence across band
alpha_Pxx = mean(meanPxx(alphaIdx));
beta_Pxx  = mean(meanPxx(betaIdx));

alpha_Pyy = mean(meanPyy(alphaIdx));
beta_Pyy  = mean(meanPyy(betaIdx));

alpha_Cxy = mean(meanCxy(alphaIdx));
beta_Cxy  = mean(meanCxy(betaIdx));

alpha_Sxy = mean(abs(meanSxy(alphaIdx)));
beta_Sxy  = mean(abs(meanSxy(betaIdx)));

fprintf('\n===== BAND SUMMARY =====\n');
fprintf('EEG/IC PSD alpha (8-14 Hz): %.4e\n', alpha_Pxx);
fprintf('EEG/IC PSD beta  (15-35 Hz): %.4e\n', beta_Pxx);

fprintf('EMG PSD alpha (8-14 Hz): %.4e\n', alpha_Pyy);
fprintf('EMG PSD beta  (15-35 Hz): %.4e\n', beta_Pyy);

fprintf('Coherence alpha (8-14 Hz): %.4f\n', alpha_Cxy);
fprintf('Coherence beta  (15-35 Hz): %.4f\n', beta_Cxy);

fprintf('|Cross-spectrum| alpha (8-14 Hz): %.4e\n', alpha_Sxy);
fprintf('|Cross-spectrum| beta  (15-35 Hz): %.4e\n', beta_Sxy);

%% ---------------- TRIALWISE BAND VALUES ----------------
alpha_Cxy_trials = mean(Cxy_all(alphaIdx,:), 1);
beta_Cxy_trials  = mean(Cxy_all(betaIdx,:), 1);

alpha_Pyy_trials = mean(Pyy_all(alphaIdx,:), 1);
beta_Pyy_trials  = mean(Pyy_all(betaIdx,:), 1);

alpha_Pxx_trials = mean(Pxx_all(alphaIdx,:), 1);
beta_Pxx_trials  = mean(Pxx_all(betaIdx,:), 1);

%% ---------------- PLOTS ----------------

% 1) EEG/IC PSD
figure;
plot(F, 10*log10(meanPxx), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Average EEG/IC power spectrum');
xlim([0 60]);
grid on;

% 2) EMG PSD
figure;
plot(F, 10*log10(meanPyy), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Average EMG power spectrum');
xlim([0 60]);
grid on;

% 3) Coherence
figure;
plot(F, meanCxy, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Average EEG/IC - EMG coherence');
xlim([0 60]);
ylim([0 1]);
grid on;

% 4) Cross-spectrum magnitude
figure;
plot(F, abs(meanSxy), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|Cross-spectrum|');
title('Average EEG/IC - EMG cross-spectrum magnitude');
xlim([0 60]);
grid on;

% 5) Trialwise coherence alpha vs beta
figure;
boxplot([alpha_Cxy_trials(:), beta_Cxy_trials(:)], 'Labels', {'Alpha','Beta'});
ylabel('Mean coherence');
title('Trialwise coherence: Alpha vs Beta');
grid on;

% 6) Trialwise EMG PSD alpha vs beta
figure;
boxplot([alpha_Pyy_trials(:), beta_Pyy_trials(:)], 'Labels', {'Alpha','Beta'});
ylabel('Mean EMG PSD');
title('Trialwise EMG PSD: Alpha vs Beta');
grid on;

% 7) Trialwise EEG/IC PSD alpha vs beta
figure;
boxplot([alpha_Pxx_trials(:), beta_Pxx_trials(:)], 'Labels', {'Alpha','Beta'});
ylabel('Mean EEG/IC PSD');
title('Trialwise EEG/IC PSD: Alpha vs Beta');
grid on;

%% ---------------- OPTIONAL: COHERENCE PEAKS ----------------
% Peak coherence in 0-60 Hz
freqMask = F >= 0 & F <= 60;
[peakCoh, peakIdx] = max(meanCxy(freqMask));
Fmask = F(freqMask);
peakFreq = Fmask(peakIdx);

fprintf('\nPeak coherence in 0-60 Hz: %.4f at %.2f Hz\n', peakCoh, peakFreq);

%% ---------------- OPTIONAL: RETURN RESULTS STRUCT ----------------
results = struct();
results.F = F;

results.meanPxx = meanPxx;
results.meanPyy = meanPyy;
results.meanCxy = meanCxy;
results.meanSxy = meanSxy;

results.Pxx_all = Pxx_all;
results.Pyy_all = Pyy_all;
results.Cxy_all = Cxy_all;
results.Sxy_all = Sxy_all;

results.alphaBand = alphaBand;
results.betaBand  = betaBand;

results.alpha_Pxx = alpha_Pxx;
results.beta_Pxx  = beta_Pxx;

results.alpha_Pyy = alpha_Pyy;
results.beta_Pyy  = beta_Pyy;

results.alpha_Cxy = alpha_Cxy;
results.beta_Cxy  = beta_Cxy;

results.alpha_Sxy = alpha_Sxy;
results.beta_Sxy  = beta_Sxy;

results.alpha_Cxy_trials = alpha_Cxy_trials;
results.beta_Cxy_trials  = beta_Cxy_trials;

results.alpha_Pyy_trials = alpha_Pyy_trials;
results.beta_Pyy_trials  = beta_Pyy_trials;

results.alpha_Pxx_trials = alpha_Pxx_trials;
results.beta_Pxx_trials  = beta_Pxx_trials;

results.peakCoh  = peakCoh;
results.peakFreq = peakFreq;

fprintf('\nDone. Results are stored in the struct: results\n');