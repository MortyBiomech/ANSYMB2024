function R = decode_demand_within_subject(S)
% DECODE_DEMAND_WITHIN_SUBJECT
% Within-subject decoding of physical-demand level (1 / 3 / 6 bar) from
% sensorimotor band power.
%   - leave-one-session-out cross-validation (4 folds)
%   - regularised (diagonal-covariance) LDA
%   - TRIAL-level features (one row per trial; cycles averaged within trial)
%   - permutation-derived chance, at the subject AND group level
%
% INPUT
%   S : struct array, one element per subject, fields:
%        .X       [nTrials x nFeat]  per-trial band-power features
%        .y       [nTrials x 1]      demand label, coded 1,2,3 for L/M/H
%        .session [nTrials x 1]      session id (e.g. 1..4) -> CV folds
%       (optional) .featNames        1 x nFeat cellstr, for reporting
%
% OUTPUT
%   R : struct with per-subject (obsAcc, subP) and group (groupMean,
%       groupP, nSigSubjects, confusion) results.
%
% NOTES
%   * The unit of analysis is the TRIAL. Build one feature row per trial by
%     averaging the 7-9 within-trial cycles (see assemble template below).
%     Do NOT pass individual cycles as separate rows -- with leave-one-
%     session-out that is safe, but with any split that mixes a trial's
%     cycles across train/test it is pseudoreplication and inflates accuracy.
%   * Balanced accuracy is used so chance = 1/3 regardless of minor imbalance.
%   * Requires Statistics and Machine Learning Toolbox (fitcdiscr).
%
% ---- HOW TO BUILD S (template; adapt to your storage) -------------------
%   bands   = [8 14; 15 30];               % alpha; beta (set your beta range)
%   chanIdx = your_sensorimotor_channel_indices;   % e.g. C3 C1 Cz C2 C4 (+CP/FC)
%   for s = 1:nSub
%       for t = 1:nTrials(s)
%           % trialEEG: [nChan x nTime x nCycle] for trial t (one pressure)
%           S(s).X(t,:)    = trial_bandpower(trialEEG, fs, chanIdx, bands);
%           S(s).y(t,1)    = demandLabel(t);     % 1/2/3
%           S(s).session(t,1) = sessionId(t);    % 1..4
%       end
%   end
%   % If you already have per-cycle power from the ERSP pipeline, just average
%   % it over the trial's cycles instead of calling trial_bandpower.
% -------------------------------------------------------------------------

rng(42,'twister');                  % reproducible permutations
nPerm   = 1000;                     % permutations per subject
classes = [1 2 3];                  % demand levels L / M / H
nC      = numel(classes);
nSub    = numel(S);

obsAcc  = nan(nSub,1);              % observed balanced accuracy
permAcc = nan(nSub,nPerm);         % null balanced accuracy
CM      = zeros(nC,nC,nSub);       % per-subject confusion (rows = true)

for s = 1:nSub
    X = S(s).X;  y = S(s).y(:);  sess = S(s).session(:);

    % observed
    [obsAcc(s), CM(:,:,s)] = loso_balacc(X, y, sess, classes);

    % permutation null: shuffle labels across ALL trials, rerun the SAME CV
    for p = 1:nPerm
        yp = y(randperm(numel(y)));
        permAcc(s,p) = loso_balacc(X, yp, sess, classes);
    end
    fprintf('  subject %2d/%d  acc = %.3f\n', s, nSub, obsAcc(s));
end

% per-subject p-values (+1 correction)
subP = (sum(permAcc >= obsAcc, 2) + 1) ./ (nPerm + 1);

% group-level permutation test on the mean balanced accuracy
obsGroup  = mean(obsAcc);
nullGroup = mean(permAcc, 1);                       % 1 x nPerm null of group mean
groupP    = (sum(nullGroup >= obsGroup) + 1) / (nPerm + 1);

% subject-averaged, row-normalised confusion (errors adjacency check)
CMnorm = CM ./ sum(CM, 2);
CMmean = mean(CMnorm, 3, 'omitnan');

R = struct('obsAcc',obsAcc, 'subP',subP, 'permAcc',permAcc, ...
           'groupMean',obsGroup, 'groupP',groupP, ...
           'nSigSubjects',sum(subP<0.05), 'chance',1/nC, ...
           'confusion',CMmean);

fprintf('\n=== Within-subject demand decoding (leave-one-session-out) ===\n');
fprintf('Chance (balanced)            = %.3f\n', 1/nC);
fprintf('Group mean balanced accuracy = %.3f  (permutation p = %.4f)\n', obsGroup, groupP);
fprintf('Subjects > own perm threshold (p<.05): %d / %d\n', sum(subP<0.05), nSub);
fprintf('Subject-averaged confusion (rows = true L/M/H, cols = predicted):\n');
disp(CMmean);
end


% ========================================================================
function [bacc, cm] = loso_balacc(X, y, sess, classes)
% Leave-one-session-out CV -> balanced accuracy + confusion matrix.
us   = unique(sess);
yhat = nan(size(y));
for f = 1:numel(us)
    te = sess == us(f);
    tr = ~te;

    % z-score on TRAIN ONLY, apply to test (prevents leakage)
    mu = mean(X(tr,:), 1);
    sd = std(X(tr,:), 0, 1);  sd(sd == 0) = 1;
    Xtr = (X(tr,:) - mu) ./ sd;
    Xte = (X(te,:) - mu) ./ sd;

    mdl = fitcdiscr(Xtr, y(tr), 'DiscrimType','diagLinear');  % regularised LDA
    yhat(te) = predict(mdl, Xte);
end
cm     = confusionmat(y, yhat, 'Order', classes);   % rows = true class
recall = diag(cm) ./ sum(cm, 2);                    % per-class recall
bacc   = mean(recall, 'omitnan');                   % balanced accuracy
end


% ========================================================================
function bp = trial_bandpower(trialEEG, fs, chanIdx, bands)
% Mean log band power per channel/band for ONE trial, averaged over cycles.
%   trialEEG : [nChan x nTime x nCycle]
%   bands    : [nBand x 2] = [lo hi] Hz
% Returns a 1 x (numel(chanIdx)*nBand) row vector.
% (Fallback feature extractor; substitute your ERSP-pipeline power if you
%  already have it. pwelch needs the Signal Processing Toolbox.)
nB  = size(bands,1);
nCh = numel(chanIdx);
P   = zeros(nCh, nB);
nCyc = size(trialEEG, 3);
for e = 1:nCyc
    for ci = 1:nCh
        x = double(squeeze(trialEEG(chanIdx(ci), :, e)));
        [pxx, fHz] = pwelch(x - mean(x), [], [], [], fs);
        for b = 1:nB
            m = fHz >= bands(b,1) & fHz < bands(b,2);
            P(ci,b) = P(ci,b) + trapz(fHz(m), pxx(m));
        end
    end
end
P  = P / nCyc;                  % mean over cycles
bp = log(P(:)' + eps);          % log power, row vector
end