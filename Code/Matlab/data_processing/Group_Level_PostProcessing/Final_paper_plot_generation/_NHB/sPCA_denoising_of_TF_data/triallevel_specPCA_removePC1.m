function out = triallevel_specPCA_removePC1(tbl, freqs, times, baseline_ms)
    % triallevel_specPCA_removePC1
    % 1) Compute trial-level ERSP (dB) from epoch-level complex TF
    % 2) Perform spectral PCA across (trial x time) observations
    % 3) Remove PC1 and return denoised trial-level ERSP maps
    
    % -------------------------
    % Inputs
    %   tbl.tf         : cell, each cell is [nFreq x nTime] complex TF for one epoch
    %   tbl.trial      : numeric trial id for each epoch
    %   freqs          : [nFreq x 1]
    %   times          : [1 x nTime] (ms)
    %   baseline_ms    : [t1 t2] in ms, e.g., [-500 0]
    % -------------------------
    
    % assert(iscell(tbl.tf), 'tbl.tf must be a cell array (one cell per epoch).');
    % nEpoch = height(tbl);
    % 
    % % --- baseline indices
    % baseIdx = times >= baseline_ms(1) & times <= baseline_ms(2);
    % assert(any(baseIdx), 'Baseline window does not overlap times vector.');
    % 
    % % --- group epochs by trial
    % trialIDs = unique(tbl.trial(:))';
    % nTrial   = numel(trialIDs);
    % 
    % % --- preallocate
    % % trial ERSP will be [nFreq x nTime x nTrial]
    % tf0 = tbl.tf{1};
    % [nFreq, nTime] = size(tf0);
    % ERSP_trial = nan(nFreq, nTime, nTrial);
    % 
    % for iT = 1:nTrial
    %     tr = trialIDs(iT);
    %     idx = (tbl.trial == tr);
    % 
    %     % stack epochs: [nFreq x nTime x nEpochsInTrial]
    %     tfStack = cat(3, tbl.tf{idx});
    % 
    %     % power, then average across epochs in this trial
    %     P = mean(abs(tfStack).^2, 3);                 % [nFreq x nTime]
    % 
    %     % baseline normalize per frequency (ERSP-like, in dB)
    %     B = mean(P(:, baseIdx), 2);                   % [nFreq x 1]
    %     B = max(B, eps);                              % avoid divide-by-zero
    %     ERSP_trial(:, :, iT) = 10*log10(P ./ B);      % [nFreq x nTime]
    % end
    
    % =========================
    % Spectral PCA across observations (trial x time)
    % =========================
    
    % Build X: rows are (time, trial), columns are frequency
    % ERSP_trial: [F x T x R] -> permute to [T x R x F] -> reshape to [(T*R) x F]
    X = reshape(permute(ERSP_trial, [2 3 1]), [], nFreq);
    
    % mean-center per frequency
    mu = mean(X, 1);
    Xc = X - mu;
    
    % PCA via SVD (no toolbox dependency)
    [~, S, V] = svd(Xc, 'econ');
    
    pc1   = V(:,1);              % [nFreq x 1] dominant spectral pattern
    score1 = Xc * pc1;            % [(T*R) x 1]
    
    % remove PC1 (rank-1 reconstruction)
    Xc_den = Xc - score1 * pc1';
    
    % optional: explained variance of PC1
    singvals = diag(S);
    expVar = (singvals.^2) / sum(singvals.^2);
    expVar1 = expVar(1);
    
    % reshape back to [F x T x R]
    X_den = Xc_den + mu;
    ERSP_den = permute(reshape(X_den, nTime, nTrial, nFreq), [3 1 2]);
    
    % package outputs
    out = struct();
    out.trialIDs    = trialIDs;
    out.freqs       = freqs;
    out.times       = times;
    out.baseline_ms = baseline_ms;
    
    out.ERSP_trial_raw = ERSP_trial;   % [F x T x R]
    out.ERSP_trial_den = ERSP_den;     % [F x T x R]
    
    out.PC1 = pc1;                      % [F x 1]
    out.expVar = expVar;                % vector
    out.expVar1 = expVar1;              % scalar

end