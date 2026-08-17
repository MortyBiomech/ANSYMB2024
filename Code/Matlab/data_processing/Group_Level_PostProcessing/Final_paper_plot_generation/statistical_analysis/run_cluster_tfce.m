% -------------------------------------------------------------------------
% CLUSTER-BASED PERMUTATION LMM + TFCE ANALYSIS FOR EEG TFR DATA
% -------------------------------------------------------------------------
%  author : Morteza Khosrotabar
%  date   : 2025-08-12
%
%  This script loops over 9 predefined brain-IC clusters and, for each:
%    1.  Builds an (freq × time) map of F-statistics (main effect of Pressure)
%        using a linear mixed-effects model:           Power ~ Pressure + (1|SubjectID)
%    2.  Enhances that map using Threshold-Free Cluster Enhancement (TFCE)
%    3.  Generates an empirical null distribution of maximum TFCE values
%        via within-IC permutation of condition labels
%    4.  Computes a cluster-wise family-wise-error-controlled significance mask
%
%  REQUIREMENTS
%  -----------------------------------------------------------------------
%  MATLAB R2018a+  (Statistics & Machine Learning Toolbox)
%  No third-party toolboxes are required; a vectorised 2-D TFCE implementation
%  is included below (function tfce2d).
%
%  USER PARAMETERS
%  -----------------------------------------------------------------------
p_value          = 0.05;     % FWE-corrected alpha level
num_permutations = 1000;     % ≥1000 recommended
tfce_H           = 2;        % TFCE height exponent
tfce_E           = 0.5;      % TFCE extent exponent
dh               = 0.1;      % TFCE integration step
fprintf('Starting cluster analysis (%d permutations, α = %.3f)\n\n', ...
        num_permutations, p_value);

% -------------------------------------------------------------------------
% MAIN LOOP OVER 9 CLUSTERS
% -------------------------------------------------------------------------
for clusterIdx = 1:9
    fprintf('==== CLUSTER %d ============================================\n', clusterIdx);

    % === (A) LOAD DATA FOR THIS CLUSTER =================================
    %  Replace the following three lines with your own loader if needed.
    %  Each variable must be a cell array with one cell per IC.
    % --------------------------------------------------------------------
    load(sprintf('cluster%d_data.mat',clusterIdx), ...          % <-- YOUR FILE
         'cluster_time_frequency_data', ...   % cell{nIC}(freq × time × trials)
         'cluster_ic_info', ...               % struct array, .subject_id
         'cluster_trial_conditions');         % cell{nIC}(1 × trials) values∈[1 3 6]

    nIC            = size(cluster_time_frequency_data, 1);
    [nFreq,nTime]  = size(cluster_time_frequency_data{1, 3}.P1{1, 1}); %#ok<NASGU>
    observed_Fmap  = zeros(nFreq,nTime);

    % === (B) LOOP OVER ALL (FREQ,TIME) SAMPLES ==========================
    for f = 1:nFreq
        for t = 1:nTime
            % Build table for this (f,t)
            tbl = build_datatable(cluster_time_frequency_data, ...
                                  cluster_trial_conditions, ...
                                  cluster_ic_info, f, t);
            % Ensure categorical predictors
            tbl.Pressure  = categorical(tbl.Pressure);
            tbl.SubjectID = categorical(tbl.SubjectID);

            % Fit LMM
            lme  = fitlme(tbl,'Power ~ Pressure + (1|SubjectID)');
            aov  = anova(lme,'DFMethod','Satterthwaite');
            observed_Fmap(f,t) = aov.FStat(2);  % main effect of Pressure
        end
    end

    % === (C) APPLY TFCE TO OBSERVED F-MAP ===============================
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);

    % === (D) PERMUTATION LOOP ===========================================
    max_TFCE_null = zeros(num_permutations,1);
    parfor perm = 1:num_permutations                      % use PARPOOL if available
        perm_Fmap = zeros(nFreq,nTime);

        % Shuffle labels *within each IC* --------------------------------
        perm_conditions = cellfun(@(c) c(randperm(numel(c))), ...
                                  cluster_trial_conditions, ...
                                  'UniformOutput',false);

        % Re-compute F-map under permuted labels
        for f = 1:nFreq
            for t = 1:nTime
                tbl = build_datatable(cluster_time_frequency_data, ...
                                      perm_conditions, ...
                                      cluster_ic_info, f, t);
                tbl.Pressure  = categorical(tbl.Pressure);
                tbl.SubjectID = categorical(tbl.SubjectID);
                lme  = fitlme(tbl,'Power ~ Pressure + (1|SubjectID)');
                aov  = anova(lme,'DFMethod','Satterthwaite');
                perm_Fmap(f,t) = aov.FStat(2);
            end
        end

        perm_TFCE         = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
        max_TFCE_null(perm) = max(perm_TFCE(:));           % store max statistic
    end

    % === (E) THRESHOLD & SIGNIFICANCE MASK ==============================
    tfce_threshold              = prctile(max_TFCE_null, 100*(1-p_value));
    significance_mask           = observed_TFCE > tfce_threshold;

    % === (F) SAVE RESULTS ===============================================
    save(sprintf('cluster%d_results.mat',clusterIdx), ...
         'observed_Fmap','observed_TFCE','significance_mask', ...
         'tfce_threshold','max_TFCE_null');

    fprintf('   Significant TFCE voxels: %d / %d\n\n', ...
            nnz(significance_mask), numel(significance_mask));
end

fprintf('Analysis complete.\n');

% ========================================================================
% FUNCTION: build_datatable
% ------------------------------------------------------------------------
function tbl = build_datatable(tf_data, cond_labels, ic_info, fIdx, tIdx)
% Assemble trial-wise data table for a given (freq,time) sample across all ICs
    power_vals  = [];
    press_vals  = [];
    subj_vals   = [];
    nIC         = numel(tf_data);

    for ic = 1:nIC
        P = squeeze(tf_data{ic}(fIdx,tIdx,:));      % column vector (trials × 1)
        C = cond_labels{ic}(:);                     % column vector (trials × 1)
        S = ic_info(ic).subject_id;                 % scalar

        power_vals = [power_vals; P];               %#ok<AGROW>
        press_vals = [press_vals; C];
        subj_vals  = [subj_vals; repmat(S,numel(P),1)];
    end

    tbl = table(power_vals, press_vals, subj_vals, ...
                'VariableNames',{'Power','Pressure','SubjectID'});
end

% ========================================================================
% FUNCTION: tfce2d
% ------------------------------------------------------------------------
function TFCE = tfce2d(statMap, H, E, dh)
% Vectorised 2-D TFCE implementation (Smith & Nichols, 2009)
%
% INPUT  statMap : 2-D matrix of statistic values (e.g., F)
%        H,E     : TFCE exponents (default H=2, E=0.5)
%        dh      : integration step over h (default 0.1)
% OUTPUT TFCE    : TFCE-transformed statistic map
%
    if nargin<4, dh = 0.1; end
    if nargin<3, E  = 0.5; end
    if nargin<2, H  = 2;   end

    statMap(statMap<0) = 0;             % TFCE defined for positive stats
    hVals  = 0:dh:max(statMap(:));
    TFCE   = zeros(size(statMap));

    % 26-connected neighbourhood (8-conn in 2-D)
    conn   = conndef(2,'maximal');
    for h = hVals
        mask         = statMap >= h;            % threshold at height h
        CC           = bwconncomp(mask,conn);   % identify clusters
        for clIdx = 1:CC.NumObjects
            voxIdx   = CC.PixelIdxList{clIdx};
            extent   = numel(voxIdx);
            TFCE(voxIdx) = TFCE(voxIdx) + ...
                            (h^H) * (extent^E) * dh;
        end
    end
end
