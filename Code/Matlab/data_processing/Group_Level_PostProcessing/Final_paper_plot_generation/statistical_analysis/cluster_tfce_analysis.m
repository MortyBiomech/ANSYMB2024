% -------------------------------------------------------------------------
% CLUSTER-BASED PERMUTATION LMM + TFCE ANALYSIS FOR EEG TFR DATA
% -------------------------------------------------------------------------
%  author : EEG Analysis Expert
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
%  DATA STRUCTURE:
%  cluster_TF_data_warped: Cell array where each row is one IC with columns:
%    Column 1: Subject ID
%    Column 2: IC ID
%    Column 3: Struct with fields P1, P3, P6 (each containing trial matrices)
%    Column 4: Additional data (frequencies, etc.)
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
    %  Replace the following line with your own loader if needed.
    %  The variable should be a cell array with structure as described above
    % --------------------------------------------------------------------
    load(sprintf('cluster%d_TF_data_warped.mat', clusterIdx), 'cluster_TF_data_warped');
    
    % Extract dimensions from first IC's first trial
    nIC = size(cluster_TF_data_warped, 1);
    first_trial = cluster_TF_data_warped{1,3}.P1{1};  % Get first trial of P1 condition
    [nFreq, nTime] = size(first_trial);
    
    observed_Fmap = zeros(nFreq, nTime);
    
    fprintf('   Processing %d ICs, %d frequencies, %d time points\n', nIC, nFreq, nTime);

    % === (B) LOOP OVER ALL (FREQ,TIME) SAMPLES ==========================
    fprintf('   Computing F-statistics across time-frequency space...\n');
    parfor f = 1:nFreq
        % if mod(f, 20) == 0
        %     fprintf('   Processing frequency %d/%d\n', f, nFreq);
        % end
        
        for t = 1:nTime
            % Build table for this (f,t) sample
            tbl = build_datatable_new(cluster_TF_data_warped, f, t);
            
            % Ensure categorical predictors
            tbl.Pressure  = categorical(tbl.Pressure);
            tbl.SubjectID = categorical(tbl.SubjectID);
            tbl.IC_ID = categorical(tbl.IC_ID);

            % Create nested IC identifier (unique across subjects)
    	    tbl.IC_nested = categorical(strcat(string(tbl.SubjectID), '_', string(tbl.IC_ID)));

            % Fit LMM and extract F-statistic
            try
                lme = fitlme(tbl, 'Power ~ Pressure + (1|IC_nested)');
                aov = anova(lme, 'DFMethod', 'Satterthwaite');
                observed_Fmap(f,t) = aov.FStat(2);  % main effect of Pressure
            catch ME
                % Handle cases where LMM fitting fails
                observed_Fmap(f,t) = 0;
                if f == 1 && t == 1
                    warning('LMM fitting failed at (%d,%d): %s', f, t, ME.message);
                end
            end
        end
    end

    % === (C) APPLY TFCE TO OBSERVED F-MAP ===============================
    fprintf('   Applying TFCE to observed F-map...\n');
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);

    % === (D) PERMUTATION LOOP ===========================================
    fprintf('   Running %d permutations...\n', num_permutations);
    max_TFCE_null = zeros(num_permutations, 1);
    
    % Use parallel processing if Parallel Computing Toolbox is available
    if exist('parfor', 'builtin') && ~isempty(gcp('nocreate'))
        use_parfor = true;
    else
        use_parfor = false;
        fprintf('   (Using regular for loop - consider starting parallel pool for speed)\n');
    end
    
    if use_parfor
        parfor perm = 1:num_permutations
            perm_Fmap = compute_permuted_fmap(cluster_TF_data_warped, nFreq, nTime);
            perm_TFCE = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
            max_TFCE_null(perm) = max(perm_TFCE(:));
        end
    else
        for perm = 1:num_permutations
            if mod(perm, 100) == 0
                fprintf('   Permutation %d/%d\n', perm, num_permutations);
            end
            perm_Fmap = compute_permuted_fmap(cluster_TF_data_warped, nFreq, nTime);
            perm_TFCE = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
            max_TFCE_null(perm) = max(perm_TFCE(:));
        end
    end

    % === (E) THRESHOLD & SIGNIFICANCE MASK ==============================
    tfce_threshold = prctile(max_TFCE_null, 100*(1-p_value));
    significance_mask = observed_TFCE > tfce_threshold;

    % === (F) SAVE RESULTS ===============================================
    save(sprintf('cluster%d_results.mat', clusterIdx), ...
         'observed_Fmap', 'observed_TFCE', 'significance_mask', ...
         'tfce_threshold', 'max_TFCE_null', 'nFreq', 'nTime', 'nIC');

    fprintf('   TFCE threshold: %.4f\n', tfce_threshold);
    fprintf('   Significant TFCE voxels: %d / %d (%.2f%%)\n\n', ...
            nnz(significance_mask), numel(significance_mask), ...
            100*nnz(significance_mask)/numel(significance_mask));
end

fprintf('Analysis complete.\n');

% ========================================================================
% FUNCTION: build_datatable_new
% ------------------------------------------------------------------------
% function tbl = build_datatable_new(cluster_data, fIdx, tIdx)
% % Assemble trial-wise data table for a given (freq,time) sample across all ICs
% % 
% % INPUT:
% %   cluster_data: Cell array where each row represents one IC
% %                 Column 1: Subject ID
% %                 Column 2: IC ID
% %                 Column 3: Struct with fields P1, P3, P6
% %   fIdx, tIdx:   Frequency and time indices
% %
% % OUTPUT:
% %   tbl: Table with columns Power, Pressure, SubjectID
% 
%     power_vals = [];
%     press_vals = [];
%     subj_vals = [];
%     ic_vals = [];
% 
%     nIC = size(cluster_data, 1);
%     condition_names = {'P1', 'P3', 'P6'};
%     condition_values = [1, 3, 6];  % Map P1->1, P3->3, P6->6
% 
%     for ic = 1:nIC
%         subject_id = cluster_data{ic, 1};  % Column 1: Subject ID
%         ic_id = cluster_data{ic, 2};       % Column 2: Subject IC id
%         tf_struct = cluster_data{ic, 3};   % Column 3: TF data struct
% 
%         % Loop through each condition (P1, P3, P6)
%         for cond_idx = 1:length(condition_names)
%             cond_name = condition_names{cond_idx};
%             cond_value = condition_values(cond_idx);
% 
%             % Check if this condition exists for this IC
%             if isfield(tf_struct, cond_name) && ~isempty(tf_struct.(cond_name))
%                 trials = tf_struct.(cond_name);  % Cell array of trial matrices
%                 nTrials = length(trials);
% 
%                 % Extract power values for this (freq, time) point across all trials
%                 trial_powers = zeros(nTrials, 1);
%                 for trial = 1:nTrials
%                     trial_matrix = trials{trial};  % 183 x 1001 matrix
%                     trial_powers(trial) = trial_matrix(fIdx, tIdx);
%                 end
% 
%                 % Append to main arrays
%                 power_vals = [power_vals; trial_powers]; %#ok<AGROW>
%                 press_vals = [press_vals; repmat(cond_value, nTrials, 1)]; %#ok<AGROW>
%                 subj_vals = [subj_vals; repmat(subject_id, nTrials, 1)]; %#ok<AGROW>
%                 ic_vals = [ic_vals; repmat(ic_id, nTrials, 1)]; %#ok<AGROW>
%             end
%         end
%     end
% 
%     % Create the table
%     tbl = table(power_vals, press_vals, subj_vals, ic_vals, ...
%                 'VariableNames', {'Power', 'Pressure', 'SubjectID', 'IC_ID'});
% end

% ========================================================================
% FUNCTION: compute_permuted_fmap
% ------------------------------------------------------------------------
function perm_Fmap = compute_permuted_fmap(cluster_data, nFreq, nTime)
% Compute F-statistic map with permuted condition labels
% 
% Permutation strategy: For each IC, randomly shuffle the condition labels
% across all trials, while maintaining the number of trials per condition

    perm_Fmap = zeros(nFreq, nTime);
    
    % Create permuted version of cluster data
    perm_cluster_data = cluster_data;
    nIC = size(cluster_data, 1);
    
    for ic = 1:nIC
        tf_struct = cluster_data{ic, 3};
        condition_names = {'P1', 'P3', 'P6'};
        
        % Collect all trials and their original labels
        all_trials = {};
        original_labels = [];
        
        for cond_idx = 1:length(condition_names)
            cond_name = condition_names{cond_idx};
            if isfield(tf_struct, cond_name) && ~isempty(tf_struct.(cond_name))
                trials = tf_struct.(cond_name);
                nTrials = length(trials);
                all_trials = [all_trials, trials]; %#ok<AGROW>
                original_labels = [original_labels, repmat(cond_idx, 1, nTrials)]; %#ok<AGROW>
            end
        end
        
        % Shuffle the labels
        shuffled_labels = original_labels(randperm(length(original_labels)));
        
        % Reassign trials to conditions based on shuffled labels
        new_struct = struct();
        for cond_idx = 1:length(condition_names)
            cond_name = condition_names{cond_idx};
            trial_indices = find(shuffled_labels == cond_idx);
            if ~isempty(trial_indices)
                new_struct.(cond_name) = all_trials(trial_indices);
            else
                new_struct.(cond_name) = {};
            end
        end
        
        perm_cluster_data{ic, 3} = new_struct;
    end
    
    % Compute F-statistics for permuted data
    for f = 1:nFreq
        for t = 1:nTime
            tbl = build_datatable_new(perm_cluster_data, f, t);
            
            % Handle case where permutation creates empty conditions
            if height(tbl) < 3 || length(unique(tbl.Pressure)) < 2
                perm_Fmap(f,t) = 0;
                continue;
            end
            
            tbl.Pressure = categorical(tbl.Pressure);
            tbl.SubjectID = categorical(tbl.SubjectID);
            
            try
                lme = fitlme(tbl, 'Power ~ Pressure + (1|SubjectID)');
                aov = anova(lme, 'DFMethod', 'Satterthwaite');
                perm_Fmap(f,t) = aov.FStat(2);
            catch
                perm_Fmap(f,t) = 0;
            end
        end
    end
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
    if nargin < 4, dh = 0.1; end
    if nargin < 3, E = 0.5; end
    if nargin < 2, H = 2; end

    % Ensure non-negative values (TFCE defined for positive stats)
    statMap(statMap < 0) = 0;
    
    if max(statMap(:)) == 0
        TFCE = zeros(size(statMap));
        return;
    end
    
    % Define height thresholds
    max_stat = max(statMap(:));
    hVals = 0:dh:max_stat;
    TFCE = zeros(size(statMap));

    % 8-connected neighbourhood for 2D
    conn = conndef(2, 'maximal');
    
    for h = hVals
        % Threshold map at height h
        mask = statMap >= h;
        
        if ~any(mask(:))
            continue;
        end
        
        % Find connected components
        CC = bwconncomp(mask, conn);
        
        % Process each cluster
        for clIdx = 1:CC.NumObjects
            voxIdx = CC.PixelIdxList{clIdx};
            extent = numel(voxIdx);
            
            % Add TFCE contribution: h^H * extent^E * dh
            TFCE(voxIdx) = TFCE(voxIdx) + (h^H) * (extent^E) * dh;
        end
    end
end