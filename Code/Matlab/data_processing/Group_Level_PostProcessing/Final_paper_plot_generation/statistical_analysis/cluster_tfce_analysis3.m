% -------------------------------------------------------------------------
% TRIAL-LEVEL SHUFFLING NESTED REPEATED MEASURES ANOVA FOR EEG TF DATA
% -------------------------------------------------------------------------
% This script performs cluster-based permutation testing on EEG time-frequency
% data using repeated measures ANOVA with trial-level permutation shuffling.
% 
% Key features:
% - Stores both original trials AND averaged data
% - Trial-level shuffling for statistically valid permutation testing
% - Simple nested ANOVA: Power ~ NestedIC + Pressure + NestedIC×Pressure
% - NestedIC as random factor (implicitly handles subject nesting)
% - Applies TFCE for spatial enhancement
% - Full resolution analysis (no segmentation)
%
% Requirements: Statistics and Machine Learning Toolbox
% -------------------------------------------------------------------------

% ANALYSIS PARAMETERS
p_value = 0.05;                 % Significance level
num_permutations = 1000;        % Number of permutations
tfce_H = 2;                     % TFCE height exponent
tfce_E = 0.5;                   % TFCE extent exponent
dh = 0.1;                       % TFCE integration step

fprintf('Starting trial-level shuffling nested repeated measures ANOVA analysis...\n');
fprintf('Model: Power ~ NestedIC + Pressure + NestedIC×Pressure (NestedIC as random)\n');
fprintf('Shuffling: Trial-level permutation for valid null distribution\n');
fprintf('Permutations: %d, Alpha: %.3f\n\n', num_permutations, p_value);

% -------------------------------------------------------------------------
% MAIN LOOP OVER 9 CLUSTERS
% -------------------------------------------------------------------------
for clusterIdx = 1:9
    % fprintf('==== CLUSTER %d (Trial-Level Shuffling) =================\n', clusterIdx);
    % 
    % % Load data for this cluster
    % load(sprintf('cluster%d_TF_data_warped.mat', clusterIdx), 'cluster_TF_data_warped');

    % Get dimensions from the data
    [nFreq, nTime] = size(cluster_TF_data_warped{1,3}.P1{1});
    nIC = size(cluster_TF_data_warped, 1);
    
    fprintf('Data dimensions: %d ICs × %d freq × %d time = %d analysis points\n', ...
            nIC, nFreq, nTime, nFreq*nTime);
    
    % Pre-compute data structure with both trials and averages
    fprintf('Preparing data structure with original trials and averages...\n');
    data_struct = prepare_trial_and_averaged_data(cluster_TF_data_warped, nFreq, nTime);
    
    % % Print data structure info
    % print_trial_level_structure_info(data_struct);
    
    % Compute observed F-map using trial averages
    fprintf('Computing observed F-statistics using averaged data...\n');
    observed_Fmap = compute_fmap_from_averages(data_struct, nFreq, nTime);
    
    % Apply TFCE to observed F-map
    fprintf('Applying TFCE to observed F-map...\n');
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);
    
    % Permutation testing with trial-level shuffling
    fprintf('Running %d permutations with trial-level shuffling...\n', num_permutations);
    max_TFCE_null = zeros(num_permutations, 1);
    
    parfor perm = 1:num_permutations
        % Create permuted data by shuffling at trial level, then averaging
        perm_averages = create_trial_level_permuted_averages(data_struct);
        
        % Compute F-map for permuted averaged data
        perm_Fmap = compute_fmap_from_averages(perm_averages, nFreq, nTime);
        
        % Apply TFCE to permuted F-map
        perm_TFCE = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
        
        % Store maximum TFCE value
        max_TFCE_null(perm) = max(perm_TFCE(:));
    end
    
    % Determine significance threshold
    tfce_threshold = prctile(max_TFCE_null, 100*(1-p_value));
    
    % Generate final significance mask
    significance_mask = observed_TFCE > tfce_threshold;
    
    % Save results for this cluster
    save(sprintf('cluster%d_trial_level_results.mat', clusterIdx), ...
         'observed_Fmap', 'observed_TFCE', 'significance_mask', ...
         'tfce_threshold', 'max_TFCE_null', 'nFreq', 'nTime', 'nIC');
    
    fprintf('Significant time-frequency points: %d / %d (%.1f%%)\n\n', ...
            nnz(significance_mask), numel(significance_mask), ...
            100*nnz(significance_mask)/numel(significance_mask));
end

fprintf('Trial-level shuffling analysis complete.\n');

% =========================================================================
% SUPPORTING FUNCTIONS
% =========================================================================

% % ========================================================================
% % FUNCTION: prepare_trial_and_averaged_data
% % ------------------------------------------------------------------------
% function data_struct = prepare_trial_and_averaged_data(cluster_TF_data_warped, nFreq, nTime)
% 
%     nIC = size(cluster_TF_data_warped, 1);
%     conditions = {'P1', 'P3', 'P6'};
% 
%     % Initialize data structure
%     data_struct = struct();
%     data_struct.subject_ids = [];
%     data_struct.ic_ids = [];
%     data_struct.nested_ic_ids = {};
%     data_struct.trial_data = {};      % Store original trials for permutation
%     data_struct.averaged_data = {};   % Store averaged data for analysis
%     data_struct.trial_counts = {};    % Store number of trials per condition
% 
%     valid_ic_count = 0;
% 
%     for ic = 1:nIC
%         subject_id = cluster_TF_data_warped{ic, 1};
%         ic_id = cluster_TF_data_warped{ic, 2};
%         tf_struct = cluster_TF_data_warped{ic, 3};
% 
%         % Create nested IC identifier
%         nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
% 
%         % Check if this IC has data for all conditions
%         has_all_conditions = true;
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             if ~isfield(tf_struct, cond_name) || isempty(tf_struct.(cond_name))
%                 has_all_conditions = false;
%                 break;
%             end
%         end
% 
%         if has_all_conditions
%             valid_ic_count = valid_ic_count + 1;
%             data_struct.subject_ids(valid_ic_count) = subject_id;
%             data_struct.ic_ids(valid_ic_count) = ic_id;
%             data_struct.nested_ic_ids{valid_ic_count} = nested_ic_id;
% 
%             % Initialize structures for this IC
%             data_struct.trial_data{valid_ic_count} = struct();
%             data_struct.averaged_data{valid_ic_count} = struct();
%             data_struct.trial_counts{valid_ic_count} = struct();
% 
%             % Process each condition
%             for c = 1:length(conditions)
%                 cond_name = conditions{c};
%                 trials = tf_struct.(cond_name);
%                 n_trials = length(trials);
% 
%                 % Store original trials (nFreq × nTime × nTrials)
%                 trial_stack = zeros(nFreq, nTime, n_trials);
%                 for trial = 1:n_trials
%                     trial_stack(:, :, trial) = trials{trial};
%                 end
%                 data_struct.trial_data{valid_ic_count}.(cond_name) = trial_stack;
%                 data_struct.trial_counts{valid_ic_count}.(cond_name) = n_trials;
% 
%                 % Compute and store averages
%                 averaged_trial = mean(trial_stack, 3, 'omitnan');
%                 data_struct.averaged_data{valid_ic_count}.(cond_name) = averaged_trial;
%             end
%         end
%     end
% 
%     % Store final count
%     data_struct.nIC = valid_ic_count;
% 
%     % Trim arrays to actual size
%     data_struct.subject_ids = data_struct.subject_ids(1:valid_ic_count);
%     data_struct.ic_ids = data_struct.ic_ids(1:valid_ic_count);
%     data_struct.nested_ic_ids = data_struct.nested_ic_ids(1:valid_ic_count);
%     data_struct.trial_data = data_struct.trial_data(1:valid_ic_count);
%     data_struct.averaged_data = data_struct.averaged_data(1:valid_ic_count);
%     data_struct.trial_counts = data_struct.trial_counts(1:valid_ic_count);
% end

% % ========================================================================
% % FUNCTION: print_trial_level_structure_info
% % ------------------------------------------------------------------------
% function print_trial_level_structure_info(data_struct)
%     fprintf('Trial-level data structure:\n');
%     fprintf('  Total valid ICs: %d\n', data_struct.nIC);
% 
%     unique_subjects = unique(data_struct.subject_ids);
%     fprintf('  Unique subjects: %d\n', length(unique_subjects));
% 
%     conditions = {'P1', 'P3', 'P6'};
%     total_trials = 0;
% 
%     % Calculate trial statistics
%     for ic = 1:min(5, data_struct.nIC)  % Show first 5 ICs
%         fprintf('  IC %s: ', data_struct.nested_ic_ids{ic});
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             n_trials = data_struct.trial_counts{ic}.(cond_name);
%             fprintf('%s=%d ', cond_name, n_trials);
%             total_trials = total_trials + n_trials;
%         end
%         fprintf('trials\n');
%     end
% 
%     if data_struct.nIC > 5
%         fprintf('  ... and %d more ICs\n', data_struct.nIC - 5);
%     end
% 
%     fprintf('  Average ICs per subject: %.1f\n', data_struct.nIC / length(unique_subjects));
%     fprintf('  Total trials across all ICs: %d\n\n', total_trials);
% end

% % ========================================================================
% % FUNCTION: create_trial_level_permuted_averages
% % ------------------------------------------------------------------------
% function perm_averages = create_trial_level_permuted_averages(data_struct)
% 
%     % Initialize permuted averages structure
%     perm_averages = struct();
%     perm_averages.subject_ids = data_struct.subject_ids;
%     perm_averages.ic_ids = data_struct.ic_ids;
%     perm_averages.nested_ic_ids = data_struct.nested_ic_ids;
%     perm_averages.nIC = data_struct.nIC;
%     perm_averages.averaged_data = cell(data_struct.nIC, 1);
% 
%     conditions = {'P1', 'P3', 'P6'};
%     condition_values = [1, 3, 6];
% 
%     % For each IC, shuffle trial labels and recompute averages
%     for ic = 1:data_struct.nIC
%         % Initialize averaged data structure for this IC
%         perm_averages.averaged_data{ic} = struct();
% 
%         % Collect all trials and their original condition labels
%         all_trials = [];
%         all_labels = [];
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             cond_val = condition_values(c);
% 
%             trials = data_struct.trial_data{ic}.(cond_name);  % nFreq × nTime × nTrials
%             n_trials = size(trials, 3);
% 
%             % Concatenate trials along 3rd dimension
%             all_trials = cat(3, all_trials, trials);
% 
%             % Create condition labels for these trials
%             all_labels = [all_labels; repmat(cond_val, n_trials, 1)]; %#ok<AGROW>
%         end
% 
%         % Shuffle the condition labels (permutation step)
%         total_trials = length(all_labels);
%         shuffled_labels = all_labels(randperm(total_trials));
% 
%         % Recompute averages based on shuffled labels
%         for c = 1:length(conditions)
%             cond_val = condition_values(c);
% 
%             % Find trials that now belong to this condition after shuffling
%             condition_trial_indices = find(shuffled_labels == cond_val);
% 
%             if ~isempty(condition_trial_indices)
%                 % Extract trials for this condition and compute average
%                 condition_trials = all_trials(:, :, condition_trial_indices);
%                 perm_averages.averaged_data{ic}.(conditions{c}) = mean(condition_trials, 3, 'omitnan');
%             else
%                 % Handle case where no trials assigned to this condition
%                 % Use NaN array as placeholder
%                 [nFreq, nTime] = size(all_trials(:, :, 1));
%                 perm_averages.averaged_data{ic}.(conditions{c}) = nan(nFreq, nTime);
%             end
%         end
%     end
% end

% % ========================================================================
% % FUNCTION: compute_fmap_from_averages
% % ------------------------------------------------------------------------
% function observed_Fmap = compute_fmap_from_averages(data_struct, nFreq, nTime)
% 
%     observed_Fmap = zeros(nFreq, nTime);
%     conditions = {'P1', 'P3', 'P6'};
% 
%     parfor f = 1:nFreq  % Parallel over frequencies
%         temp_row = zeros(1, nTime);
%         for t = 1:nTime
%             % Build data table for this (f,t) point using averaged data
%             tbl = build_table_from_averages(data_struct, f, t, conditions);
% 
%             if height(tbl) >= 6  % Need minimum data points for ANOVA
%                 try
%                     % Convert to categorical variables
%                     tbl.NestedIC = categorical(tbl.NestedIC);
%                     tbl.Pressure = categorical(tbl.Pressure);
% 
%                     % Perform simple nested repeated measures ANOVA
%                     [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
%                                           'model', 'interaction', ...
%                                           'random', 1, ...  % NestedIC as random factor
%                                           'display', 'off');
% 
%                     temp_row(t) = cell2mat(T(3, 6));
% 
%                 catch ME
%                     % Handle ANOVA failures
%                     if contains(ME.message, 'singular') || contains(ME.message, 'rank deficient')
%                         try
%                             % Try simpler linear model if interaction fails
%                             [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
%                                                   'model', 'linear', ...
%                                                   'random', 1, ...
%                                                   'display', 'off');
%                             temp_row(t) = cell2mat(T(3, 6));
%                         catch
%                             temp_row(t) = 0;
%                         end
%                     else
%                         temp_row(t) = 0;
%                     end
%                 end
%             end
%         end
%         observed_Fmap(f, :) = temp_row;
%     end
% end

% % ========================================================================
% % FUNCTION: build_table_from_averages
% % ------------------------------------------------------------------------
% function tbl = build_table_from_averages(data_struct, f, t, conditions)
% 
%     power_vals = [];
%     pressure_vals = [];
%     nested_ic_vals = {};
% 
%     nIC = data_struct.nIC;
%     condition_values = [1, 3, 6];
% 
%     for ic = 1:nIC
%         nested_ic_id = data_struct.nested_ic_ids{ic};
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             cond_val = condition_values(c);
% 
%             % Extract averaged power value for this IC at (f,t)
%             if isfield(data_struct.averaged_data{ic}, cond_name)
%                 power_val = data_struct.averaged_data{ic}.(cond_name)(f, t);
% 
%                 % Only include finite values
%                 if isfinite(power_val)
%                     power_vals = [power_vals; power_val]; %#ok<AGROW>
%                     pressure_vals = [pressure_vals; cond_val]; %#ok<AGROW>
%                     nested_ic_vals = [nested_ic_vals; {nested_ic_id}]; %#ok<AGROW>
%                 end
%             end
%         end
%     end
% 
%     tbl = table(power_vals, pressure_vals, nested_ic_vals, ...
%                 'VariableNames', {'Power', 'Pressure', 'NestedIC'});
% end

% % ========================================================================
% % FUNCTION: tfce2d - Threshold-Free Cluster Enhancement
% % ------------------------------------------------------------------------
% function TFCE = tfce2d(statMap, H, E, dh)
%     % Vectorised 2-D TFCE implementation (Smith & Nichols, 2009)
% 
%     if nargin<4, dh = 0.1; end
%     if nargin<3, E = 0.5; end
%     if nargin<2, H = 2; end
% 
%     statMap(statMap<0) = 0;             % TFCE defined for positive stats
%     statMap(~isfinite(statMap)) = 0;    % Handle NaN and Inf values
% 
%     if max(statMap(:)) == 0
%         TFCE = zeros(size(statMap));
%         return;
%     end
% 
%     hVals = 0:dh:max(statMap(:));
%     TFCE = zeros(size(statMap));
% 
%     % 8-connected neighbourhood in 2D
%     conn = conndef(2,'maximal');
% 
%     for h = hVals
%         if h == 0, continue; end  % Skip h=0 to avoid empty clusters
% 
%         mask = statMap >= h;            % Threshold at height h
%         if ~any(mask(:)), continue; end % Skip if no voxels survive threshold
% 
%         CC = bwconncomp(mask, conn);    % Identify clusters
% 
%         for clIdx = 1:CC.NumObjects
%             voxIdx = CC.PixelIdxList{clIdx};
%             extent = numel(voxIdx);
%             if extent > 0  % Safety check
%                 TFCE(voxIdx) = TFCE(voxIdx) + (h^H) * (extent^E) * dh;
%             end
%         end
%     end
% end
