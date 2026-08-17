% -------------------------------------------------------------------------
% SEGMENTED TIME-FREQUENCY LMM ANALYSIS WITH NESTED IC STRUCTURE
% -------------------------------------------------------------------------
% This script performs cluster-based permutation testing on EEG time-frequency
% data using Linear Mixed-Effects Models with nested IC structure and
% Threshold-Free Cluster Enhancement (TFCE).
%
% Key features:
% - Reduces computational load by averaging within time-frequency segments
% - Uses nested IC identifiers (SubjectID_ICID) as random effects
% - LMM: Power ~ Pressure + (1|SubjectID) + (1|NestedIC)
% - Applies TFCE for spatial enhancement
% - Permutation testing for family-wise error control
%
% Requirements: Statistics and Machine Learning Toolbox
% -------------------------------------------------------------------------

% ANALYSIS PARAMETERS
p_value = 0.05;                 % Significance level
num_permutations = 1000;        % Number of permutations
tfce_H = 2;                     % TFCE height exponent
tfce_E = 0.5;                   % TFCE extent exponent
dh = 0.1;                       % TFCE integration step

% RESOLUTION REDUCTION PARAMETERS
time_reduction_factor = 0.1;    % Use 10% of time points
freq_reduction_factor = 0.1;    % Use 10% of frequency points

fprintf('Starting nested IC segmented cluster analysis...\n');
fprintf('Model: Power ~ Pressure + (1|SubjectID) + (1|NestedIC)\n');
fprintf('Time reduction: %.1f%%, Frequency reduction: %.1f%%\n', ...
        time_reduction_factor*100, freq_reduction_factor*100);
fprintf('Permutations: %d, Alpha: %.3f\n\n', num_permutations, p_value);

% -------------------------------------------------------------------------
% MAIN LOOP OVER 9 CLUSTERS
% -------------------------------------------------------------------------
for clusterIdx = 1:9
    % fprintf('==== CLUSTER %d (Nested IC Analysis) ================\n', clusterIdx);
    % 
    % % Load data for this cluster
    % load(sprintf('cluster%d_TF_data_warped.mat', clusterIdx), 'cluster_TF_data_warped');
    
    % Get original dimensions
    [orig_nFreq, orig_nTime] = size(cluster_TF_data_warped{1,3}.P1{1});
    
    % Calculate segment dimensions
    [seg_nFreq, seg_nTime, freq_segments, time_segments] = ...
        calculate_segment_dimensions(orig_nFreq, orig_nTime, ...
                                   freq_reduction_factor, time_reduction_factor);
    
    fprintf('Original: %d freq × %d time = %d points\n', orig_nFreq, orig_nTime, orig_nFreq*orig_nTime);
    fprintf('Segmented: %d freq × %d time = %d points (%.1fx reduction)\n', ...
            seg_nFreq, seg_nTime, seg_nFreq*seg_nTime, (orig_nFreq*orig_nTime)/(seg_nFreq*seg_nTime));
    
    % Pre-compute segmented data with IC identifiers
    fprintf('Pre-computing segmented data with nested IC identifiers...\n');
    segmented_data = precompute_segmented_data_with_ic(cluster_TF_data_warped, freq_segments, time_segments);
    
    % Print nested IC structure info
    print_nested_ic_info(segmented_data);
    
    % Compute observed F-map with nested IC model
    fprintf('Computing observed F-statistics with nested IC model...\n');
    observed_Fmap = compute_segmented_fmap_nested(segmented_data, seg_nFreq, seg_nTime);
    
    % Apply TFCE to observed F-map
    fprintf('Applying TFCE to observed F-map...\n');
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);
    
    % Permutation testing with nested structure
    fprintf('Running %d permutations with nested IC structure...\n', num_permutations);
    max_TFCE_null = zeros(num_permutations, 1);
    
    parfor perm = 1:num_permutations
        % Create permuted data (shuffle within each IC)
        perm_data = create_permuted_segmented_data_nested(segmented_data);
        
        % Compute F-map for permuted data
        perm_Fmap = compute_segmented_fmap_nested(perm_data, seg_nFreq, seg_nTime);
        
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
    save(sprintf('cluster%d_nested_results.mat', clusterIdx), ...
         'observed_Fmap', 'observed_TFCE', 'significance_mask', ...
         'tfce_threshold', 'max_TFCE_null', 'freq_segments', 'time_segments', ...
         'seg_nFreq', 'seg_nTime', 'orig_nFreq', 'orig_nTime');
    
    fprintf('Significant segments: %d / %d (%.1f%%)\n\n', ...
            nnz(significance_mask), numel(significance_mask), ...
            100*nnz(significance_mask)/numel(significance_mask));
end

fprintf('Nested IC segmented analysis complete.\n');

% =========================================================================
% SUPPORTING FUNCTIONS
% =========================================================================

% ========================================================================
% FUNCTION: calculate_segment_dimensions
% ------------------------------------------------------------------------
function [seg_nFreq, seg_nTime, freq_segments, time_segments] = ...
    calculate_segment_dimensions(orig_nFreq, orig_nTime, freq_factor, time_factor)
    
    % Calculate number of segments
    seg_nFreq = round(orig_nFreq * freq_factor);
    seg_nTime = round(orig_nTime * time_factor);
    
    % Ensure minimum of 1 segment
    seg_nFreq = max(1, seg_nFreq);
    seg_nTime = max(1, seg_nTime);
    
    % Create segment boundaries
    freq_edges = round(linspace(1, orig_nFreq+1, seg_nFreq+1));
    time_edges = round(linspace(1, orig_nTime+1, seg_nTime+1));
    
    % Create segment definitions
    freq_segments = cell(seg_nFreq, 1);
    time_segments = cell(seg_nTime, 1);
    
    for f = 1:seg_nFreq
        freq_segments{f} = freq_edges(f):(freq_edges(f+1)-1);
    end
    
    for t = 1:seg_nTime
        time_segments{t} = time_edges(t):(time_edges(t+1)-1);
    end
end

% ========================================================================
% FUNCTION: precompute_segmented_data_with_ic
% ------------------------------------------------------------------------
function segmented_data = precompute_segmented_data_with_ic(cluster_TF_data_warped, freq_segments, time_segments)
    
    nIC = size(cluster_TF_data_warped, 1);
    seg_nFreq = length(freq_segments);
    seg_nTime = length(time_segments);
    
    % Initialize segmented data structure
    segmented_data = cell(nIC, 1);
    
    for ic = 1:nIC
        subject_id = cluster_TF_data_warped{ic, 1};
        ic_id = cluster_TF_data_warped{ic, 2};  % Extract IC ID from column 2
        tf_struct = cluster_TF_data_warped{ic, 3};
        
        % Initialize structure for this IC
        segmented_data{ic} = struct();
        segmented_data{ic}.subject_id = subject_id;
        segmented_data{ic}.ic_id = ic_id;  % Store IC ID
        segmented_data{ic}.conditions = struct();
        
        conditions = {'P1', 'P3', 'P6'};
        condition_values = [1, 3, 6];
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            if isfield(tf_struct, cond_name) && ~isempty(tf_struct.(cond_name))
                n_trials = length(tf_struct.(cond_name));
                
                % Pre-allocate segmented trial data
                segmented_trials = zeros(seg_nFreq, seg_nTime, n_trials);
                
                for trial = 1:n_trials
                    original_data = tf_struct.(cond_name){trial}; % orig_nFreq × orig_nTime
                    
                    % Average within each segment
                    for f_seg = 1:seg_nFreq
                        for t_seg = 1:seg_nTime
                            freq_idx = freq_segments{f_seg};
                            time_idx = time_segments{t_seg};
                            
                            % Extract segment and compute mean
                            segment_data = original_data(freq_idx, time_idx);
                            segmented_trials(f_seg, t_seg, trial) = mean(segment_data(:), 'omitnan');
                        end
                    end
                end
                
                segmented_data{ic}.conditions.(cond_name) = segmented_trials;
                segmented_data{ic}.condition_values(c) = cond_val;
            end
        end
    end
end

% ========================================================================
% FUNCTION: print_nested_ic_info
% ------------------------------------------------------------------------
function print_nested_ic_info(segmented_data)
    fprintf('Nested IC Structure Summary:\n');
    
    unique_subjects = [];
    nested_ic_count = 0;
    
    for ic = 1:length(segmented_data)
        subject_id = segmented_data{ic}.subject_id;
        ic_id = segmented_data{ic}.ic_id;
        unique_subjects = [unique_subjects, subject_id]; %#ok<AGROW>
        nested_ic_count = nested_ic_count + 1;
        
        if ic <= 5  % Print first few examples
            fprintf('  IC %d: Subject %d, IC %d → Nested ID: %d_%d\n', ...
                   ic, subject_id, ic_id, subject_id, ic_id);
        end
    end
    
    unique_subjects = unique(unique_subjects);
    fprintf('  Total: %d ICs nested within %d subjects\n', nested_ic_count, length(unique_subjects));

end

% % ========================================================================
% % FUNCTION: compute_segmented_fmap_nested
% % ------------------------------------------------------------------------
% function observed_Fmap = compute_segmented_fmap_nested(segmented_data, seg_nFreq, seg_nTime)
% 
%     observed_Fmap = zeros(seg_nFreq, seg_nTime);
% 
%     parfor f = 1:seg_nFreq  % Parallel over frequency segments
%         temp_row = zeros(1, seg_nTime);
%         for t = 1:seg_nTime
%             tbl = build_segmented_datatable_nested(segmented_data, f, t);
%             if height(tbl) > 6  % Need minimum data points for LMM
%                 % Convert to categorical variables
%                 tbl.Pressure = categorical(tbl.Pressure);
%                 tbl.SubjectID = categorical(tbl.SubjectID);
%                 tbl.NestedIC = categorical(tbl.NestedIC);  % Nested IC variable
% 
%                 try
%                     % Fit LMM with nested IC structure
%                     lme = fitlme(tbl, 'Power ~ Pressure + (1|SubjectID) + (1|NestedIC)');
% 
%                     % Extract F-statistic for main effect of Pressure
%                     aov = anova(lme, 'DFMethod', 'Satterthwaite');
%                     temp_row(t) = aov.FStat(2); % Main effect of Pressure
%                 catch ME
%                     % Handle convergence failures
%                     if contains(ME.message, 'singular')
%                         % Try simpler model if convergence issues
%                         try
%                             lme = fitlme(tbl, 'Power ~ Pressure + (1|SubjectID)');
%                             aov = anova(lme, 'DFMethod', 'Satterthwaite');
%                             temp_row(t) = aov.FStat(2);
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
% % FUNCTION: build_segmented_datatable_nested
% % ------------------------------------------------------------------------
% function tbl = build_segmented_datatable_nested(segmented_data, f_seg, t_seg)
% 
%     power_vals = [];
%     press_vals = [];
%     subj_vals = [];
%     nested_ic_vals = [];  % Nested IC identifier
% 
%     nIC = length(segmented_data);
%     conditions = {'P1', 'P3', 'P6'};
%     condition_values = [1, 3, 6];
% 
%     for ic = 1:nIC
%         subject_id = segmented_data{ic}.subject_id;
%         ic_id = segmented_data{ic}.ic_id;
% 
%         % Create nested IC identifier: SubjectID_ICID
%         nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             cond_val = condition_values(c);
% 
%             if isfield(segmented_data{ic}.conditions, cond_name)
%                 trial_data = segmented_data{ic}.conditions.(cond_name);
%                 n_trials = size(trial_data, 3);
% 
%                 % Extract power values for this segment across all trials
%                 segment_powers = squeeze(trial_data(f_seg, t_seg, :));
% 
%                 % Remove NaN values
%                 valid_idx = ~isnan(segment_powers);
%                 if any(valid_idx)
%                     segment_powers = segment_powers(valid_idx);
%                     n_valid = sum(valid_idx);
% 
%                     % Add to vectors
%                     power_vals = [power_vals; segment_powers]; %#ok<AGROW>
%                     press_vals = [press_vals; repmat(cond_val, n_valid, 1)]; %#ok<AGROW>
%                     subj_vals = [subj_vals; repmat(subject_id, n_valid, 1)]; %#ok<AGROW>
%                     nested_ic_vals = [nested_ic_vals; repmat({nested_ic_id}, n_valid, 1)]; %#ok<AGROW>
%                 end
%             end
%         end
%     end
% 
%     tbl = table(power_vals, press_vals, subj_vals, nested_ic_vals, ...
%                 'VariableNames', {'Power', 'Pressure', 'SubjectID', 'NestedIC'});
% end

% % ========================================================================
% % FUNCTION: create_permuted_segmented_data_nested
% % ------------------------------------------------------------------------
% function perm_data = create_permuted_segmented_data_nested(segmented_data)
% 
%     perm_data = segmented_data; % Copy structure including IC IDs
%     nIC = length(segmented_data);
%     conditions = {'P1', 'P3', 'P6'};
% 
%     for ic = 1:nIC
%         % Get all trial data for this IC
%         all_trials = [];
%         trial_counts = [];
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             if isfield(segmented_data{ic}.conditions, cond_name)
%                 cond_data = segmented_data{ic}.conditions.(cond_name);
%                 all_trials = cat(3, all_trials, cond_data);
%                 trial_counts(c) = size(cond_data, 3); %#ok<AGROW>
%             else
%                 trial_counts(c) = 0; %#ok<AGROW>
%             end
%         end
% 
%         % Permute all trials within this IC (maintaining IC structure)
%         total_trials = size(all_trials, 3);
%         if total_trials > 0
%             perm_order = randperm(total_trials);
%             perm_trials = all_trials(:, :, perm_order);
% 
%             % Redistribute to conditions with original counts
%             start_idx = 1;
%             for c = 1:length(conditions)
%                 cond_name = conditions{c};
%                 if trial_counts(c) > 0
%                     end_idx = start_idx + trial_counts(c) - 1;
%                     perm_data{ic}.conditions.(cond_name) = perm_trials(:, :, start_idx:end_idx);
%                     start_idx = end_idx + 1;
%                 end
%             end
%         end
%     end
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
%     hVals = 0:dh:max(statMap(:));
%     TFCE = zeros(size(statMap));
% 
%     % 8-connected neighbourhood in 2D
%     conn = conndef(2,'maximal');
% 
%     for h = hVals
%         mask = statMap >= h;            % Threshold at height h
%         CC = bwconncomp(mask, conn);    % Identify clusters
% 
%         for clIdx = 1:CC.NumObjects
%             voxIdx = CC.PixelIdxList{clIdx};
%             extent = numel(voxIdx);
%             TFCE(voxIdx) = TFCE(voxIdx) + (h^H) * (extent^E) * dh;
%         end
%     end
% end
