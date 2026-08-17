% -------------------------------------------------------------------------
% SEGMENTED TRIAL-LEVEL SHUFFLING NESTED REPEATED MEASURES ANOVA
% -------------------------------------------------------------------------
% This script performs cluster-based permutation testing with dramatic
% computational reduction through time-frequency segmentation:
% - Time: divided into 5% segments (20 segments total)  
% - Frequency: grouped every 3 frequencies and averaged
% - Trial-level shuffling for statistical validity
% - Simple nested ANOVA with TFCE enhancement
%
% Expected speedup: ~300x reduction in analysis points
% Requirements: Statistics and Machine Learning Toolbox
% -------------------------------------------------------------------------

% ANALYSIS PARAMETERS
p_value = 0.01;                 % Significance level
num_permutations = 1000;        % Number of permutations
tfce_H = 2;                     % TFCE height exponent
tfce_E = 0.5;                   % TFCE extent exponent
dh = 0.1;                       % TFCE integration step

% SEGMENTATION PARAMETERS
time_segment_percent = 0.02;    % Use 5% time segments (20 segments total)
freq_group_size = 2;            % Group every 3 frequencies together

fprintf('Starting segmented trial-level shuffling analysis...\n');
fprintf('Time segmentation: %.0f%% segments, Frequency grouping: every %d frequencies\n', ...
        time_segment_percent*100, freq_group_size);
fprintf('Model: Power ~ NestedIC + Pressure + NestedIC×Pressure (NestedIC as random)\n');
fprintf('Permutations: %d, Alpha: %.3f\n\n', num_permutations, p_value);

% -------------------------------------------------------------------------
% MAIN LOOP OVER 9 CLUSTERS
% -------------------------------------------------------------------------
for clusterIdx = 1:9
    fprintf('==== CLUSTER %d (Segmented Analysis) ====================\n', clusterIdx);
    
    % Load data for this cluster
    load(sprintf('cluster%d_TF_data_warped.mat', clusterIdx), 'cluster_TF_data_warped');
    
    % Get original dimensions
    [orig_nFreq, orig_nTime] = size(cluster_TF_data_warped{1,3}.P1{1});
    nIC = size(cluster_TF_data_warped, 1);
    
    % Calculate segmentation parameters
    [seg_nFreq, seg_nTime, freq_groups, time_segments] = ...
        calculate_segmentation_params(orig_nFreq, orig_nTime, freq_group_size, time_segment_percent);
    
    fprintf('Original: %d freq × %d time = %d points\n', orig_nFreq, orig_nTime, orig_nFreq*orig_nTime);
    fprintf('Segmented: %d freq groups × %d time segments = %d points\n', ...
            seg_nFreq, seg_nTime, seg_nFreq*seg_nTime);
    fprintf('Computational reduction: %.0fx faster\n\n', (orig_nFreq*orig_nTime)/(seg_nFreq*seg_nTime));
    
    % Prepare segmented data structure with trials and averages
    fprintf('Preparing segmented data with trial-level information...\n');
    segmented_data = prepare_segmented_trial_data(cluster_TF_data_warped, freq_groups, time_segments);
    
    % % Print structure info
    % print_segmented_structure_info(segmented_data, seg_nFreq, seg_nTime);
    
    % Compute observed F-map on segmented data
    fprintf('Computing observed F-statistics on segmented data...\n');
    tic
    observed_Fmap = compute_segmented_fmap(segmented_data, seg_nFreq, seg_nTime);
    timeC = toc;

    % Apply TFCE to observed F-map
    fprintf('Applying TFCE to observed F-map...\n');
    observed_TFCE = tfce2d(observed_Fmap, tfce_H, tfce_E, dh);
    
    % Permutation testing with trial-level shuffling
    fprintf('Running %d permutations with trial-level shuffling...\n', num_permutations);
    max_TFCE_null = zeros(num_permutations, 1);
    
    % Use smaller parallel pool for memory efficiency
    current_pool = gcp('nocreate');
    if isempty(current_pool)
        parpool('threads', min(4, feature('numcores')));
    end
    
    parfor perm = 1:num_permutations
        % Create trial-level permuted segmented averages
        perm_seg_averages = create_trial_permuted_segmented_data(segmented_data);
        
        % Compute F-map for permuted data
        perm_Fmap = compute_segmented_fmap(perm_seg_averages, seg_nFreq, seg_nTime);
        
        % Apply TFCE
        perm_TFCE = tfce2d(perm_Fmap, tfce_H, tfce_E, dh);
        
        % Store maximum TFCE value
        max_TFCE_null(perm) = max(perm_TFCE(:));
    end
    
    % Determine significance threshold
    tfce_threshold = prctile(max_TFCE_null, 100*(1-p_value));
    
    % Generate final significance mask
    significance_mask = observed_TFCE > tfce_threshold;


    % Keep only clusters (8-connected) that contain ≥8 voxels
    conn = conndef(2,'maximal');        % 8-connected neighbourhood in 2-D
    CC   = bwconncomp(significance_mask, conn);
    
    final_mask = false(size(significance_mask));
    cluster_sizes = cellfun(@numel, CC.PixelIdxList);
    
    for cl = 1:CC.NumObjects
        if cluster_sizes(cl) >= 8          % cluster size threshold
            final_mask(CC.PixelIdxList{cl}) = true;
        end
    end


    
    % Save results with segmentation information
    save(sprintf('cluster%d_segmented_trial_results.mat', clusterIdx), ...
         'observed_Fmap', 'observed_TFCE', 'significance_mask', ...
         'tfce_threshold', 'max_TFCE_null', 'seg_nFreq', 'seg_nTime', ...
         'freq_groups', 'time_segments', 'orig_nFreq', 'orig_nTime');
    
    fprintf('Significant segments: %d / %d (%.1f%%)\n\n', ...
            nnz(significance_mask), numel(significance_mask), ...
            100*nnz(significance_mask)/numel(significance_mask));
end

fprintf('Segmented trial-level analysis complete.\n');

% =========================================================================
% SUPPORTING FUNCTIONS
% =========================================================================

% % ========================================================================
% % FUNCTION: calculate_segmentation_params
% % ------------------------------------------------------------------------
% function [seg_nFreq, seg_nTime, freq_groups, time_segments] = ...
%     calculate_segmentation_params(orig_nFreq, orig_nTime, freq_group_size, time_segment_percent)
% 
%     % Calculate number of segments
%     seg_nTime = round(1 / time_segment_percent);  % 5% = 20 segments
%     seg_nFreq = ceil(orig_nFreq / freq_group_size);  % Group every 3 frequencies
% 
%     % Create time segments (5% each)
%     time_edges = round(linspace(1, orig_nTime+1, seg_nTime+1));
%     time_segments = cell(seg_nTime, 1);
%     for t = 1:seg_nTime
%         time_segments{t} = time_edges(t):(time_edges(t+1)-1);
%     end
% 
%     % Create frequency groups (every 3 frequencies)
%     freq_groups = cell(seg_nFreq, 1);
%     for f = 1:seg_nFreq
%         start_freq = (f-1) * freq_group_size + 1;
%         end_freq = min(f * freq_group_size, orig_nFreq);
%         freq_groups{f} = start_freq:end_freq;
%     end
% end

% % ========================================================================
% % FUNCTION: prepare_segmented_trial_data
% % ------------------------------------------------------------------------
% function segmented_data = prepare_segmented_trial_data(cluster_TF_data_warped, freq_groups, time_segments)
% 
%     nIC = size(cluster_TF_data_warped, 1);
%     conditions = {'P1', 'P3', 'P6'};
%     seg_nFreq = length(freq_groups);
%     seg_nTime = length(time_segments);
% 
%     % Initialize segmented data structure
%     segmented_data = struct();
%     segmented_data.subject_ids = [];
%     segmented_data.ic_ids = [];
%     segmented_data.nested_ic_ids = {};
%     segmented_data.trial_data = {};      % Segmented trial data
%     segmented_data.averaged_data = {};   % Segmented averaged data
%     segmented_data.trial_counts = {};
% 
%     valid_ic_count = 0;
% 
%     for ic = 1:nIC
%         subject_id = cluster_TF_data_warped{ic, 1};
%         ic_id = cluster_TF_data_warped{ic, 2};
%         tf_struct = cluster_TF_data_warped{ic, 3};
% 
%         nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
% 
%         % Check if IC has all conditions
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
%             segmented_data.subject_ids(valid_ic_count) = subject_id;
%             segmented_data.ic_ids(valid_ic_count) = ic_id;
%             segmented_data.nested_ic_ids{valid_ic_count} = nested_ic_id;
% 
%             % Initialize structures for this IC
%             segmented_data.trial_data{valid_ic_count} = struct();
%             segmented_data.averaged_data{valid_ic_count} = struct();
%             segmented_data.trial_counts{valid_ic_count} = struct();
% 
%             % Process each condition
%             for c = 1:length(conditions)
%                 cond_name = conditions{c};
%                 trials = tf_struct.(cond_name);
%                 n_trials = length(trials);
% 
%                 % Create segmented trial data (seg_nFreq × seg_nTime × n_trials)
%                 segmented_trials = zeros(seg_nFreq, seg_nTime, n_trials);
% 
%                 for trial = 1:n_trials
%                     original_trial = trials{trial};  % orig_nFreq × orig_nTime
% 
%                     % Segment this trial
%                     for f_seg = 1:seg_nFreq
%                         for t_seg = 1:seg_nTime
%                             freq_indices = freq_groups{f_seg};
%                             time_indices = time_segments{t_seg};
% 
%                             % Extract segment and compute average
%                             segment_data = original_trial(freq_indices, time_indices);
%                             segmented_trials(f_seg, t_seg, trial) = mean(segment_data(:), 'omitnan');
%                         end
%                     end
%                 end
% 
%                 % Store segmented trial data
%                 segmented_data.trial_data{valid_ic_count}.(cond_name) = segmented_trials;
%                 segmented_data.trial_counts{valid_ic_count}.(cond_name) = n_trials;
% 
%                 % Compute and store segmented averages
%                 segmented_average = mean(segmented_trials, 3, 'omitnan');
%                 segmented_data.averaged_data{valid_ic_count}.(cond_name) = segmented_average;
%             end
%         end
%     end
% 
%     % Store final count and trim arrays
%     segmented_data.nIC = valid_ic_count;
%     segmented_data.subject_ids = segmented_data.subject_ids(1:valid_ic_count);
%     segmented_data.ic_ids = segmented_data.ic_ids(1:valid_ic_count);
%     segmented_data.nested_ic_ids = segmented_data.nested_ic_ids(1:valid_ic_count);
%     segmented_data.trial_data = segmented_data.trial_data(1:valid_ic_count);
%     segmented_data.averaged_data = segmented_data.averaged_data(1:valid_ic_count);
%     segmented_data.trial_counts = segmented_data.trial_counts(1:valid_ic_count);
% end

% ========================================================================
% FUNCTION: print_segmented_structure_info
% ------------------------------------------------------------------------
function print_segmented_structure_info(segmented_data, seg_nFreq, seg_nTime)
    fprintf('Segmented trial-level data structure:\n');
    fprintf('  Total valid ICs: %d\n', segmented_data.nIC);
    fprintf('  Segmented dimensions: %d freq groups × %d time segments\n', seg_nFreq, seg_nTime);
    
    unique_subjects = unique(segmented_data.subject_ids);
    fprintf('  Unique subjects: %d\n', length(unique_subjects));
    
    conditions = {'P1', 'P3', 'P6'};
    
    % Show example trial counts for first few ICs
    for ic = 1:min(3, segmented_data.nIC)
        fprintf('  IC %s: ', segmented_data.nested_ic_ids{ic});
        for c = 1:length(conditions)
            cond_name = conditions{c};
            n_trials = segmented_data.trial_counts{ic}.(cond_name);
            fprintf('%s=%d ', cond_name, n_trials);
        end
        fprintf('trials\n');
    end
    
    fprintf('  Average ICs per subject: %.1f\n\n', segmented_data.nIC / length(unique_subjects));
end

% % ========================================================================
% % FUNCTION: compute_segmented_fmap
% % ------------------------------------------------------------------------
% function observed_Fmap = compute_segmented_fmap(segmented_data, seg_nFreq, seg_nTime)
% 
%     observed_Fmap = zeros(seg_nFreq, seg_nTime);
%     conditions = {'P1', 'P3', 'P6'};
% 
%     % Process in smaller chunks if needed for memory efficiency
%     for f = 1:seg_nFreq
%         for t = 1:seg_nTime
%             % Build table for this segmented (f,t) point
%             tbl = build_segmented_table(segmented_data, f, t, conditions);
% 
%             if height(tbl) >= 6
%                 try
%                     tbl.NestedIC = categorical(tbl.NestedIC);
%                     tbl.Pressure = categorical(tbl.Pressure);
% 
%                     [~, ~, stats] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
%                                           'model', 'interaction', ...
%                                           'random', 1, ...
%                                           'display', 'off');
% 
%                     if length(stats.fstat) >= 2
%                         observed_Fmap(f, t) = stats.fstat(2);
%                     end
% 
%                 catch
%                     observed_Fmap(f, t) = 0;
%                 end
%             end
%         end
% 
%         % Progress indicator
%         if mod(f, max(1, round(seg_nFreq/10))) == 0
%             fprintf('  Completed %d/%d frequency groups\n', f, seg_nFreq);
%         end
%     end
% end

% % ========================================================================
% % FUNCTION: build_segmented_table
% % ------------------------------------------------------------------------
% function tbl = build_segmented_table(segmented_data, f_seg, t_seg, conditions)
% 
%     power_vals = [];
%     pressure_vals = [];
%     nested_ic_vals = {};
% 
%     nIC = segmented_data.nIC;
%     condition_values = [1, 3, 6];
% 
%     for ic = 1:nIC
%         nested_ic_id = segmented_data.nested_ic_ids{ic};
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             cond_val = condition_values(c);
% 
%             % Extract segmented averaged power value
%             if isfield(segmented_data.averaged_data{ic}, cond_name)
%                 power_val = segmented_data.averaged_data{ic}.(cond_name)(f_seg, t_seg);
% 
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
% % FUNCTION: create_trial_permuted_segmented_data
% % ------------------------------------------------------------------------
% function perm_segmented_data = create_trial_permuted_segmented_data(segmented_data)
% 
%     % Initialize permuted structure
%     perm_segmented_data = segmented_data; % Copy metadata
%     perm_segmented_data.averaged_data = cell(segmented_data.nIC, 1);
% 
%     conditions = {'P1', 'P3', 'P6'};
%     condition_values = [1, 3, 6];
% 
%     % For each IC, perform trial-level shuffling on segmented data
%     for ic = 1:segmented_data.nIC
%         perm_segmented_data.averaged_data{ic} = struct();
% 
%         % Collect all segmented trials and labels
%         all_segmented_trials = [];
%         all_labels = [];
% 
%         for c = 1:length(conditions)
%             cond_name = conditions{c};
%             cond_val = condition_values(c);
% 
%             trials = segmented_data.trial_data{ic}.(cond_name);  % seg_nFreq × seg_nTime × nTrials
%             n_trials = size(trials, 3);
% 
%             all_segmented_trials = cat(3, all_segmented_trials, trials);
%             all_labels = [all_labels; repmat(cond_val, n_trials, 1)]; %#ok<AGROW>
%         end
% 
%         % Shuffle labels
%         total_trials = length(all_labels);
%         shuffled_labels = all_labels(randperm(total_trials));
% 
%         % Recompute averages with shuffled labels
%         for c = 1:length(conditions)
%             cond_val = condition_values(c);
%             condition_trial_indices = find(shuffled_labels == cond_val);
% 
%             if ~isempty(condition_trial_indices)
%                 condition_trials = all_segmented_trials(:, :, condition_trial_indices);
%                 perm_segmented_data.averaged_data{ic}.(conditions{c}) = mean(condition_trials, 3, 'omitnan');
%             else
%                 % Handle empty case
%                 [seg_nFreq, seg_nTime] = size(all_segmented_trials(:, :, 1));
%                 perm_segmented_data.averaged_data{ic}.(conditions{c}) = nan(seg_nFreq, seg_nTime);
%             end
%         end
%     end
% end

% ========================================================================
% FUNCTION: tfce2d - Threshold-Free Cluster Enhancement
% ------------------------------------------------------------------------
function TFCE = tfce2d(statMap, H, E, dh)
    if nargin<4, dh = 0.1; end
    if nargin<3, E = 0.5; end
    if nargin<2, H = 2; end
    
    statMap(statMap<0) = 0;
    statMap(~isfinite(statMap)) = 0;
    
    if max(statMap(:)) == 0
        TFCE = zeros(size(statMap));
        return;
    end
    
    hVals = 0:dh:max(statMap(:));
    TFCE = zeros(size(statMap));
    conn = conndef(2,'maximal');
    
    for h = hVals
        if h == 0, continue; end
        mask = statMap >= h;
        if ~any(mask(:)), continue; end
        
        CC = bwconncomp(mask, conn);
        for clIdx = 1:CC.NumObjects
            voxIdx = CC.PixelIdxList{clIdx};
            extent = numel(voxIdx);
            if extent > 0
                TFCE(voxIdx) = TFCE(voxIdx) + (h^H) * (extent^E) * dh;
            end
        end
    end
end
