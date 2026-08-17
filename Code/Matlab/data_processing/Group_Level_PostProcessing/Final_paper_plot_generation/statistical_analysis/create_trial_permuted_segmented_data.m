% ========================================================================
% FUNCTION: create_trial_permuted_segmented_data
% ------------------------------------------------------------------------
function perm_segmented_data = create_trial_permuted_segmented_data(segmented_data)
    
    % Initialize permuted structure
    perm_segmented_data = segmented_data; % Copy metadata
    perm_segmented_data.averaged_data = cell(segmented_data.nIC, 1);
    
    conditions = {'P1', 'P3', 'P6'};
    condition_values = [1, 3, 6];
    
    % For each IC, perform trial-level shuffling on segmented data
    for ic = 1:segmented_data.nIC
        perm_segmented_data.averaged_data{ic} = struct();
        
        % Collect all segmented trials and labels
        all_segmented_trials = [];
        all_labels = [];
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            trials = segmented_data.trial_data{ic}.(cond_name);  % seg_nFreq × seg_nTime × nTrials
            n_trials = size(trials, 3);
            
            all_segmented_trials = cat(3, all_segmented_trials, trials);
            all_labels = [all_labels; repmat(cond_val, n_trials, 1)]; %#ok<AGROW>
        end
        
        % Shuffle labels
        total_trials = length(all_labels);
        shuffled_labels = all_labels(randperm(total_trials));
        
        % Recompute averages with shuffled labels
        for c = 1:length(conditions)
            cond_val = condition_values(c);
            condition_trial_indices = find(shuffled_labels == cond_val);
            
            if ~isempty(condition_trial_indices)
                condition_trials = all_segmented_trials(:, :, condition_trial_indices);
                perm_segmented_data.averaged_data{ic}.(conditions{c}) = mean(condition_trials, 3, 'omitnan');
            else
                % Handle empty case
                [seg_nFreq, seg_nTime] = size(all_segmented_trials(:, :, 1));
                perm_segmented_data.averaged_data{ic}.(conditions{c}) = nan(seg_nFreq, seg_nTime);
            end
        end
    end
end