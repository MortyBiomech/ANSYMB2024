% ========================================================================
% FUNCTION: create_trial_level_permuted_averages
% ------------------------------------------------------------------------
function perm_averages = create_trial_level_permuted_averages(data_struct)
    
    % Initialize permuted averages structure
    perm_averages = struct();
    perm_averages.subject_ids = data_struct.subject_ids;
    perm_averages.ic_ids = data_struct.ic_ids;
    perm_averages.nested_ic_ids = data_struct.nested_ic_ids;
    perm_averages.nIC = data_struct.nIC;
    perm_averages.averaged_data = cell(data_struct.nIC, 1);
    
    conditions = {'P1', 'P3', 'P6'};
    condition_values = [1, 3, 6];
    
    % For each IC, shuffle trial labels and recompute averages
    for ic = 1:data_struct.nIC
        % Initialize averaged data structure for this IC
        perm_averages.averaged_data{ic} = struct();
        
        % Collect all trials and their original condition labels
        all_trials = [];
        all_labels = [];
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            trials = data_struct.trial_data{ic}.(cond_name);  % nFreq × nTime × nTrials
            n_trials = size(trials, 3);
            
            % Concatenate trials along 3rd dimension
            all_trials = cat(3, all_trials, trials);
            
            % Create condition labels for these trials
            all_labels = [all_labels; repmat(cond_val, n_trials, 1)]; %#ok<AGROW>
        end
        
        % Shuffle the condition labels (permutation step)
        total_trials = length(all_labels);
        shuffled_labels = all_labels(randperm(total_trials));
        
        % Recompute averages based on shuffled labels
        for c = 1:length(conditions)
            cond_val = condition_values(c);
            
            % Find trials that now belong to this condition after shuffling
            condition_trial_indices = find(shuffled_labels == cond_val);
            
            if ~isempty(condition_trial_indices)
                % Extract trials for this condition and compute average
                condition_trials = all_trials(:, :, condition_trial_indices);
                perm_averages.averaged_data{ic}.(conditions{c}) = mean(condition_trials, 3, 'omitnan');
            else
                % Handle case where no trials assigned to this condition
                % Use NaN array as placeholder
                [nFreq, nTime] = size(all_trials(:, :, 1));
                perm_averages.averaged_data{ic}.(conditions{c}) = nan(nFreq, nTime);
            end
        end
    end
end