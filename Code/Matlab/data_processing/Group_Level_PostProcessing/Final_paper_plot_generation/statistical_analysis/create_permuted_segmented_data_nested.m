% ========================================================================
% FUNCTION: create_permuted_segmented_data_nested
% ------------------------------------------------------------------------
function perm_data = create_permuted_segmented_data_nested(segmented_data)
    
    perm_data = segmented_data; % Copy structure including IC IDs
    nIC = length(segmented_data);
    conditions = {'P1', 'P3', 'P6'};
    
    for ic = 1:nIC
        % Get all trial data for this IC
        all_trials = [];
        trial_counts = [];
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            if isfield(segmented_data{ic}.conditions, cond_name)
                cond_data = segmented_data{ic}.conditions.(cond_name);
                all_trials = cat(3, all_trials, cond_data);
                trial_counts(c) = size(cond_data, 3); %#ok<AGROW>
            else
                trial_counts(c) = 0; %#ok<AGROW>
            end
        end
        
        % Permute all trials within this IC (maintaining IC structure)
        total_trials = size(all_trials, 3);
        if total_trials > 0
            perm_order = randperm(total_trials);
            perm_trials = all_trials(:, :, perm_order);
            
            % Redistribute to conditions with original counts
            start_idx = 1;
            for c = 1:length(conditions)
                cond_name = conditions{c};
                if trial_counts(c) > 0
                    end_idx = start_idx + trial_counts(c) - 1;
                    perm_data{ic}.conditions.(cond_name) = perm_trials(:, :, start_idx:end_idx);
                    start_idx = end_idx + 1;
                end
            end
        end
    end
end