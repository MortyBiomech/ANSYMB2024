% ========================================================================
% FUNCTION: compute_trial_averaged_data
% ------------------------------------------------------------------------
function averaged_data = compute_trial_averaged_data(cluster_TF_data_warped, nFreq, nTime)
    
    nIC = size(cluster_TF_data_warped, 1);
    conditions = {'P1', 'P3', 'P6'};
    
    % Initialize averaged data structure
    averaged_data = struct();
    averaged_data.subject_ids = [];
    averaged_data.ic_ids = [];
    averaged_data.nested_ic_ids = {};  % Nested IC identifiers (SubjectID_ICID)
    averaged_data.data = struct();
    
    % Initialize data arrays for each condition
    for c = 1:length(conditions)
        cond_name = conditions{c};
        averaged_data.data.(cond_name) = [];
    end
    
    valid_ic_count = 0;
    
    for ic = 1:nIC
        subject_id = cluster_TF_data_warped{ic, 1};
        ic_id = cluster_TF_data_warped{ic, 2};
        tf_struct = cluster_TF_data_warped{ic, 3};
        
        % Create nested IC identifier: SubjectID_ICID
        nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
        
        % Check if this IC has data for all conditions
        has_all_conditions = true;
        for c = 1:length(conditions)
            cond_name = conditions{c};
            if ~isfield(tf_struct, cond_name) || isempty(tf_struct.(cond_name))
                has_all_conditions = false;
                break;
            end
        end
        
        if has_all_conditions
            valid_ic_count = valid_ic_count + 1;
            averaged_data.subject_ids(valid_ic_count) = subject_id;
            averaged_data.ic_ids(valid_ic_count) = ic_id;
            averaged_data.nested_ic_ids{valid_ic_count} = nested_ic_id;
            
            % Average trials for each condition
            for c = 1:length(conditions)
                cond_name = conditions{c};
                trials = tf_struct.(cond_name);
                n_trials = length(trials);
                
                % Stack all trials for this condition
                trial_stack = zeros(nFreq, nTime, n_trials);
                for trial = 1:n_trials
                    trial_stack(:, :, trial) = trials{trial};
                end
                
                % Compute mean across trials
                averaged_trial = mean(trial_stack, 3, 'omitnan');
                
                % Store in data structure
                if isempty(averaged_data.data.(cond_name))
                    averaged_data.data.(cond_name) = zeros(nFreq, nTime, nIC);
                end
                averaged_data.data.(cond_name)(:, :, valid_ic_count) = averaged_trial;
            end
        end
    end
    
    % Trim data to actual number of valid ICs
    averaged_data.subject_ids = averaged_data.subject_ids(1:valid_ic_count);
    averaged_data.ic_ids = averaged_data.ic_ids(1:valid_ic_count);
    averaged_data.nested_ic_ids = averaged_data.nested_ic_ids(1:valid_ic_count);
    averaged_data.nIC = valid_ic_count;
    
    for c = 1:length(conditions)
        cond_name = conditions{c};
        averaged_data.data.(cond_name) = averaged_data.data.(cond_name)(:, :, 1:valid_ic_count);
    end
end