% ========================================================================
% FUNCTION: prepare_trial_and_averaged_data
% ------------------------------------------------------------------------
function data_struct = prepare_trial_and_averaged_data(cluster_TF_data_warped, nFreq, nTime)
    
    nIC = size(cluster_TF_data_warped, 1);
    conditions = {'P1', 'P3', 'P6'};
    
    % Initialize data structure
    data_struct = struct();
    data_struct.subject_ids = [];
    data_struct.ic_ids = [];
    data_struct.nested_ic_ids = {};
    data_struct.trial_data = {};      % Store original trials for permutation
    data_struct.averaged_data = {};   % Store averaged data for analysis
    data_struct.trial_counts = {};    % Store number of trials per condition
    
    valid_ic_count = 0;
    
    for ic = 1:nIC
        subject_id = cluster_TF_data_warped{ic, 1};
        ic_id = cluster_TF_data_warped{ic, 2};
        tf_struct = cluster_TF_data_warped{ic, 3};
        
        % Create nested IC identifier
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
            data_struct.subject_ids(valid_ic_count) = subject_id;
            data_struct.ic_ids(valid_ic_count) = ic_id;
            data_struct.nested_ic_ids{valid_ic_count} = nested_ic_id;
            
            % Initialize structures for this IC
            data_struct.trial_data{valid_ic_count} = struct();
            data_struct.averaged_data{valid_ic_count} = struct();
            data_struct.trial_counts{valid_ic_count} = struct();
            
            % Process each condition
            for c = 1:length(conditions)
                cond_name = conditions{c};
                trials = tf_struct.(cond_name);
                n_trials = length(trials);
                
                % Store original trials (nFreq × nTime × nTrials)
                trial_stack = zeros(nFreq, nTime, n_trials);
                for trial = 1:n_trials
                    trial_stack(:, :, trial) = trials{trial};
                end
                data_struct.trial_data{valid_ic_count}.(cond_name) = trial_stack;
                data_struct.trial_counts{valid_ic_count}.(cond_name) = n_trials;
                
                % Compute and store averages
                averaged_trial = mean(trial_stack, 3, 'omitnan');
                data_struct.averaged_data{valid_ic_count}.(cond_name) = averaged_trial;
            end
        end
    end
    
    % Store final count
    data_struct.nIC = valid_ic_count;
    
    % Trim arrays to actual size
    data_struct.subject_ids = data_struct.subject_ids(1:valid_ic_count);
    data_struct.ic_ids = data_struct.ic_ids(1:valid_ic_count);
    data_struct.nested_ic_ids = data_struct.nested_ic_ids(1:valid_ic_count);
    data_struct.trial_data = data_struct.trial_data(1:valid_ic_count);
    data_struct.averaged_data = data_struct.averaged_data(1:valid_ic_count);
    data_struct.trial_counts = data_struct.trial_counts(1:valid_ic_count);
end