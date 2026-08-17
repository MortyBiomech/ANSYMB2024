% ========================================================================
% FUNCTION: prepare_segmented_trial_data
% ------------------------------------------------------------------------
function segmented_data = prepare_segmented_trial_data(cluster_TF_data_warped, freq_groups, time_segments)
    
    nIC = size(cluster_TF_data_warped, 1);
    conditions = {'P1', 'P3', 'P6'};
    seg_nFreq = length(freq_groups);
    seg_nTime = length(time_segments);
    
    % Initialize segmented data structure
    segmented_data = struct();
    segmented_data.subject_ids = [];
    segmented_data.ic_ids = [];
    segmented_data.nested_ic_ids = {};
    segmented_data.trial_data = {};      % Segmented trial data
    segmented_data.averaged_data = {};   % Segmented averaged data
    segmented_data.trial_counts = {};
    
    valid_ic_count = 0;
    
    for ic = 1:nIC
        subject_id = cluster_TF_data_warped{ic, 1};
        ic_id = cluster_TF_data_warped{ic, 2};
        tf_struct = cluster_TF_data_warped{ic, 3};
        
        nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
        
        % Check if IC has all conditions
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
            segmented_data.subject_ids(valid_ic_count) = subject_id;
            segmented_data.ic_ids(valid_ic_count) = ic_id;
            segmented_data.nested_ic_ids{valid_ic_count} = nested_ic_id;
            
            % Initialize structures for this IC
            segmented_data.trial_data{valid_ic_count} = struct();
            segmented_data.averaged_data{valid_ic_count} = struct();
            segmented_data.trial_counts{valid_ic_count} = struct();
            
            % Process each condition
            for c = 1:length(conditions)
                cond_name = conditions{c};
                trials = tf_struct.(cond_name);
                n_trials = length(trials);
                
                % Create segmented trial data (seg_nFreq × seg_nTime × n_trials)
                segmented_trials = zeros(seg_nFreq, seg_nTime, n_trials);
                
                for trial = 1:n_trials
                    original_trial = trials{trial};  % orig_nFreq × orig_nTime
                    
                    % Segment this trial
                    for f_seg = 1:seg_nFreq
                        for t_seg = 1:seg_nTime
                            freq_indices = freq_groups{f_seg};
                            time_indices = time_segments{t_seg};
                            
                            % Extract segment and compute average
                            segment_data = original_trial(freq_indices, time_indices);
                            segmented_trials(f_seg, t_seg, trial) = mean(segment_data(:), 'omitnan');
                        end
                    end
                end
                
                % Store segmented trial data
                segmented_data.trial_data{valid_ic_count}.(cond_name) = segmented_trials;
                segmented_data.trial_counts{valid_ic_count}.(cond_name) = n_trials;
                
                % Compute and store segmented averages
                segmented_average = mean(segmented_trials, 3, 'omitnan');
                segmented_data.averaged_data{valid_ic_count}.(cond_name) = segmented_average;
            end
        end
    end
    
    % Store final count and trim arrays
    segmented_data.nIC = valid_ic_count;
    segmented_data.subject_ids = segmented_data.subject_ids(1:valid_ic_count);
    segmented_data.ic_ids = segmented_data.ic_ids(1:valid_ic_count);
    segmented_data.nested_ic_ids = segmented_data.nested_ic_ids(1:valid_ic_count);
    segmented_data.trial_data = segmented_data.trial_data(1:valid_ic_count);
    segmented_data.averaged_data = segmented_data.averaged_data(1:valid_ic_count);
    segmented_data.trial_counts = segmented_data.trial_counts(1:valid_ic_count);
end