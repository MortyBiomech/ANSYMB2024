% ========================================================================
% FUNCTION: build_segmented_datatable_nested
% ------------------------------------------------------------------------
function tbl = build_segmented_datatable_nested(segmented_data, f_seg, t_seg)
    
    power_vals = [];
    press_vals = [];
    subj_vals = [];
    nested_ic_vals = [];  % Nested IC identifier
    
    nIC = length(segmented_data);
    conditions = {'P1', 'P3', 'P6'};
    condition_values = [1, 3, 6];
    
    for ic = 1:nIC
        subject_id = segmented_data{ic}.subject_id;
        ic_id = segmented_data{ic}.ic_id;
        
        % Create nested IC identifier: SubjectID_ICID
        nested_ic_id = sprintf('%d_%d', subject_id, ic_id);
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            if isfield(segmented_data{ic}.conditions, cond_name)
                trial_data = segmented_data{ic}.conditions.(cond_name);
                n_trials = size(trial_data, 3);
                
                % Extract power values for this segment across all trials
                segment_powers = squeeze(trial_data(f_seg, t_seg, :));
                
                % Remove NaN values
                valid_idx = ~isnan(segment_powers);
                if any(valid_idx)
                    segment_powers = segment_powers(valid_idx);
                    n_valid = sum(valid_idx);
                    
                    % Add to vectors
                    power_vals = [power_vals; segment_powers]; %#ok<AGROW>
                    press_vals = [press_vals; repmat(cond_val, n_valid, 1)]; %#ok<AGROW>
                    subj_vals = [subj_vals; repmat(subject_id, n_valid, 1)]; %#ok<AGROW>
                    nested_ic_vals = [nested_ic_vals; repmat({nested_ic_id}, n_valid, 1)]; %#ok<AGROW>
                end
            end
        end
    end
    
    tbl = table(power_vals, press_vals, subj_vals, nested_ic_vals, ...
                'VariableNames', {'Power', 'Pressure', 'SubjectID', 'NestedIC'});
end