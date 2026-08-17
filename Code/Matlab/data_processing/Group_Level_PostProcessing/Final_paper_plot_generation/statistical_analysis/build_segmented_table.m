% ========================================================================
% FUNCTION: build_segmented_table
% ------------------------------------------------------------------------
function tbl = build_segmented_table(segmented_data, f_seg, t_seg, conditions)
    
    power_vals = [];
    pressure_vals = [];
    nested_ic_vals = {};
    
    nIC = segmented_data.nIC;
    condition_values = [1, 3, 6];
    
    for ic = 1:nIC
        nested_ic_id = segmented_data.nested_ic_ids{ic};
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            % Extract segmented averaged power value
            if isfield(segmented_data.averaged_data{ic}, cond_name)
                power_val = segmented_data.averaged_data{ic}.(cond_name)(f_seg, t_seg);
                
                if isfinite(power_val)
                    power_vals = [power_vals; power_val]; %#ok<AGROW>
                    pressure_vals = [pressure_vals; cond_val]; %#ok<AGROW>
                    nested_ic_vals = [nested_ic_vals; {nested_ic_id}]; %#ok<AGROW>
                end
            end
        end
    end
    
    tbl = table(power_vals, pressure_vals, nested_ic_vals, ...
                'VariableNames', {'Power', 'Pressure', 'NestedIC'});
end