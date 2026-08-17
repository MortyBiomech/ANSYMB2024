% ========================================================================
% FUNCTION: build_rmanova_table
% ------------------------------------------------------------------------
function tbl = build_rmanova_table(averaged_data, f, t, conditions)
    
    power_vals = [];
    pressure_vals = [];
    nested_ic_vals = {};
    
    nIC = averaged_data.nIC;
    condition_values = [1, 3, 6];
    
    for ic = 1:nIC
        nested_ic_id = averaged_data.nested_ic_ids{ic};
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            if isfield(averaged_data.data, cond_name)
                % Extract averaged power value for this IC at (f,t)
                power_val = averaged_data.data.(cond_name)(f, t, ic);
                
                % Only include non-NaN values
                if ~isnan(power_val)
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