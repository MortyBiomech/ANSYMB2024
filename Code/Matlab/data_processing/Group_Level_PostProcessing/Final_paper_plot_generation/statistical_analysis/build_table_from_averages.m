% ========================================================================
% FUNCTION: build_table_from_averages
% ------------------------------------------------------------------------
function tbl = build_table_from_averages(data_struct, f, t, conditions)
    
    power_vals = [];
    pressure_vals = [];
    nested_ic_vals = {};
    
    nIC = data_struct.nIC;
    condition_values = [1, 3, 6];
    
    for ic = 1:nIC
        nested_ic_id = data_struct.nested_ic_ids{ic};
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            cond_val = condition_values(c);
            
            % Extract averaged power value for this IC at (f,t)
            if isfield(data_struct.averaged_data{ic}, cond_name)
                power_val = data_struct.averaged_data{ic}.(cond_name)(f, t);
                
                % Only include finite values
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