% ========================================================================
% FUNCTION: compute_fmap_from_averages
% ------------------------------------------------------------------------
function observed_Fmap = compute_fmap_from_averages(data_struct, nFreq, nTime)
    
    observed_Fmap = zeros(nFreq, nTime);
    conditions = {'P1', 'P3', 'P6'};
    
    parfor f = 1:nFreq  % Parallel over frequencies
        temp_row = zeros(1, nTime);
        for t = 1:nTime
            % Build data table for this (f,t) point using averaged data
            tbl = build_table_from_averages(data_struct, f, t, conditions);
            
            if height(tbl) >= 6  % Need minimum data points for ANOVA
                try
                    % Convert to categorical variables
                    tbl.NestedIC = categorical(tbl.NestedIC);
                    tbl.Pressure = categorical(tbl.Pressure);
                    
                    % Perform simple nested repeated measures ANOVA
                    [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
                                          'model', 'interaction', ...
                                          'random', 1, ...  % NestedIC as random factor
                                          'display', 'off');

                    temp_row(t) = cell2mat(T(3, 6));
                    
                catch ME
                    % Handle ANOVA failures
                    if contains(ME.message, 'singular') || contains(ME.message, 'rank deficient')
                        try
                            % Try simpler linear model if interaction fails
                            [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
                                                  'model', 'linear', ...
                                                  'random', 1, ...
                                                  'display', 'off');
                            temp_row(t) = cell2mat(T(3, 6));
                        catch
                            temp_row(t) = 0;
                        end
                    else
                        temp_row(t) = 0;
                    end
                end
            end
        end
        observed_Fmap(f, :) = temp_row;
    end
end