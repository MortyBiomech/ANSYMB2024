% ========================================================================
% FUNCTION: compute_rmanova_fmap
% ------------------------------------------------------------------------
function observed_Fmap = compute_rmanova_fmap(averaged_data, nFreq, nTime)
    
    observed_Fmap = zeros(nFreq, nTime);
    conditions = {'P1', 'P3', 'P6'};
    
    parfor f = 1:nFreq  % Parallel over frequencies
        temp_row = zeros(1, nTime);
        for t = 1:nTime
            % Build data table for this (f,t) point
            tbl = build_rmanova_table(averaged_data, f, t, conditions);
            
            if height(tbl) >= 6  % Need minimum data points for ANOVA
                try
                    % Convert to categorical for proper factor treatment
                    tbl.NestedIC = categorical(tbl.NestedIC);
                    tbl.Pressure = categorical(tbl.Pressure);
                    
                    % Simple nested repeated measures ANOVA
                    % Model: Power ~ NestedIC + Pressure + NestedIC×Pressure
                    % NestedIC as random factor
                    [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
                                          'model', 'interaction', ...
                                          'random', 1, ...  % NestedIC as random factor
                                          'display', 'off');
                    
                    % Extract F-statistic for main effect of Pressure (factor 2)
                    temp_row(t) = cell2mat(T(3, 6));  % Main effect of Pressure
                    
                catch ME
                    % Handle ANOVA failures gracefully
                    temp_row(t) = 0;
                end
            end
        end
        observed_Fmap(f, :) = temp_row;
    end
end