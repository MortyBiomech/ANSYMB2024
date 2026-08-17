% ========================================================================
% FUNCTION: compute_segmented_fmap_nested
% ------------------------------------------------------------------------
function observed_Fmap = compute_segmented_fmap_nested(segmented_data, seg_nFreq, seg_nTime)
    
    observed_Fmap = zeros(seg_nFreq, seg_nTime);
    
    parfor f = 1:seg_nFreq  % Parallel over frequency segments
        temp_row = zeros(1, seg_nTime);
        for t = 1:seg_nTime
            tbl = build_segmented_datatable_nested(segmented_data, f, t);
            if height(tbl) > 6  % Need minimum data points for LMM
                % Convert to categorical variables
                tbl.Pressure = categorical(tbl.Pressure);
                tbl.SubjectID = categorical(tbl.SubjectID);
                tbl.NestedIC = categorical(tbl.NestedIC);  % Nested IC variable
                
                try
                    % Fit LMM with nested IC structure
                    lme = fitlme(tbl, 'Power ~ Pressure + (1|SubjectID) + (1|NestedIC)');
                    
                    % Extract F-statistic for main effect of Pressure
                    aov = anova(lme, 'DFMethod', 'Satterthwaite');
                    temp_row(t) = aov.FStat(2); % Main effect of Pressure
                catch ME
                    % Handle convergence failures
                    if contains(ME.message, 'singular')
                        % Try simpler model if convergence issues
                        try
                            lme = fitlme(tbl, 'Power ~ Pressure + (1|SubjectID)');
                            aov = anova(lme, 'DFMethod', 'Satterthwaite');
                            temp_row(t) = aov.FStat(2);
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