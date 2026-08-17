% ========================================================================
% FUNCTION: compute_segmented_fmap
% ------------------------------------------------------------------------
function observed_Fmap = compute_segmented_fmap(segmented_data, seg_nFreq, seg_nTime)
    
    observed_Fmap = zeros(seg_nFreq, seg_nTime);
    conditions = {'P1', 'P3', 'P6'};
    
    % Process in smaller chunks if needed for memory efficiency
    for f = 1:seg_nFreq
        for t = 1:seg_nTime
            % Build table for this segmented (f,t) point
            tbl = build_segmented_table(segmented_data, f, t, conditions);
            
            if height(tbl) >= 6
                try
                    tbl.NestedIC = categorical(tbl.NestedIC);
                    tbl.Pressure = categorical(tbl.Pressure);
                    
                    [~, T, ~] = anovan(tbl.Power, {tbl.NestedIC, tbl.Pressure}, ...
                                          'model', 'interaction', ...
                                          'random', 1, ...
                                          'display', 'off');

                    observed_Fmap(f, t) = cell2mat(T(3, 6));
                    
                catch
                    observed_Fmap(f, t) = 0;
                end
            end
        end
        
        % % Progress indicator
        % if mod(f, max(1, round(seg_nFreq/10))) == 0
        %     fprintf('  Completed %d/%d frequency groups\n', f, seg_nFreq);
        % end
    end
end