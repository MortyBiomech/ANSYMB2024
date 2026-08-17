function [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed_psd(...
    trial_max_psd, isExperiment, P1, P3, P6, threshold_factor)
    
    % Only consider experimental trials for each pressure group
    P1_exp = intersect(P1, isExperiment);
    P3_exp = intersect(P3, isExperiment);
    P6_exp = intersect(P6, isExperiment);
    
    % Extract max PSDs for each experimental pressure group
    max_psd_P1 = trial_max_psd(P1_exp);
    max_psd_P3 = trial_max_psd(P3_exp);
    max_psd_P6 = trial_max_psd(P6_exp);
    
    % Find outliers within experimental trials only
    outliers_P1 = isoutlier(max_psd_P1, "median", "ThresholdFactor", threshold_factor);
    outliers_P3 = isoutlier(max_psd_P3, "median", "ThresholdFactor", threshold_factor);
    outliers_P6 = isoutlier(max_psd_P6, "median", "ThresholdFactor", threshold_factor);
    
    % Remove outliers from experimental trials only
    P1_clean = P1_exp(~outliers_P1);
    P3_clean = P3_exp(~outliers_P3);
    P6_clean = P6_exp(~outliers_P6);
    
    fprintf('Outliers removed (experimental trials only): P1: %d/%d, P3: %d/%d, P6: %d/%d\n', ...
        sum(outliers_P1), length(P1_exp), sum(outliers_P3), length(P3_exp), sum(outliers_P6), length(P6_exp));
end