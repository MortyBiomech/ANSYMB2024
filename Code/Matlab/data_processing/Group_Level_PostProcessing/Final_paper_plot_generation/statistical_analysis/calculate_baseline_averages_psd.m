function baseline_averages = calculate_baseline_averages_psd(psd_baseline_mean, ...
    isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode)
    % Calculate baseline averages based on the specified mode
    
    if isempty(psd_baseline_mean)
        baseline_averages = struct('P1', [], 'P3', [], 'P6', []);
        return;
    end
    
    switch baseline_mode
        case 'common'
            % Option A: One common baseline from all experimental trials (excluding outliers)
            all_clean_trials = [P1_clean, P3_clean, P6_clean];
            baseline_trials = intersect(all_clean_trials, isExperiment);
            
            if ~isempty(baseline_trials)
                common_baseline = mean(cat(2, psd_baseline_mean{baseline_trials}), 2);
                baseline_averages.P1 = common_baseline;
                baseline_averages.P3 = common_baseline;
                baseline_averages.P6 = common_baseline;
                fprintf('Using common baseline from %d trials\n', length(baseline_trials));
            else
                % Fallback to all baselines if no clean trials found
                common_baseline = mean(cat(2, psd_baseline_mean{:}), 2);
                baseline_averages.P1 = common_baseline;
                baseline_averages.P3 = common_baseline;
                baseline_averages.P6 = common_baseline;
                fprintf('Using common baseline from all %d baseline trials (fallback)\n', length(psd_baseline_mean));
            end
            
        case 'per_condition'
            % Option B: Separate baselines for each pressure condition (excluding outliers)
            P1_baseline_trials = intersect(P1_clean, isExperiment);
            P3_baseline_trials = intersect(P3_clean, isExperiment);
            P6_baseline_trials = intersect(P6_clean, isExperiment);
            
            if ~isempty(P1_baseline_trials)
                baseline_averages.P1 = mean(cat(2, psd_baseline_mean{P1_baseline_trials}), 2);
                fprintf('P1 baseline from %d trials\n', length(P1_baseline_trials));
            else
                baseline_averages.P1 = mean(cat(2, psd_baseline_mean{:}), 2);
                fprintf('P1 baseline: fallback to all trials\n');
            end
            
            if ~isempty(P3_baseline_trials)
                baseline_averages.P3 = mean(cat(2, psd_baseline_mean{P3_baseline_trials}), 2);
                fprintf('P3 baseline from %d trials\n', length(P3_baseline_trials));
            else
                baseline_averages.P3 = mean(cat(2, psd_baseline_mean{:}), 2);
                fprintf('P3 baseline: fallback to all trials\n');
            end
            
            if ~isempty(P6_baseline_trials)
                baseline_averages.P6 = mean(cat(2, psd_baseline_mean{P6_baseline_trials}), 2);
                fprintf('P6 baseline from %d trials\n', length(P6_baseline_trials));
            else
                baseline_averages.P6 = mean(cat(2, psd_baseline_mean{:}), 2);
                fprintf('P6 baseline: fallback to all trials\n');
            end
            
        otherwise
            error('Invalid baseline_mode. Use ''common'' or ''per_condition''');
    end
end