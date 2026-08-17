function [TF_struct, freqs, freqs_baseline] = main_TF_calculation_per_IC(subject, ic, trial_info, subject_data, EEGLAB_source, EEGLAB_events)
    % Highly optimized time-frequency calculation per independent component
    % Key optimization: CWT computed only ONCE per trial
    
    %% Configuration parameters
    SAMPLING_RATE = 500;
    FREQ_LIMITS = [1 70];
    VOICES_PER_OCTAVE = 32;
    BASELINE_WINDOW = 1000; % samples before event
    OUTLIER_THRESHOLD = 3;
    
    % Baseline correction options:
    % 'common' - One baseline average from all experimental trials (excluding outliers)
    % 'per_condition' - Separate baseline averages for P1, P3, P6 (excluding outliers)
    BASELINE_MODE = 'common';  % Switch between 'common' and 'per_condition'
    
    %% Baseline calculations
    [power_baseline_mean, freqs_baseline] = calculate_baseline_power(...
        subject, EEGLAB_events, EEGLAB_source, SAMPLING_RATE, FREQ_LIMITS, ...
        VOICES_PER_OCTAVE, BASELINE_WINDOW);
    
    %% Get experimental trials and pressure groups
    [P1, P3, P6, isExperiment] = get_pressure_groups(subject, trial_info);
    
    %% Single-pass processing: CWT + outlier detection + grand means + final processing
    [TF_struct, freqs] = process_all_trials_single_pass(subject, ic, trial_info, ...
        subject_data, isExperiment, P1, P3, P6, power_baseline_mean, freqs_baseline, ...
        SAMPLING_RATE, FREQ_LIMITS, VOICES_PER_OCTAVE, OUTLIER_THRESHOLD, BASELINE_MODE);
    
end

%% Helper Functions

function [power_baseline_mean, freqs_baseline] = calculate_baseline_power(...
    subject, EEGLAB_events, EEGLAB_source, sampling_rate, freq_limits, ...
    voices_per_octave, baseline_window)
    
    if subject >= 10
        event_type = 'SM_Start_Move';
    else
        event_type = 'SB_Start_Beep';
    end
    
    % Find event indices
    event_ids = strcmp({EEGLAB_events.type}, event_type);
    latencies = [EEGLAB_events(event_ids).latency];
    
    % Extract baseline segments
    n_events = length(latencies);
    baselines = cell(1, n_events);
    
    for i = 1:n_events
        start_idx = latencies(i) - baseline_window;
        end_idx = latencies(i) - 1;
        baselines{i} = EEGLAB_source(start_idx:end_idx);
    end
    
    % Calculate CWT for all baselines
    [cwt_coeffs, freqs_baseline] = cellfun(@(x) cwt(x, 'amor', sampling_rate, ...
        'FrequencyLimits', freq_limits, 'VoicesPerOctave', voices_per_octave), ...
        baselines, 'UniformOutput', false);
    
    % Calculate power and mean
    power_baseline = cellfun(@(x) abs(x).^2, cwt_coeffs, 'UniformOutput', false);
    power_baseline_mean = cellfun(@(x) mean(x, 2), power_baseline, 'UniformOutput', false);
    
    % Use first frequency vector (all should be identical)
    freqs_baseline = freqs_baseline{1};
end

function [P1, P3, P6, isExperiment] = get_pressure_groups(subject, trial_info)
    % Get experimental trials
    if subject < 10
        isExperiment = 1:length(trial_info);
    else
        isExperiment = find(cellfun(@(x) strcmp(x.General.Description, 'Experiment'), trial_info));
    end
    
    % Extract pressure levels more efficiently
    pressures = cellfun(@(x) x.General.Pressure, trial_info);
    
    P1 = intersect(find(pressures == 1), isExperiment);
    P3 = intersect(find(pressures == 3), isExperiment);
    P6 = intersect(find(pressures == 6), isExperiment);
end

function [TF_struct, freqs] = process_all_trials_single_pass(subject, ic, trial_info, ...
    subject_data, isExperiment, P1, P3, P6, power_baseline_mean, freqs_baseline, ...
    sampling_rate, freq_limits, voices_per_octave, outlier_threshold, baseline_mode)
    
    % Pre-allocate storage for all trials
    n_trials = length(subject_data);
    trial_powers = cell(n_trials, 1);           % Raw power data
    trial_times = cell(n_trials, 1);            % Time vectors
    trial_max_power = zeros(n_trials, 1);       % For outlier detection
    % trial_warped_data = cell(n_trials, 1);      % Final processed data
    freqs = [];
    
    %% SINGLE PASS: Compute CWT once per trial and extract all needed info
    fprintf('Processing %d out of %d trials with single CWT calculation per trial...\n', ...
        sum(ismember(1:n_trials, isExperiment)), n_trials);
    
    for trial = 1:n_trials
        if ~ismember(trial, isExperiment)
            continue;
        end
        
        % Extract signal and times
        [signal, times] = extract_trial_data(subject, subject_data{trial}, ic);
        
        % *** SINGLE CWT COMPUTATION PER TRIAL ***
        [cwt_coeffs, freqs] = cwt(signal, 'amor', sampling_rate, ...
            'FrequencyLimits', freq_limits, 'VoicesPerOctave', voices_per_octave);
        power = abs(cwt_coeffs).^2;
        
        % Store data for later processing
        trial_powers{trial} = power;
        trial_times{trial} = times;
        trial_max_power(trial) = max(power(:));
        
        % if mod(trial, 10) == 0
        %     fprintf('Processed %d trials\n', trial);
        % end
    end
    
    %% Outlier detection using pre-computed max powers
    [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed(...
        trial_max_power, isExperiment, P1, P3, P6, outlier_threshold);
    
    % %% Optional: Calculate and plot grand means
    % if nargout > 0  % Only if outputs requested
    %     plot_grand_means_from_precomputed(trial_powers, freqs, isExperiment, ...
    %         P1_clean, P3_clean, P6_clean, power_baseline_mean, freqs_baseline, subject, ic);
    % end
    
    %% Final processing: Time warping and baseline correction
    baseline_averages = calculate_baseline_averages(power_baseline_mean, ...
        isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode);
    
    % Pre-allocate final storage
    TF_P1 = cell(3, length(P1_clean));
    TF_P3 = cell(3, length(P3_clean));
    TF_P6 = cell(3, length(P6_clean));
    
    P1_idx = 1; P3_idx = 1; P6_idx = 1;
    
    fprintf('Applying time warping and baseline correction (mode: %s)...\n', baseline_mode);
    
    for trial = 1:n_trials
        if ~ismember(trial, isExperiment)
            continue;
        end
        
        % Skip if this trial was marked as outlier
        if ~(ismember(trial, P1_clean) || ismember(trial, P3_clean) || ismember(trial, P6_clean))
            continue;
        end
        
        % Use pre-computed power data
        power = trial_powers{trial};
        times = trial_times{trial};
        
        % Time warping
        [power_warped, median_flx, median_ext] = perform_time_warping(...
            power, times, trial_info{trial});
        
        % Baseline correction
        if ismember(trial, P1_clean)
            power_warped = apply_baseline_correction(power_warped, ...
                baseline_averages.P1, freqs_baseline);
        elseif ismember(trial, P3_clean)
            power_warped = apply_baseline_correction(power_warped, ...
                baseline_averages.P3, freqs_baseline);
        elseif ismember(trial, P6_clean)
            power_warped = apply_baseline_correction(power_warped, ...
                baseline_averages.P6, freqs_baseline);
        end
        
        % Average across epochs
        trial_TF = mean(cat(3, power_warped{:}), 3);
        
        % Store in appropriate pressure group
        if ismember(trial, P1_clean)
            TF_P1{1, P1_idx} = trial_TF;
            TF_P1{2, P1_idx} = median_flx;
            TF_P1{3, P1_idx} = median_ext;
            P1_idx = P1_idx + 1;
        elseif ismember(trial, P3_clean)
            TF_P3{1, P3_idx} = trial_TF;
            TF_P3{2, P3_idx} = median_flx;
            TF_P3{3, P3_idx} = median_ext;
            P3_idx = P3_idx + 1;
        elseif ismember(trial, P6_clean)
            TF_P6{1, P6_idx} = trial_TF;
            TF_P6{2, P6_idx} = median_flx;
            TF_P6{3, P6_idx} = median_ext;
            P6_idx = P6_idx + 1;
        end
    end
    
    % Remove empty columns
    TF_P1 = remove_empty_columns(TF_P1);
    TF_P3 = remove_empty_columns(TF_P3);
    TF_P6 = remove_empty_columns(TF_P6);
    
    % Create output structure
    TF_struct = struct('P1', {TF_P1}, 'P3', {TF_P3}, 'P6', {TF_P6});
    
    fprintf('Processing complete. P1: %d trials, P3: %d trials, P6: %d trials\n', ...
        size(TF_P1, 2), size(TF_P3, 2), size(TF_P6, 2));
end

function baseline_averages = calculate_baseline_averages(power_baseline_mean, ...
    isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode)
    % Calculate baseline averages based on the specified mode
    
    if isempty(power_baseline_mean)
        baseline_averages = struct('P1', [], 'P3', [], 'P6', []);
        return;
    end
    
    switch baseline_mode
        case 'common'
            % Option A: One common baseline from all experimental trials (excluding outliers)
            all_clean_trials = [P1_clean, P3_clean, P6_clean];
            baseline_trials = intersect(all_clean_trials, isExperiment);
            
            if ~isempty(baseline_trials)
                common_baseline = mean(cat(2, power_baseline_mean{baseline_trials}), 2);
                baseline_averages.P1 = common_baseline;
                baseline_averages.P3 = common_baseline;
                baseline_averages.P6 = common_baseline;
                fprintf('Using common baseline from %d trials\n', length(baseline_trials));
            else
                % Fallback to all baselines if no clean trials found
                common_baseline = mean(cat(2, power_baseline_mean{:}), 2);
                baseline_averages.P1 = common_baseline;
                baseline_averages.P3 = common_baseline;
                baseline_averages.P6 = common_baseline;
                fprintf('Using common baseline from all %d baseline trials (fallback)\n', length(power_baseline_mean));
            end
            
        case 'per_condition'
            % Option B: Separate baselines for each pressure condition (excluding outliers)
            P1_baseline_trials = intersect(P1_clean, isExperiment);
            P3_baseline_trials = intersect(P3_clean, isExperiment);
            P6_baseline_trials = intersect(P6_clean, isExperiment);
            
            if ~isempty(P1_baseline_trials)
                baseline_averages.P1 = mean(cat(2, power_baseline_mean{P1_baseline_trials}), 2);
                fprintf('P1 baseline from %d trials\n', length(P1_baseline_trials));
            else
                baseline_averages.P1 = mean(cat(2, power_baseline_mean{:}), 2);
                fprintf('P1 baseline: fallback to all trials\n');
            end
            
            if ~isempty(P3_baseline_trials)
                baseline_averages.P3 = mean(cat(2, power_baseline_mean{P3_baseline_trials}), 2);
                fprintf('P3 baseline from %d trials\n', length(P3_baseline_trials));
            else
                baseline_averages.P3 = mean(cat(2, power_baseline_mean{:}), 2);
                fprintf('P3 baseline: fallback to all trials\n');
            end
            
            if ~isempty(P6_baseline_trials)
                baseline_averages.P6 = mean(cat(2, power_baseline_mean{P6_baseline_trials}), 2);
                fprintf('P6 baseline from %d trials\n', length(P6_baseline_trials));
            else
                baseline_averages.P6 = mean(cat(2, power_baseline_mean{:}), 2);
                fprintf('P6 baseline: fallback to all trials\n');
            end
            
        otherwise
            error('Invalid baseline_mode. Use ''common'' or ''per_condition''');
    end
end

function [signal, times] = extract_trial_data(subject, trial_data, ic)
    % Extract signal and time data from trial
    if subject < 10
        sources = trial_data.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
        times = trial_data.EEG_stream.Preprocessed.Time_Domain.Times;
    else
        sources = trial_data.EEG_stream.Preprocessed.Sources;
        times = trial_data.EEG_stream.Preprocessed.Times;
    end
    
    % Concatenate all source segments for this component
    signal = cell2mat(cellfun(@(x) x(ic, :), sources, 'UniformOutput', false));
end

function [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed(...
    trial_max_power, isExperiment, P1, P3, P6, threshold_factor)
    
    % Only consider experimental trials for each pressure group
    P1_exp = intersect(P1, isExperiment);
    P3_exp = intersect(P3, isExperiment);
    P6_exp = intersect(P6, isExperiment);
    
    % Extract max powers for each experimental pressure group
    max_power_P1 = trial_max_power(P1_exp);
    max_power_P3 = trial_max_power(P3_exp);
    max_power_P6 = trial_max_power(P6_exp);
    
    % Find outliers within experimental trials only
    outliers_P1 = isoutlier(max_power_P1, "median", "ThresholdFactor", threshold_factor);
    outliers_P3 = isoutlier(max_power_P3, "median", "ThresholdFactor", threshold_factor);
    outliers_P6 = isoutlier(max_power_P6, "median", "ThresholdFactor", threshold_factor);
    
    % Remove outliers from experimental trials only
    P1_clean = P1_exp(~outliers_P1);
    P3_clean = P3_exp(~outliers_P3);
    P6_clean = P6_exp(~outliers_P6);
    
    fprintf('Outliers removed (experimental trials only): P1: %d/%d, P3: %d/%d, P6: %d/%d\n', ...
        sum(outliers_P1), length(P1_exp), sum(outliers_P3), length(P3_exp), sum(outliers_P6), length(P6_exp));
end

% function plot_grand_means_from_precomputed(trial_powers, freqs, isExperiment, ...
%     P1, P3, P6, power_baseline_mean, freqs_baseline, subject, ic)
% 
%     % Collect power data for each pressure group
%     powers_P1 = [];
%     powers_P3 = [];
%     powers_P6 = [];
% 
%     for trial = 1:length(trial_powers)
%         if ~ismember(trial, isExperiment) || isempty(trial_powers{trial})
%             continue;
%         end
% 
%         power = trial_powers{trial};
% 
%         if ismember(trial, P1)
%             powers_P1 = [powers_P1, power]; %#ok<AGROW>
%         elseif ismember(trial, P3)
%             powers_P3 = [powers_P3, power]; %#ok<AGROW>
%         elseif ismember(trial, P6)
%             powers_P6 = [powers_P6, power]; %#ok<AGROW>
%         end
%     end
% 
%     % Calculate grand means for each condition
%     grand_mean_P1 = mean(powers_P1, 2);
%     grand_mean_P3 = mean(powers_P3, 2);
%     grand_mean_P6 = mean(powers_P6, 2);
% 
%     % Calculate baseline means if available
%     baseline_P1 = [];
%     baseline_P3 = [];
%     baseline_P6 = [];
%     baseline_total = [];
% 
%     if ~isempty(power_baseline_mean)
%         % Get baseline powers for each condition
%         P1_baseline_trials = intersect(P1, isExperiment);
%         P3_baseline_trials = intersect(P3, isExperiment);
%         P6_baseline_trials = intersect(P6, isExperiment);
% 
%         if ~isempty(P1_baseline_trials)
%             baseline_P1 = mean(cat(2, power_baseline_mean{P1_baseline_trials}), 2);
%         end
%         if ~isempty(P3_baseline_trials)
%             baseline_P3 = mean(cat(2, power_baseline_mean{P3_baseline_trials}), 2);
%         end
%         if ~isempty(P6_baseline_trials)
%             baseline_P6 = mean(cat(2, power_baseline_mean{P6_baseline_trials}), 2);
%         end
% 
%         % Total baseline from all experimental trials
%         all_experimental_trials = [P1_baseline_trials, P3_baseline_trials, P6_baseline_trials];
%         if ~isempty(all_experimental_trials)
%             baseline_total = mean(cat(2, power_baseline_mean{all_experimental_trials}), 2);
%         else
%             % Fallback to all baseline trials
%             baseline_total = mean(cat(2, power_baseline_mean{:}), 2);
%         end
%     end
% 
%     % Define colors for each condition
%     color_P1 = [0.2, 0.4, 0.8];  % Blue
%     color_P3 = [0.8, 0.4, 0.2];  % Orange/Red
%     color_P6 = [0.2, 0.7, 0.3];  % Green
%     color_total = [0, 0, 0];     % Black
% 
%     % Create the plot
%     figure('Position', [100, 100, 800, 600]);
%     hold on;
% 
%     % Plot grand means (solid lines)
%     h1 = plot(freqs, 10*log10(grand_mean_P1), 'Color', color_P1, 'LineWidth', 2, ...
%         'DisplayName', 'P1 Grand Mean');
%     h2 = plot(freqs, 10*log10(grand_mean_P3), 'Color', color_P3, 'LineWidth', 2, ...
%         'DisplayName', 'P3 Grand Mean');
%     h3 = plot(freqs, 10*log10(grand_mean_P6), 'Color', color_P6, 'LineWidth', 2, ...
%         'DisplayName', 'P6 Grand Mean');
% 
%     % Plot baselines (dashed lines with matching colors) - using freqs_baseline for x-axis
%     if ~isempty(baseline_P1)
%         h4 = plot(freqs_baseline, 10*log10(baseline_P1), '--', 'Color', color_P1, 'LineWidth', 1, ...
%             'DisplayName', 'P1 Baseline');
%     end
%     if ~isempty(baseline_P3)
%         h5 = plot(freqs_baseline, 10*log10(baseline_P3), '--', 'Color', color_P3, 'LineWidth', 1, ...
%             'DisplayName', 'P3 Baseline');
%     end
%     if ~isempty(baseline_P6)
%         h6 = plot(freqs_baseline, 10*log10(baseline_P6), '--', 'Color', color_P6, 'LineWidth', 1, ...
%             'DisplayName', 'P6 Baseline');
%     end
% 
%     % Plot total baseline (black dashed) - using freqs_baseline for x-axis
%     if ~isempty(baseline_total)
%         h7 = plot(freqs_baseline, 10*log10(baseline_total), '--', 'Color', color_total, 'LineWidth', 1, ...
%             'DisplayName', 'Total Baseline');
%     end
% 
%     % Formatting
%     ylabel('Power (dB)', 'FontSize', 14);
%     xlabel('Frequency (Hz)', 'FontSize', 14);
%     title(sprintf('Grand Means and Baselines - Subject %d, IC %d', subject, ic), 'FontSize', 16, 'FontWeight', 'normal');
%     legend('Location', 'best', 'FontSize', 10);
%     grid on;
%     grid minor;
% 
%     % Set axis properties - use combined frequency range
%     set(gca, 'FontSize', 10);
%     if ~isempty(freqs_baseline)
%         xlim([min([freqs; freqs_baseline]), max([freqs; freqs_baseline])]);
%     else
%         xlim([min(freqs), max(freqs)]);
%     end
% 
%     % Add some styling
%     box on;
%     set(gca, 'LineWidth', 1);
% 
%     hold off;
% 
%     % % Print summary information
%     % fprintf('Plotted grand means and baselines:\n');
%     % fprintf('  P1: %d trials, baseline from %d trials\n', size(powers_P1,2), length(intersect(P1, isExperiment)));
%     % fprintf('  P3: %d trials, baseline from %d trials\n', size(powers_P3,2), length(intersect(P3, isExperiment)));
%     % fprintf('  P6: %d trials, baseline from %d trials\n', size(powers_P6,2), length(intersect(P6, isExperiment)));
% end

function [power_warped, median_flx, median_ext] = perform_time_warping(power, times, trial_info)
    % Extract event information
    ev = trial_info.Events.EEG_stream.Preprocessed;
    flxS = ev.flextoflex_start_indx;
    ext = ev.flexion_end_indx;
    flxE = ev.flextoflex_end_indx;
    
    events = [reshape(flxS, [], 1), ...
              reshape(ext(1:length(flxS)), [], 1), ...
              reshape(flxE, [], 1)];
    
    % Calculate median lengths
    flex_lens = events(:, 2) - events(:, 1) + 1;
    ext_lens = events(:, 3) - (events(:, 2) + 1) + 1;
    median_flx = round(median(flex_lens));
    median_ext = round(median(ext_lens));
    
    % Split power into epochs
    time_sizes = cellfun(@(x) size(x, 2), times);
    c = [0, cumsum(time_sizes)];
    n_epochs = length(times);
    
    power_epochs = cell(1, n_epochs);
    for i = 1:n_epochs
        power_epochs{i} = power(:, c(i)+1:c(i+1));
    end
    
    % Time warp each epoch
    power_warped = cell(size(power_epochs));
    for i = 1:n_epochs
        original_time = times{i};
        flx_len = flex_lens(i);
        
        % Flexion segment
        flx_seg = power_epochs{i}(:, 1:flx_len);
        flx_t = original_time(1:flx_len);
        flx_target_t = linspace(flx_t(1), flx_t(end), median_flx);
        flx_interp = interp1(flx_t', flx_seg', flx_target_t, 'linear', 'extrap')';
        
        % Extension segment
        ext_seg = power_epochs{i}(:, flx_len+1:end);
        ext_t = original_time(flx_len+1:end);
        ext_target_t = linspace(ext_t(1), ext_t(end), median_ext);
        ext_interp = interp1(ext_t', ext_seg', ext_target_t, 'linear', 'extrap')';
        
        % Combine segments
        power_warped{i} = [flx_interp, ext_interp];
    end
end

function power_warped = apply_baseline_correction(power_warped, baseline_avg, freqs_baseline)
    if ~isempty(baseline_avg)
        % Apply dB change baseline correction
        for i = 1:length(power_warped)
            freq_subset = 1:length(freqs_baseline);
            power_warped{i} = 10 * log10(power_warped{i}(freq_subset, :) ./ baseline_avg);
        end
    end
end

function data_clean = remove_empty_columns(data)
    non_empty_cols = ~all(cellfun(@isempty, data), 1);
    data_clean = data(:, non_empty_cols);
end