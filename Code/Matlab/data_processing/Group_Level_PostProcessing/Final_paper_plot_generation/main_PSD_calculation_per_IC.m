function [PSD_struct, freqs, freqs_baseline] = main_PSD_calculation_per_IC(subject, ic, trial_info, subject_data, EEGLAB_source, EEGLAB_events)
    % Power spectral density calculation per independent component using pwelch
    % Key difference: PSD computed using Welch's method instead of CWT
    
    %% Configuration parameters
    SAMPLING_RATE = 500;
    FREQ_LIMITS = [1 70];
    BASELINE_WINDOW = 1000; % samples before event
    OUTLIER_THRESHOLD = 3;
    
    % Pwelch parameters
    WINDOW_LENGTH = SAMPLING_RATE * 2; % 3 seconds window
    OVERLAP = WINDOW_LENGTH * 0.9; % 90% overlap
    NFFT = 2^nextpow2(WINDOW_LENGTH); % FFT length
    
    % Baseline correction options:
    % 'common' - One baseline average from all experimental trials (excluding outliers)
    % 'per_condition' - Separate baseline averages for P1, P3, P6 (excluding outliers)
    BASELINE_MODE = 'common';  % Switch between 'common' and 'per_condition'
    
    %% Baseline calculations
    [psd_baseline_mean, freqs_baseline] = calculate_baseline_psd(...
        subject, EEGLAB_events, EEGLAB_source, SAMPLING_RATE, FREQ_LIMITS, ...
        WINDOW_LENGTH, OVERLAP, NFFT, BASELINE_WINDOW);
    
    %% Get experimental trials and pressure groups
    [P1, P3, P6, isExperiment] = get_pressure_groups(subject, trial_info);
    
    %% Single-pass processing: PSD + outlier detection + grand means + final processing
    [PSD_struct, freqs] = process_all_trials_psd_single_pass(subject, ic, trial_info, ...
        subject_data, isExperiment, P1, P3, P6, psd_baseline_mean, freqs_baseline, ...
        SAMPLING_RATE, FREQ_LIMITS, WINDOW_LENGTH, OVERLAP, NFFT, OUTLIER_THRESHOLD, BASELINE_MODE);
    
end

%% Helper Functions
% function [psd_baseline_mean, freqs_baseline] = calculate_baseline_psd(...
%     subject, EEGLAB_events, EEGLAB_source, sampling_rate, freq_limits, ...
%     window_length, overlap, nfft, baseline_window)
% 
%     if subject >= 10
%         event_type = 'SM_Start_Move';
%     else
%         event_type = 'SB_Start_Beep';
%     end
% 
%     % Find event indices
%     event_ids = strcmp({EEGLAB_events.type}, event_type);
%     latencies = [EEGLAB_events(event_ids).latency];
% 
%     % Extract baseline segments
%     n_events = length(latencies);
%     baselines = cell(1, n_events);
% 
%     for i = 1:n_events
%         start_idx = latencies(i) - baseline_window;
%         end_idx = latencies(i) - 1;
%         baselines{i} = EEGLAB_source(start_idx:end_idx);
%     end
% 
%     % Calculate PSD for all baselines using pwelch
%     psd_results = cell(1, n_events);
%     for i = 1:n_events
%         [pxx, freqs_baseline] = pwelch(baselines{i}, window_length, overlap, nfft, sampling_rate);
%         % Filter frequencies to desired range
%         freq_mask = freqs_baseline >= freq_limits(1) & freqs_baseline <= freq_limits(2);
%         psd_results{i} = pxx(freq_mask);
%         freqs_baseline = freqs_baseline(freq_mask);
%     end
% 
% 
%     psd_baseline_mean = psd_results;
% end

% function [P1, P3, P6, isExperiment] = get_pressure_groups(subject, trial_info)
%     % Get experimental trials
%     if subject < 10
%         isExperiment = 1:length(trial_info);
%     else
%         isExperiment = find(cellfun(@(x) strcmp(x.General.Description, 'Experiment'), trial_info));
%     end
% 
%     % Extract pressure levels more efficiently
%     pressures = cellfun(@(x) x.General.Pressure, trial_info);
% 
%     P1 = intersect(find(pressures == 1), isExperiment);
%     P3 = intersect(find(pressures == 3), isExperiment);
%     P6 = intersect(find(pressures == 6), isExperiment);
% end

% function [PSD_struct, freqs] = process_all_trials_psd_single_pass(subject, ic, trial_info, ...
%     subject_data, isExperiment, P1, P3, P6, psd_baseline_mean, freqs_baseline, ...
%     sampling_rate, freq_limits, window_length, overlap, nfft, outlier_threshold, baseline_mode)
% 
%     % Pre-allocate storage for all trials
%     n_trials = length(subject_data);
%     trial_psds = cell(n_trials, 1);           % Raw PSD data
%     trial_max_psd = zeros(n_trials, 1);       % For outlier detection
%     freqs = [];
% 
%     %% SINGLE PASS: Compute PSD once per trial and extract all needed info
%     fprintf('Processing %d out of %d trials with single PSD calculation per trial...\n', ...
%         sum(ismember(1:n_trials, isExperiment)), n_trials);
% 
%     for trial = 1:n_trials
%         if ~ismember(trial, isExperiment)
%             continue;
%         end
% 
%         % Extract signal
%         signal = extract_trial_data_psd(subject, subject_data{trial}, ic);
% 
%         % *** SINGLE PSD COMPUTATION PER TRIAL USING PWELCH ***
%         [pxx, freqs] = pwelch(signal, window_length, overlap, nfft, sampling_rate);
% 
%         % Filter frequencies to desired range
%         freq_mask = freqs >= freq_limits(1) & freqs <= freq_limits(2);
%         psd = pxx(freq_mask);
%         freqs = freqs(freq_mask);
% 
%         % Store data for later processing
%         trial_psds{trial} = psd;
%         trial_max_psd(trial) = max(psd);
%     end
% 
%     %% Outlier detection using pre-computed max PSDs
%     [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed_psd(...
%         trial_max_psd, isExperiment, P1, P3, P6, outlier_threshold);
% 
%     %% Final processing: Baseline correction
%     baseline_averages = calculate_baseline_averages_psd(psd_baseline_mean, ...
%         isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode);
% 
%     % Pre-allocate final storage
%     PSD_P1 = cell(1, length(P1_clean));
%     PSD_P3 = cell(1, length(P3_clean));
%     PSD_P6 = cell(1, length(P6_clean));
% 
%     P1_idx = 1; P3_idx = 1; P6_idx = 1;
% 
%     fprintf('Applying baseline correction (mode: %s)...\n', baseline_mode);
% 
%     for trial = 13:n_trials
%         if ~ismember(trial, isExperiment)
%             continue;
%         end
% 
%         % Skip if this trial was marked as outlier
%         if ~(ismember(trial, P1_clean) || ismember(trial, P3_clean) || ismember(trial, P6_clean))
%             continue;
%         end
% 
%         % Use pre-computed PSD data
%         psd = trial_psds{trial};
% 
%         % Baseline correction
%         if ismember(trial, P1_clean)
%             psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P1);
%         elseif ismember(trial, P3_clean)
%             psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P3);
%         elseif ismember(trial, P6_clean)
%             psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P6);
%         end
% 
%         % Store in appropriate pressure group
%         if ismember(trial, P1_clean)
%             PSD_P1{P1_idx} = psd_corrected;
%             P1_idx = P1_idx + 1;
%         elseif ismember(trial, P3_clean)
%             PSD_P3{P3_idx} = psd_corrected;
%             P3_idx = P3_idx + 1;
%         elseif ismember(trial, P6_clean)
%             PSD_P6{P6_idx} = psd_corrected;
%             P6_idx = P6_idx + 1;
%         end
%     end
% 
%     % Remove empty columns
%     PSD_P1 = remove_empty_cells(PSD_P1);
%     PSD_P3 = remove_empty_cells(PSD_P3);
%     PSD_P6 = remove_empty_cells(PSD_P6);
% 
%     % Create output structure
%     PSD_struct = struct('P1', {PSD_P1}, 'P3', {PSD_P3}, 'P6', {PSD_P6});
% 
%     fprintf('Processing complete. P1: %d trials, P3: %d trials, P6: %d trials\n', ...
%         length(PSD_P1), length(PSD_P3), length(PSD_P6));
% end

% function baseline_averages = calculate_baseline_averages_psd(psd_baseline_mean, ...
%     isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode)
%     % Calculate baseline averages based on the specified mode
% 
%     if isempty(psd_baseline_mean)
%         baseline_averages = struct('P1', [], 'P3', [], 'P6', []);
%         return;
%     end
% 
%     switch baseline_mode
%         case 'common'
%             % Option A: One common baseline from all experimental trials (excluding outliers)
%             all_clean_trials = [P1_clean, P3_clean, P6_clean];
%             baseline_trials = intersect(all_clean_trials, isExperiment);
% 
%             if ~isempty(baseline_trials)
%                 common_baseline = mean(cat(2, psd_baseline_mean{baseline_trials}), 2);
%                 baseline_averages.P1 = common_baseline;
%                 baseline_averages.P3 = common_baseline;
%                 baseline_averages.P6 = common_baseline;
%                 fprintf('Using common baseline from %d trials\n', length(baseline_trials));
%             else
%                 % Fallback to all baselines if no clean trials found
%                 common_baseline = mean(cat(2, psd_baseline_mean{:}), 2);
%                 baseline_averages.P1 = common_baseline;
%                 baseline_averages.P3 = common_baseline;
%                 baseline_averages.P6 = common_baseline;
%                 fprintf('Using common baseline from all %d baseline trials (fallback)\n', length(psd_baseline_mean));
%             end
% 
%         case 'per_condition'
%             % Option B: Separate baselines for each pressure condition (excluding outliers)
%             P1_baseline_trials = intersect(P1_clean, isExperiment);
%             P3_baseline_trials = intersect(P3_clean, isExperiment);
%             P6_baseline_trials = intersect(P6_clean, isExperiment);
% 
%             if ~isempty(P1_baseline_trials)
%                 baseline_averages.P1 = mean(cat(2, psd_baseline_mean{P1_baseline_trials}), 2);
%                 fprintf('P1 baseline from %d trials\n', length(P1_baseline_trials));
%             else
%                 baseline_averages.P1 = mean(cat(2, psd_baseline_mean{:}), 2);
%                 fprintf('P1 baseline: fallback to all trials\n');
%             end
% 
%             if ~isempty(P3_baseline_trials)
%                 baseline_averages.P3 = mean(cat(2, psd_baseline_mean{P3_baseline_trials}), 2);
%                 fprintf('P3 baseline from %d trials\n', length(P3_baseline_trials));
%             else
%                 baseline_averages.P3 = mean(cat(2, psd_baseline_mean{:}), 2);
%                 fprintf('P3 baseline: fallback to all trials\n');
%             end
% 
%             if ~isempty(P6_baseline_trials)
%                 baseline_averages.P6 = mean(cat(2, psd_baseline_mean{P6_baseline_trials}), 2);
%                 fprintf('P6 baseline from %d trials\n', length(P6_baseline_trials));
%             else
%                 baseline_averages.P6 = mean(cat(2, psd_baseline_mean{:}), 2);
%                 fprintf('P6 baseline: fallback to all trials\n');
%             end
% 
%         otherwise
%             error('Invalid baseline_mode. Use ''common'' or ''per_condition''');
%     end
% end

% function signal = extract_trial_data_psd(subject, trial_data, ic)
%     % Extract signal data from trial for PSD analysis
%     if subject < 10
%         sources = trial_data.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
%     else
%         sources = trial_data.EEG_stream.Preprocessed.Sources;
%     end
% 
%     % Concatenate all source segments for this component
%     signal = cell2mat(cellfun(@(x) x(ic, :), sources, 'UniformOutput', false));
% end

% function [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed_psd(...
%     trial_max_psd, isExperiment, P1, P3, P6, threshold_factor)
% 
%     % Only consider experimental trials for each pressure group
%     P1_exp = intersect(P1, isExperiment);
%     P3_exp = intersect(P3, isExperiment);
%     P6_exp = intersect(P6, isExperiment);
% 
%     % Extract max PSDs for each experimental pressure group
%     max_psd_P1 = trial_max_psd(P1_exp);
%     max_psd_P3 = trial_max_psd(P3_exp);
%     max_psd_P6 = trial_max_psd(P6_exp);
% 
%     % Find outliers within experimental trials only
%     outliers_P1 = isoutlier(max_psd_P1, "median", "ThresholdFactor", threshold_factor);
%     outliers_P3 = isoutlier(max_psd_P3, "median", "ThresholdFactor", threshold_factor);
%     outliers_P6 = isoutlier(max_psd_P6, "median", "ThresholdFactor", threshold_factor);
% 
%     % Remove outliers from experimental trials only
%     P1_clean = P1_exp(~outliers_P1);
%     P3_clean = P3_exp(~outliers_P3);
%     P6_clean = P6_exp(~outliers_P6);
% 
%     fprintf('Outliers removed (experimental trials only): P1: %d/%d, P3: %d/%d, P6: %d/%d\n', ...
%         sum(outliers_P1), length(P1_exp), sum(outliers_P3), length(P3_exp), sum(outliers_P6), length(P6_exp));
% end

% function psd_corrected = apply_baseline_correction_psd(psd, baseline_avg)
%     if ~isempty(baseline_avg)
%         % Apply dB change baseline correction for PSD
%         psd_corrected = 10 * log10(psd ./ baseline_avg);
%     else
%         psd_corrected = psd;
%     end
% end

% function data_clean = remove_empty_cells(data)
%     non_empty_idx = ~cellfun(@isempty, data);
%     data_clean = data(non_empty_idx);
% end
