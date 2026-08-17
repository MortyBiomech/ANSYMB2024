function [PSD_struct, freqs] = process_all_trials_psd_single_pass(subject, ic, trial_info, ...
    subject_data, isExperiment, P1, P3, P6, psd_baseline_mean, freqs_baseline, ...
    sampling_rate, freq_limits, window_length, overlap, nfft, outlier_threshold, baseline_mode)
    
    % Pre-allocate storage for all trials
    n_trials = length(subject_data);
    trial_psds = cell(n_trials, 1);           % Raw PSD data
    trial_max_psd = zeros(n_trials, 1);       % For outlier detection
    freqs = [];
    
    %% SINGLE PASS: Compute PSD once per trial and extract all needed info
    fprintf('Processing %d out of %d trials with single PSD calculation per trial...\n', ...
        sum(ismember(1:n_trials, isExperiment)), n_trials);
    
    for trial = 1:n_trials
        if ~ismember(trial, isExperiment)
            continue;
        end
        
        % Extract signal
        signal = extract_trial_data_psd(subject, subject_data{trial}, ic);
        
        % *** SINGLE PSD COMPUTATION PER TRIAL USING PWELCH ***
        [pxx, freqs] = pwelch(signal, window_length, overlap, nfft, sampling_rate);
        pxx = 10*log10(pxx);
        
        % Filter frequencies to desired range
        freq_mask = freqs >= freq_limits(1) & freqs <= freq_limits(2);
        psd = pxx(freq_mask);
        freqs = freqs(freq_mask);
        
        % Store data for later processing
        trial_psds{trial} = psd;
        trial_max_psd(trial) = max(psd);
    end
    
    %% Outlier detection using pre-computed max PSDs
    [P1_clean, P3_clean, P6_clean] = remove_outliers_from_precomputed_psd(...
        trial_max_psd, isExperiment, P1, P3, P6, outlier_threshold);
    
    %% Final processing: Baseline correction
    baseline_averages = calculate_baseline_averages_psd(psd_baseline_mean, ...
        isExperiment, P1_clean, P3_clean, P6_clean, baseline_mode);
    
    % Pre-allocate final storage
    PSD_P1 = cell(1, length(P1_clean));
    PSD_P3 = cell(1, length(P3_clean));
    PSD_P6 = cell(1, length(P6_clean));
    
    P1_idx = 1; P3_idx = 1; P6_idx = 1;
    
    fprintf('Applying baseline correction (mode: %s)...\n', baseline_mode);
    
    for trial = 1:n_trials
        if ~ismember(trial, isExperiment)
            continue;
        end
        
        % Skip if this trial was marked as outlier
        if ~(ismember(trial, P1_clean) || ismember(trial, P3_clean) || ismember(trial, P6_clean))
            continue;
        end
        
        % Use pre-computed PSD data
        psd = trial_psds{trial};
        
        % Baseline correction
        if ismember(trial, P1_clean)
            psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P1);
        elseif ismember(trial, P3_clean)
            psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P3);
        elseif ismember(trial, P6_clean)
            psd_corrected = apply_baseline_correction_psd(psd, baseline_averages.P6);
        end
        
        % Store in appropriate pressure group
        if ismember(trial, P1_clean)
            PSD_P1{P1_idx} = psd_corrected;
            P1_idx = P1_idx + 1;
        elseif ismember(trial, P3_clean)
            PSD_P3{P3_idx} = psd_corrected;
            P3_idx = P3_idx + 1;
        elseif ismember(trial, P6_clean)
            PSD_P6{P6_idx} = psd_corrected;
            P6_idx = P6_idx + 1;
        end
    end
    
    % Remove empty columns
    PSD_P1 = remove_empty_cells(PSD_P1);
    PSD_P3 = remove_empty_cells(PSD_P3);
    PSD_P6 = remove_empty_cells(PSD_P6);
    
    % Create output structure
    PSD_struct = struct('P1', {PSD_P1}, 'P3', {PSD_P3}, 'P6', {PSD_P6}, 'baseline', {baseline_averages});
    
    fprintf('Processing complete. P1: %d trials, P3: %d trials, P6: %d trials\n', ...
        length(PSD_P1), length(PSD_P3), length(PSD_P6));
end