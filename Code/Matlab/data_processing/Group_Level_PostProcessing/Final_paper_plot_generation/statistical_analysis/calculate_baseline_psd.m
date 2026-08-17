function [psd_baseline_mean, freqs_baseline] = calculate_baseline_psd(...
    subject, EEGLAB_events, EEGLAB_source, sampling_rate, freq_limits, ...
    window_length, overlap, nfft, baseline_window)
    
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
    
    % Calculate PSD for all baselines using pwelch
    psd_results = cell(1, n_events);
    for i = 1:n_events
        [pxx, freqs_baseline] = pwelch(baselines{i}, window_length, overlap, nfft, sampling_rate);
        % Filter frequencies to desired range
        freq_mask = freqs_baseline >= freq_limits(1) & freqs_baseline <= freq_limits(2);
        psd_results{i} = pxx(freq_mask);
        freqs_baseline = freqs_baseline(freq_mask);
    end
    
    
    psd_baseline_mean = psd_results;
end