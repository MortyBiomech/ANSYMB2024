function signal = extract_trial_data_psd(subject, trial_data, ic)
    % Extract signal data from trial for PSD analysis
    if subject < 10
        sources = trial_data.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
    else
        sources = trial_data.EEG_stream.Preprocessed.Sources;
    end
    
    % Concatenate all source segments for this component
    signal = cell2mat(cellfun(@(x) x(ic, :), sources, 'UniformOutput', false));
end