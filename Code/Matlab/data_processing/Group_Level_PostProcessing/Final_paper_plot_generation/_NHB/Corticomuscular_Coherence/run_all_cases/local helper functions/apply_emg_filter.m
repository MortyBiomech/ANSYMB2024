function EMG_filt = apply_emg_filter(EMG_raw, fs, fcond)
% APPLY_EMG_FILTER  Apply band-pass, high-pass, or no filter to raw EMG.
%   Input / output: [nChan x nSamples]
    switch fcond.type
        case 'bandpass'
            [b, a]   = butter(4, fcond.cutoff / (fs/2), 'bandpass');
            EMG_filt = filtfilt(b, a, double(EMG_raw'))';
        case 'highpass'
            [b, a]   = butter(4, fcond.cutoff / (fs/2), 'high');
            EMG_filt = filtfilt(b, a, double(EMG_raw'))';
        case 'none'
            EMG_filt = double(EMG_raw);
        otherwise
            error('apply_emg_filter: unknown filter type "%s".', fcond.type);
    end
end