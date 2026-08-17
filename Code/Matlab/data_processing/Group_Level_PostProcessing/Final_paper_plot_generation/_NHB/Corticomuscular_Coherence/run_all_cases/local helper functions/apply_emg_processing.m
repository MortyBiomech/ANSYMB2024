function EMG_proc = apply_emg_processing(EMG_filt, proc_label)
% APPLY_EMG_PROCESSING  Rectify or demodulate filtered EMG.
%   The Hilbert transform is applied along the time dimension.
%   Input / output: [nChan x nSamples]
    switch proc_label
        case 'abs_rect'
            % Full-wave rectification
            EMG_proc = abs(EMG_filt);

        case 'hilbert_amp'
            % Instantaneous amplitude: A(t) = |z(t)|  where z = analytic signal
            z        = hilbert(double(EMG_filt'))';   % hilbert expects columns
            EMG_proc = abs(z);

        case 'hilbert_demod'
            % Phase-only (demodulated): cos(phi(t)) = Re{ z(t) / |z(t)| }
            % Removes amplitude; unit amplitude, carries phase information only.
            z        = hilbert(double(EMG_filt'))';
            EMG_proc = real(z ./ abs(z));

        case 'none'
            EMG_proc = double(EMG_filt);

        otherwise
            error('apply_emg_processing: unknown proc label "%s".', proc_label);
    end
end