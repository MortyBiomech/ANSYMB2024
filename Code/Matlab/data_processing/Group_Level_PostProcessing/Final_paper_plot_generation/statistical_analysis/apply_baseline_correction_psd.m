function psd_corrected = apply_baseline_correction_psd(psd, baseline_avg)
    % if ~isempty(baseline_avg)
    %     % Apply dB change baseline correction for PSD
    %     psd_corrected = 10 * log10(psd ./ baseline_avg);
    % else
    %     psd_corrected = psd;
    % end

    psd_corrected = psd;
   
end