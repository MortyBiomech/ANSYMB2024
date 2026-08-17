function [RMS_nonNorm, RMS_Norm] = ...
    RMS_features_generator_Flx_Ext_separate(condition_indices, data, IC, ...
    frequencies, per_trial_or_all_epochs, pressure_score)

    RMS_nonNorm = struct();
    RMS_Norm = struct();

    freq_bands(1, :) = [2, 4];   % Delta
    freq_bands(2, :) = [4, 8];   % Theta
    freq_bands(3, :) = [8, 13];  % Alpha
    freq_bands(4, :) = [13, 30]; % Beta
    freq_bands(5, :) = [30, 50]; % Gamma


    freq_bands_indx = zeros(5, 2);
    for i = 1:5 % 5 frequency bands
        freq_bands_indx(i, 1) = find(frequencies > freq_bands(i, 1), 1);
        freq_bands_indx(i, 2) = find(frequencies < freq_bands(i, 2), 1, 'last');
    end

    if strcmp(pressure_score, 'pressure')
        X = {'P1', 'P3', 'P6'}; % three pressure conditions
    elseif strcmp(pressure_score, 'score')
        X = {'S1', 'S2', 'S3'}; % three score clusters
    end
    
    for Pn = X
        
        P = Pn{1, 1};
        RMS_nonNorm.(P) = [];
        RMS_Norm.(P)    = [];
        for i = 1:length(condition_indices.(P))
        
            signal_Flx = squeeze(data{1, condition_indices.(P)(i)}.EEG_stream...
                .Preprocessed.Freq_Domain.Sources.Flexion(IC, :, :))';
            rms_signal_Flx = ...
                [ rms(signal_Flx(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2), ...
                      rms(signal_Flx(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2), ...
                      rms(signal_Flx(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2), ...
                      rms(signal_Flx(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2), ...
                      rms(signal_Flx(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2) ];

            denominator_Flx = rms(signal_Flx(:, freq_bands_indx(1):freq_bands_indx(end)), 2);
            rms_signal_Flx_norm = ...
            [ rms(signal_Flx(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2)./ denominator_Flx, ...
                      rms(signal_Flx(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2)./ denominator_Flx, ...
                      rms(signal_Flx(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2)./ denominator_Flx, ...
                      rms(signal_Flx(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2)./ denominator_Flx, ...
                      rms(signal_Flx(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2)./ denominator_Flx ];



            signal_Ext = squeeze(data{1, condition_indices.(P)(i)}.EEG_stream...
                .Preprocessed.Freq_Domain.Sources.Extension(IC, :, :))';
            rms_signal_Ext = ...
                [ rms(signal_Ext(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2), ...
                      rms(signal_Ext(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2), ...
                      rms(signal_Ext(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2), ...
                      rms(signal_Ext(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2), ...
                      rms(signal_Ext(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2) ];

            denominator_Ext = rms(signal_Ext(:, freq_bands_indx(1):freq_bands_indx(end)), 2);
            rms_signal_Ext_norm = ...
            [ rms(signal_Ext(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2)./ denominator_Ext, ...
                      rms(signal_Ext(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2)./ denominator_Ext, ...
                      rms(signal_Ext(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2)./ denominator_Ext, ...
                      rms(signal_Ext(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2)./ denominator_Ext, ...
                      rms(signal_Ext(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2)./ denominator_Ext ];

            

            mean_rms_signal_Flx = mean(rms_signal_Flx, 1);
            mean_rms_signal_Ext = mean(rms_signal_Ext, 1);
            mean_rms_signal_Flx_norm = mean(rms_signal_Flx_norm, 1);
            mean_rms_signal_Ext_norm = mean(rms_signal_Ext_norm, 1);
            
            if strcmp(per_trial_or_all_epochs, 'all_epochs')
                
                RMS_nonNorm.(P) = cat(1, RMS_nonNorm.(P), [rms_signal_Flx, rms_signal_Ext]);
                RMS_Norm.(P) = cat(1, RMS_Norm.(P), [rms_signal_Flx_norm, rms_signal_Ext_norm]);
        
            elseif strcmp(per_trial_or_all_epochs, 'per_trial')
                
                RMS_nonNorm.(P) = cat(1, RMS_nonNorm.(P), [mean_rms_signal_Flx, mean_rms_signal_Ext]);
                RMS_Norm.(P) = cat(1, RMS_Norm.(P), [mean_rms_signal_Flx_norm, mean_rms_signal_Ext_norm]);
        
            end
    
            
        end

    end

    RMS_nonNorm = {RMS_nonNorm};
    RMS_Norm = {RMS_Norm};
    
end