function PSD_integral = ...
    PSD_integral_feature_generator_Flx_Ext_separate(condition_indices, data, IC, ...
    frequencies, per_trial_or_all_epochs, pressure_score)

    PSD_integral = struct();

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
        PSD_integral.(P) = [];

        for i = 1:length(condition_indices.(P))
        
            signal_Flx = squeeze(data{1, condition_indices.(P)(i)}.EEG_stream...
                .Preprocessed.Freq_Domain.Sources.Flexion(IC, :, :))';
            PSD_integral_signal_Flx = ...
                [ trapz(frequencies(freq_bands_indx(1, 1):freq_bands_indx(1, 2)), ...
                        signal_Flx(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(2, 1):freq_bands_indx(2, 2)), ...
                        signal_Flx(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(3, 1):freq_bands_indx(3, 2)), ...
                        signal_Flx(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(4, 1):freq_bands_indx(4, 2)), ...
                        signal_Flx(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(5, 1):freq_bands_indx(5, 2)), ...
                        signal_Flx(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2) ];

            


            signal_Ext = squeeze(data{1, condition_indices.(P)(i)}.EEG_stream...
                .Preprocessed.Freq_Domain.Sources.Extension(IC, :, :))';
            PSD_integral_signal_Ext = ...
                [ trapz(frequencies(freq_bands_indx(1, 1):freq_bands_indx(1, 2)), ...
                        signal_Ext(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(2, 1):freq_bands_indx(2, 2)), ...
                        signal_Ext(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(3, 1):freq_bands_indx(3, 2)), ...
                        signal_Ext(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(4, 1):freq_bands_indx(4, 2)), ...
                        signal_Ext(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2), ...
                  trapz(frequencies(freq_bands_indx(5, 1):freq_bands_indx(5, 2)), ...
                        signal_Ext(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2) ];

            
            

            mean_PSD_integral_signal_Flx = mean(PSD_integral_signal_Flx, 1);
            mean_PSD_integral_signal_Ext = mean(PSD_integral_signal_Ext, 1);
           

            if strcmp(per_trial_or_all_epochs, 'all_epochs')
                
                PSD_integral.(P) = cat(1, PSD_integral.(P), [PSD_integral_signal_Flx, PSD_integral_signal_Ext]);
        
            elseif strcmp(per_trial_or_all_epochs, 'per_trial')
                
                PSD_integral.(P) = cat(1, PSD_integral.(P), [mean_PSD_integral_signal_Flx, mean_PSD_integral_signal_Ext]);
        
            end
    
            
        end

    end

    PSD_integral = {PSD_integral};
    
end