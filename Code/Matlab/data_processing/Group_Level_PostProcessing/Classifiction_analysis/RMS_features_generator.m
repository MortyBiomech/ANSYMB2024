function [RMS_nonNorm, RMS_Norm] = ...
    RMS_features_generator(condition_indices, data, IC, frequencies)

    RMS_nonNorm = struct();
    RMS_Norm = struct();

    freq_bands(1, :) = [2, 4];   % Delta
    freq_bands(2, :) = [4, 8];   % Theta
    freq_bands(3, :) = [8, 13];  % Alpha
    freq_bands(4, :) = [13, 30]; % Beta
    freq_bands(5, :) = [30, 40]; % Gamma


    freq_bands_indx = zeros(5, 2);
    for i = 1:5 % 5 frequency bands
        freq_bands_indx(i, 1) = find(frequencies > freq_bands(i, 1), 1);
        freq_bands_indx(i, 2) = find(frequencies < freq_bands(i, 2), 1, 'last');
    end

    
    for Pn = {'P1', 'P3', 'P6'}
        
        P = Pn{1, 1};
        RMS_nonNorm.(P) = [];
        RMS_Norm.(P)    = [];
        for i = 1:length(condition_indices.(P))
        
            signal = squeeze(data{1, condition_indices.(P)(i)}.EEG_stream...
                .Preprocessed.Freq_Domain.Sources(IC, :, :))';
            
            RMS_nonNorm.(P) = cat(1, RMS_nonNorm.(P), ...
                [ rms(signal(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2), ...
                  rms(signal(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2), ...
                  rms(signal(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2), ...
                  rms(signal(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2), ...
                  rms(signal(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2) ]);
        
    
            denominator = rms(signal(:, freq_bands_indx(1):freq_bands_indx(end)), 2);
            RMS_Norm.(P) = cat(1, RMS_Norm.(P), ...
                [ rms(signal(:, freq_bands_indx(1, 1):freq_bands_indx(1, 2)), 2)./ denominator, ...
                  rms(signal(:, freq_bands_indx(2, 1):freq_bands_indx(2, 2)), 2)./ denominator, ...
                  rms(signal(:, freq_bands_indx(3, 1):freq_bands_indx(3, 2)), 2)./ denominator, ...
                  rms(signal(:, freq_bands_indx(4, 1):freq_bands_indx(4, 2)), 2)./ denominator, ...
                  rms(signal(:, freq_bands_indx(5, 1):freq_bands_indx(5, 2)), 2)./ denominator ]);
        
        end

    end

    RMS_nonNorm = {RMS_nonNorm};
    RMS_Norm = {RMS_Norm};
    
end