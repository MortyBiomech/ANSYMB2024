function [data1, frequencies] = calculate_PSD(data0, Trials_Info)

    Fs = 500; % Sampling frequency
    nfft = 2048; % number of FFT points

    base_struct = struct('EEG_stream', struct('Preprocessed', ...
        struct('Freq_Domain', struct('Sources', []))));
    data1 = repmat({base_struct}, 1, length(data0));

    signal = data0{1, 1}.EEG_stream.Preprocessed.Sources;
    window = floor(size(signal{1, 1}, 2)/2);
    noverlap = floor(0.9*window);
    [~, freqs] = ...
        pwelch(signal{1,1}(1,:)', window, noverlap, nfft, Fs);
    frequencies = freqs';

    for i = 1:length(data0)
        if strcmp(Trials_Info{1, i}.General.Description, 'Experiment')
            signal = data0{1, i}.EEG_stream.Preprocessed.Sources;
            if ~isempty(signal)
                window = floor(cellfun(@(x) size(x, 2), signal)/2); % length of each segment
                noverlap = floor(0.9*window); % number of samples to overlap between segments
            
                psdCell = cellfun(@(x, win, ol) ...
                    pwelch(x', win, ol, nfft, Fs)', ...
                    signal, num2cell(window), num2cell(noverlap), ...
                    'UniformOutput', false);
                data1{1, i}.EEG_stream.Preprocessed.Freq_Domain.Sources = ...
                    cat(3, psdCell{:});
            end
        end
    end

end