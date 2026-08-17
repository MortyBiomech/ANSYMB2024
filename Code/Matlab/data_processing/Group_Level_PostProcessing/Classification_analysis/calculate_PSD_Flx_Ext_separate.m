function [data1, frequencies] = calculate_PSD_Flx_Ext_separate(data0, Trials_Info, subject)


    Fs = 500; % Sampling frequency
    nfft = 2048; % number of FFT points

    base_struct = struct('EEG_stream', struct('Preprocessed', ...
        struct('Freq_Domain', struct('Sources', ...
        struct('Flexion', [], 'Extension', [], 'FlextoFlex', [])))));
    data1 = repmat({base_struct}, 1, length(data0));

    if subject < 10
        signal = data0{1, 1}.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized;
        fieldPath  = {'EEG_stream', 'Preprocessed', 'Time_Domain', 'Sources', 'Not_Length_Normalized'};
    else
        signal = data0{1, 1}.EEG_stream.Preprocessed.Sources;
        fieldPath  = {'EEG_stream', 'Preprocessed', 'Sources'};
    end
    window = floor(size(signal{1, 1}, 2)/2);
    noverlap = floor(0.7*window);
    [~, freqs] = ...
        pwelch(signal{1,1}(1,:)', window, noverlap, nfft, Fs);
    frequencies = freqs';

    for i = 1:length(data0)

        if (subject>=10) && ~strcmp(Trials_Info{1, i}.General.Description, 'Experiment')
            continue
        end
            
        signal = getfield(data0{1, i}, fieldPath{:});
        flexion_start_events = reshape(num2cell(Trials_Info{1, i}.Events.EEG_stream.Preprocessed.flextoflex_start_indx), size(signal));
        extension_start_events = reshape(num2cell(Trials_Info{1, i}.Events.EEG_stream.Preprocessed.extension_start_indx(1:length(flexion_start_events))), size(signal));
        
        c = cell2mat(extension_start_events) - cell2mat(flexion_start_events);
        if max(c) < 100 
            continue
        end
        
        
        signal_Flx = cellfun(@(x, a, b) x(:, 1:b-a+1), signal, flexion_start_events, extension_start_events, 'UniformOutput', false);
        signal_Ext = cellfun(@(x, a, b) x(:, b-a+2:end), signal, flexion_start_events, extension_start_events, 'UniformOutput', false);

        if ~isempty(signal)

            %% Flexion Part
            % window = floor(cellfun(@(x) size(x, 2), signal_Flx)/4); 
            % noverlap = floor(0.9*window); 
            window = floor(cellfun(@(x) size(x, 2), signal_Flx)/2); 
            noverlap = floor(0.7*window);
        
            psdCell = cellfun(@(x, win, ol) ...
                pwelch(x', win, ol, nfft, Fs)', ...
                signal_Flx, num2cell(window), num2cell(noverlap), ...
                'UniformOutput', false);
            data1{1, i}.EEG_stream.Preprocessed.Freq_Domain.Sources.Flexion = ...
                cat(3, psdCell{:});


            %% Extension Part
            % window = floor(cellfun(@(x) size(x, 2), signal_Ext)/4); 
            % noverlap = floor(0.9*window);
            window = floor(cellfun(@(x) size(x, 2), signal_Ext)/2); 
            noverlap = floor(0.7*window); 
        
            psdCell = cellfun(@(x, win, ol) ...
                pwelch(x', win, ol, nfft, Fs)', ...
                signal_Ext, num2cell(window), num2cell(noverlap), ...
                'UniformOutput', false);
            data1{1, i}.EEG_stream.Preprocessed.Freq_Domain.Sources.Extension = ...
                cat(3, psdCell{:});


            %% Flex_to_Flex Part
            window = floor(cellfun(@(x) size(x, 2), signal)/4); 
            noverlap = floor(0.7*window); 
        
            psdCell = cellfun(@(x, win, ol) ...
                pwelch(x', win, ol, nfft, Fs)', ...
                signal, num2cell(window), num2cell(noverlap), ...
                'UniformOutput', false);
            data1{1, i}.EEG_stream.Preprocessed.Freq_Domain.Sources.FlextoFlex = ...
                cat(3, psdCell{:});


        end

    end


end