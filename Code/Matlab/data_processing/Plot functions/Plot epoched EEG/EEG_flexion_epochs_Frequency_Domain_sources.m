function EEG_flexion_epochs_Frequency_Domain_sources(data, std_sem)

    All_colours = struct('dark_blue', [0, 0.4470, 0.7410], 'light_blue', [0.3010, 0.7450, 0.9330], ...
        'dark_orange', [0.8500, 0.3250, 0.0980], 'light_orange', [0.9290, 0.6940, 0.1250], ...
        'dark_green', [0.4660, 0.6740, 0.1880], 'light_green', [0.5960, 0.8740, 0.5410]);


    %% EEG_flexion - Frequency-Domain: sources
    %%% making 3d matrix 
    numsources = size(data.EEG_P1_flexion.Frequency_Domain.sources{1,1,1} , 1);
    
    numPointsP1   = size(data.EEG_P1_flexion.Frequency_Domain.sources{1,1,1} , 2);
    numEpochsP1   = size(data.EEG_P1_flexion.Frequency_Domain.sources , 3);
    signal_P1_3d = zeros(numsources, numPointsP1, numEpochsP1);
    countP1 = 0;
    for i = 1:numEpochsP1
        tempdata = data.EEG_P1_flexion.Frequency_Domain.sources{1,1, i};
        if ~isempty(tempdata)
            signal_P1_3d(:,:, i) = 10*log10(tempdata);
        else
            countP1 = countP1 + 1;
        end 
    end
    
    numPointsP3   = size(data.EEG_P3_flexion.Frequency_Domain.sources{1,1,1} , 2);
    numEpochsP3   = size(data.EEG_P3_flexion.Frequency_Domain.sources , 3);
    signal_P3_3d = zeros(numsources, numPointsP3, numEpochsP3);
    countP3 = 0;
    for i = 1:numEpochsP3
        tempdata = data.EEG_P3_flexion.Frequency_Domain.sources{1,1, i};
        if ~isempty(tempdata)
            signal_P3_3d(:,:, i) = 10*log10(tempdata);
        else
            countP3 = countP3 + 1;
        end
    end
    
    numPointsP6   = size(data.EEG_P6_flexion.Frequency_Domain.sources{1,1,1} , 2);
    numEpochsP6   = size(data.EEG_P6_flexion.Frequency_Domain.sources , 3);
    signal_P6_3d = zeros(numsources, numPointsP6, numEpochsP6);
    countP6 = 0;
    for i = 1:numEpochsP6
        tempdata = data.EEG_P6_flexion.Frequency_Domain.sources{1,1, i};
        if ~isempty(tempdata)
            signal_P6_3d(:,:, i) = 10*log10(tempdata);
        else
            countP6 = countP6 + 1;
        end
    end
    
    
    %
    XP1 = linspace(0, 100, numPointsP1);
    XP3 = linspace(0, 100, numPointsP3);
    XP6 = linspace(0, 100, numPointsP6);
    
    mean_signalP1 = mean(signal_P1_3d, 3, 'omitnan');
    mean_signalP3 = mean(signal_P3_3d, 3, 'omitnan');
    mean_signalP6 = mean(signal_P6_3d, 3, 'omitnan');
    
    std_signalP1 = std(signal_P1_3d, 0, 3, 'omitnan');
    std_signalP3 = std(signal_P3_3d, 0, 3, 'omitnan');
    std_signalP6 = std(signal_P6_3d, 0, 3, 'omitnan');

    sem_signalP1 = std_signalP1 ./ sqrt( sum(~isnan(signal_P1_3d), 3) - countP1);
    sem_signalP3 = std_signalP3 ./ sqrt( sum(~isnan(signal_P3_3d), 3) - countP3);
    sem_signalP6 = std_signalP6 ./ sqrt( sum(~isnan(signal_P6_3d), 3) - countP6);

    if strcmp(std_sem, 'STD')
        std_sem_signalP1 = std_signalP1;
        std_sem_signalP3 = std_signalP3;
        std_sem_signalP6 = std_signalP6;
    elseif strcmp(std_sem, 'SEM')
        std_sem_signalP1 = sem_signalP1;
        std_sem_signalP3 = sem_signalP3;
        std_sem_signalP6 = sem_signalP6;
    end
    
    % plot the results
    h1 = []; h2 = []; h3 = [];
    XLabel_n = 'Frequency [Hz]';
    YLabel_n = 'Power 10*log_{10}(\muV^2/Hz)';
    XLim_n = [0.5 50];
    figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
    tiledlayout(8,8)
    sgtitle(sprintf('Frequency-Domain (Flexion, Mean $\\pm$ %s): %d epochs P1, %d epochs P3, %d epochs P6', ...
        std_sem, numEpochsP1, numEpochsP3, numEpochsP6), 'Interpreter', 'latex');
    for i = 1:61
        nexttile; hold on
    
        % fill([XP1 fliplr(XP1)], [mean_signalP1(i, :) + std_signalP1(i, :), ...
        %     fliplr(mean_signalP1(i, :) - std_signalP1(i, :))], ...
        %     All_colours.light_blue, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h1 = plot(XP1, mean_signalP1(i, :), 'Color', All_colours.dark_blue);
    
        % fill([XP3 fliplr(XP3)], [mean_signalP3(i, :) + std_signalP3(i, :), ...
        %     fliplr(mean_signalP3(i, :) - std_signalP3(i, :))], ...
        %     All_colours.light_orange, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h2 = plot(XP3, mean_signalP3(i, :), 'Color', All_colours.dark_orange);
    
        % fill([XP6 fliplr(XP6)], [mean_signalP6(i, :) + std_signalP6(i, :), ...
        %     fliplr(mean_signalP6(i, :) - std_signalP6(i, :))], ...
        %     All_colours.light_green, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h3 = plot(XP6, mean_signalP6(i, :), 'Color', All_colours.dark_green);
        
        set(gca, 'XLim', [0.5 50])

        title(['IC ', num2str(i)]);
        ax = gca;
        set(ax, 'ButtonDownFcn', ...
            @(src, event)showDetails({XP1, XP3, XP6}, ...
            {mean_signalP1, mean_signalP3, mean_signalP6}, ...
            {std_sem_signalP1, std_sem_signalP3, std_sem_signalP6}, ...
            i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
    end
    
    % Create a legend for the entire tiled layout
    lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south'; % Position the legend at the bottom

    %% Callback function to display detailed plot
    function showDetails(X, mean_signals, std_signals, ch, title_n, XLabel_n, YLabel_n, XLim_n, All_colours)
        figure; hold on;
        alpha = 0.3;
    
        fill([X{1,1} fliplr(X{1,1})], [mean_signals{1,1}(ch, :) + std_signals{1,1}(ch, :), ...
            fliplr(mean_signals{1,1}(ch, :) - std_signals{1,1}(ch, :))], ...
            All_colours.light_blue, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        fill([X{1,2} fliplr(X{1,2})], [mean_signals{1,2}(ch, :) + std_signals{1,2}(ch, :), ...
            fliplr(mean_signals{1,2}(ch, :) - std_signals{1,2}(ch, :))], ...
            All_colours.light_orange, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        fill([X{1,3} fliplr(X{1,3})], [mean_signals{1,3}(ch, :) + std_signals{1,3}(ch, :), ...
            fliplr(mean_signals{1,3}(ch, :) - std_signals{1,3}(ch, :))], ...
            All_colours.light_green, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    
    
        h1 = plot(X{1,1}, mean_signals{1,1}(ch, :), 'Color', All_colours.dark_blue, 'LineWidth', 2); 
    
        h2 = plot(X{1,2}, mean_signals{1,2}(ch, :), 'Color', All_colours.dark_orange, 'LineWidth', 2); 
        
        h3 = plot(X{1,3}, mean_signals{1,3}(ch, :), 'Color', All_colours.dark_green, 'LineWidth', 2); 
        
        set(gca, 'XLim', XLim_n)
    
        title(title_n)
        xlabel(XLabel_n);
        ylabel(YLabel_n);
        grid on;
    
        legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal')
    
        zoom on
    end

end