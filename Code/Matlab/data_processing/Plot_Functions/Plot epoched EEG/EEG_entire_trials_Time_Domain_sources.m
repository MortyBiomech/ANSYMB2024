function EEG_entire_trials_Time_Domain_sources(data, std_sem, P1, P3, P6)

    All_colours = struct('dark_blue', [0, 0.4470, 0.7410], 'light_blue', [0.3010, 0.7450, 0.9330], ...
        'dark_orange', [0.8500, 0.3250, 0.0980], 'light_orange', [0.9290, 0.6940, 0.1250], ...
        'dark_green', [0.4660, 0.6740, 0.1880], 'light_green', [0.5960, 0.8740, 0.5410]);


    %% EEG_trials - Time-Domain: channels
    % making 3d matrix for each pressure condiotion
    % P1
    numChannels = size(data{1, P1.trials(1)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized, 1);
    numPointsP1 = size(data{1, P1.trials(1)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized, 2);
    numEpochsP1 = length(P1.trials);
    signal_P1_3d = zeros(numChannels, numPointsP1, numEpochsP1);
    for i = 1:numEpochsP1
        tempdata = data{1, P1.trials(i)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized;
        signal_P1_3d(:,:, i) = tempdata;
    end
    
    % P3
    numPointsP3 = size(data{1, P3.trials(1)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized, 2);
    numEpochsP3 = length(P3.trials);
    signal_P3_3d = zeros(numChannels, numPointsP3, numEpochsP3);
    for i = 1:numEpochsP3
        tempdata = data{1, P3.trials(i)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized;
        signal_P3_3d(:,:, i) = tempdata;
    end
    
    % P6
    numPointsP6 = size(data{1, P6.trials(1)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized, 2);
    numEpochsP6 = length(P6.trials);
    signal_P6_3d = zeros(numChannels, numPointsP6, numEpochsP6);
    for i = 1:numEpochsP6
        tempdata = data{1, P6.trials(i)}.EEG_stream.Preprocessed.Time_Domain.Sources.Length_Normalized;
        signal_P6_3d(:,:, i) = tempdata;
    end
    
    
    % perecent of epoch for x-axis
    XP1 = linspace(0, 100, numPointsP1);
    XP3 = linspace(0, 100, numPointsP3);
    XP6 = linspace(0, 100, numPointsP6);
    
    mean_signalP1 = mean(signal_P1_3d, 3);
    mean_signalP3 = mean(signal_P3_3d, 3);
    mean_signalP6 = mean(signal_P6_3d, 3);
    
    std_signalP1 = std(signal_P1_3d, 0, 3);
    std_signalP3 = std(signal_P3_3d, 0, 3);
    std_signalP6 = std(signal_P6_3d, 0, 3);

    sem_signalP1 = std_signalP1 ./ sqrt(numEpochsP1);
    sem_signalP3 = std_signalP3 ./ sqrt(numEpochsP3);
    sem_signalP6 = std_signalP6 ./ sqrt(numEpochsP6);

    if strcmp(std_sem, 'STD')
        stats_signalP1 = std_signalP1;
        stats_signalP3 = std_signalP3;
        stats_signalP6 = std_signalP6;
    elseif strcmp(std_sem, 'SEM')
        stats_signalP1 = sem_signalP1;
        stats_signalP3 = sem_signalP3;
        stats_signalP6 = sem_signalP6;
    else
        disp('second input must be either "std" or "sem".')
    end
    
    %% plot the results
    h1 = []; h2 = []; h3 = [];
    XLabel_n = 'Cycle [%]';
    YLabel_n = 'Amplitude [\muV]';
    XLim_n = [0 100];
    figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
    tiledlayout(8,8)
    sgtitle(sprintf('Time-Domain (Entire Trial, Mean $\\pm$ %s): %d trials P1, %d trials P3, %d trials P6', ...
        std_sem, numEpochsP1, numEpochsP3, numEpochsP6), 'Interpreter', 'latex');
    for i = 1:numChannels
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
        
        title(['IC ', num2str(i)]);
        ax = gca;
        set(ax, 'ButtonDownFcn', ...
            @(src, event)showDetails({XP1, XP3, XP6}, ...
            {mean_signalP1, mean_signalP3, mean_signalP6}, ...
            {stats_signalP1, stats_signalP3, stats_signalP6}, ...
            i, ['IC ', num2str(i)], XLabel_n, YLabel_n, XLim_n, All_colours));
    end
    
    % Create a legend for the entire tiled layout
    lgd = legend([h1, h2, h3], {'P1', 'P3', 'P6'}, 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south'; % Position the legend at the bottom

    %% Callback function to display detailed plot
    function showDetails(X, mean_signals, stats_signals, ch, title_n, XLabel_n, YLabel_n, XLim_n, All_colours)
        figure; hold on;
        alpha = 0.3;
        means_alpha = 0.4;
    
        fill([X{1,1} fliplr(X{1,1})], [mean_signals{1,1}(ch, :) + stats_signals{1,1}(ch, :), ...
            fliplr(mean_signals{1,1}(ch, :) - stats_signals{1,1}(ch, :))], ...
            All_colours.light_blue, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        fill([X{1,2} fliplr(X{1,2})], [mean_signals{1,2}(ch, :) + stats_signals{1,2}(ch, :), ...
            fliplr(mean_signals{1,2}(ch, :) - stats_signals{1,2}(ch, :))], ...
            All_colours.light_orange, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        fill([X{1,3} fliplr(X{1,3})], [mean_signals{1,3}(ch, :) + stats_signals{1,3}(ch, :), ...
            fliplr(mean_signals{1,3}(ch, :) - stats_signals{1,3}(ch, :))], ...
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