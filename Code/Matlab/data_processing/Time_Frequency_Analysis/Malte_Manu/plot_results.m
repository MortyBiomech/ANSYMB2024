function plot_results(outputDir)
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    mainDir = pwd; 
    timewarpedDir = fullfile(mainDir, 'timewarped');
    matFiles = dir(fullfile(timewarpedDir, '**', '*.mat'));

    numClusters = min(length(matFiles), 8); 

    for i = 1:numClusters
        filePath = fullfile(matFiles(i).folder, matFiles(i).name);
        dataStruct = load(filePath);
        varName = fieldnames(dataStruct);
        cellArray = dataStruct.(varName{1});
        dataFile = cellArray(:,3);

        [~, clusterName, ~] = fileparts(matFiles(i).name);

        Pressure1 = [];
        Pressure3 = [];
        Pressure6 = [];

        for k = 1:length(dataFile)
            n_trials = size(dataFile{k}, 2);

            for trialIdx = 1:n_trials
                
                if ~strcmpi(dataFile{k}{4, trialIdx}, 'Experiment')
                    continue;
                end

                n_epochs = length(dataFile{k}{1, trialIdx});
                if n_epochs == 0
                    continue;
                end

                pressureValue = dataFile{k}{2, trialIdx};

                for n_epoch = 1:n_epochs
                    epochData = dataFile{k}{1, trialIdx}{n_epoch}; 
                    
                    if size(epochData, 1) ~= 50
                        warning('Skipping epochData with incorrect height: %dx%d', size(epochData,1), size(epochData,2));
                        continue;
                    end

                    if isempty(Pressure1) && pressureValue == 1
                        Pressure1 = epochData;
                        continue;
                    elseif isempty(Pressure3) && pressureValue == 3
                        Pressure3 = epochData;
                        continue;
                    elseif isempty(Pressure6) && pressureValue == 6
                        Pressure6 = epochData;
                        continue;
                    end

                    if pressureValue == 1 && size(epochData, 2) == size(Pressure1, 2)
                        Pressure1 = cat(3, Pressure1, epochData);
                    elseif pressureValue == 3 && size(epochData, 2) == size(Pressure3, 2)
                        Pressure3 = cat(3, Pressure3, epochData);
                    elseif pressureValue == 6 && size(epochData, 2) == size(Pressure6, 2)
                        Pressure6 = cat(3, Pressure6, epochData);
                    else
                        warning('Skipping epochData with incorrect width: %dx%d', size(epochData,1), size(epochData,2));
                    end
                end
            end
        end

        compute_dB = @(x) 2 * ((10 * log10(x) - min(10 * log10(x(:)))) / (max(10 * log10(x(:))) - min(10 * log10(x(:))))) - 1;

        plot_and_save(Pressure1, clusterName, 1, compute_dB, outputDir);
        plot_and_save(Pressure3, clusterName, 3, compute_dB, outputDir);
        plot_and_save(Pressure6, clusterName, 6, compute_dB, outputDir);
    end
end

function plot_and_save(Pressure, clusterName, pressureValue, compute_dB, outputDir)
    if ~isempty(Pressure)
        varName = sprintf('%s__Pressure%d', clusterName, pressureValue);
        eval([varName ' = Pressure;']);
        save(fullfile(outputDir, [varName, '.mat']), varName);
        disp("Saved Pressure-3D-Matrix: " + varName)
        
        Pressure_mean = mean(Pressure, 3);
        varNameMean = sprintf('%s__Pressure%d_mean', clusterName, pressureValue);
        eval([varNameMean ' = Pressure_mean;']);
        save(fullfile(outputDir, [varNameMean, '.mat']), varNameMean);
        disp("Saved Pressure-Mean-2D-Matrix: " + varNameMean)
        
        Pressure_dB = compute_dB(Pressure_mean);
        varName_dB = sprintf('%s__Pressure%d_dB', clusterName, pressureValue);
        eval([varName_dB ' = Pressure_dB;']);
        save(fullfile(outputDir, [varName_dB, '.mat']), varName_dB);
        disp("Saved Pressure-Mean-dB-Normalized-2D-Matrix: " + varName_dB)
        
        figure;
        imagesc(Pressure_dB);
        colormap jet;
        colorbar;
        title(sprintf('%s - Pressure %d dB Normalized', clusterName, pressureValue));
        xlabel('Time Points');
        ylabel('Channels');
        saveas(gcf, fullfile(outputDir, sprintf('%s__Pressure%d_dB.png', clusterName, pressureValue)));
        disp("Saved Plot: " + sprintf('%s__Pressure%d_dB.png', clusterName, pressureValue));
    end
end
