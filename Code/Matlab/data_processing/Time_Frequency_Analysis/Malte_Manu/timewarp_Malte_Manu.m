function timewarp_Malte_Manu(outputDir)
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    
    mainDir = pwd; % Setzt mainDir auf das aktuelle Arbeitsverzeichnis

    analysisDir = fullfile(mainDir, 'analysis');
    matFiles = dir(fullfile(analysisDir, '**', '*.mat'));
    % Cluster-Variablen initialisieren
    numClusters = min(length(matFiles), 8); % Maximal 8 Cluster
    
    lengthFlexion = [];
    lengthExtension = [];
    for i = 1:numClusters
        filePath = fullfile(matFiles(i).folder, matFiles(i).name);
        dataStruct = load(filePath);
        varName = fieldnames(dataStruct);
        cellArray = dataStruct.(varName{1});
        dataFile = cellArray(:,3);

        for k = 1:length(dataFile)
            n_trials = size(dataFile{k}, 2);

            for trialIdx = 1:n_trials
                % Extrahiere Event-Zeitpunkte
                events = dataFile{k}{3, trialIdx};
                flexion_start_indx = events.flexion_start_indx;
                extension_start_indx = events.extension_start_indx;
                extension_end_indx = events.extension_end_indx;
    
                % Bestimme die minimale Länge
                minLength = min([length(flexion_start_indx), length(extension_start_indx), length(extension_end_indx)]);
                
                if minLength == 0
                    continue;
                end

                % Kürze alle Arrays auf die minimale Länge
                flexion_start_indx = flexion_start_indx(1:minLength);
                extension_start_indx = extension_start_indx(1:minLength);
                extension_end_indx = extension_end_indx(1:minLength);
    
                for l = 1:minLength
                    origLengthFlexion = extension_start_indx(l) - flexion_start_indx(l);
                    origLengthExtension = extension_end_indx(l) - extension_start_indx(l) + 1;
    
                    lengthFlexion = [lengthFlexion origLengthFlexion];
                    lengthExtension = [lengthExtension origLengthExtension];
                end                
            end
        end
    end

    medianFlexion = median(lengthFlexion);
    medianExtension = median(lengthExtension);

    disp("Median Flexion = " + medianFlexion)
    disp("Median Extension = " + medianExtension)

    for i = 1:numClusters
        filePath = fullfile(matFiles(i).folder, matFiles(i).name);
        dataStruct = load(filePath);
        varName = fieldnames(dataStruct);
        cellArray = dataStruct.(varName{1});
        dataFile = cellArray(:,3);

        for k = 1:length(dataFile)
            n_trials = size(dataFile{k}, 2);
            allTrials = cell(1, n_trials); % Cell-Array für alle Trials

            for trialIdx = 1:n_trials

                disp("Cluster: " + i + " -> #IC " + k + " -> trial: " + trialIdx)
    
                % Extrahiere Event-Zeitpunkte
                events = dataFile{k}{3, trialIdx};
                flexion_start_indx = events.flexion_start_indx;
                extension_start_indx = events.extension_start_indx;
                extension_end_indx = events.extension_end_indx;
    
                % Bestimme die minimale Länge
                minLength = min([length(flexion_start_indx), length(extension_start_indx), length(extension_end_indx)]);
                
                if minLength == 0
                    continue;
                end

                % Kürze alle Arrays auf die minimale Länge
                flexion_start_indx = flexion_start_indx(1:minLength);
                extension_start_indx = extension_start_indx(1:minLength);
                extension_end_indx = extension_end_indx(1:minLength);

                allEpochs = cell(1, minLength); % Cell-Array für alle Epochen
                
                for epoch = 1:minLength
                    try
                        epochData = dataFile{k}{1, trialIdx}{epoch};
                    catch
                        continue;
                    end
        
                    origLengthFlexion = extension_start_indx(epoch) - flexion_start_indx(epoch);
                    origLengthExtension = extension_end_indx(epoch) - extension_start_indx(epoch) + 1;
                    origLengthFull = extension_end_indx(epoch) - flexion_start_indx(epoch);
                    
                    if origLengthFlexion <= 1 || origLengthExtension <= 1 || origLengthFull <= 1
                        continue;
                    end
    
                    flexionData = epochData(:, 1:origLengthFlexion);
                    extensionData = epochData(:, origLengthFlexion:origLengthFull);
        
                    if origLengthFlexion < medianFlexion || origLengthFlexion > medianFlexion
                        tempDataFlexion = do_timewarp(flexionData, origLengthFlexion, medianFlexion);
                    else
                        tempDataFlexion = flexionData;
                    end
        
                    if origLengthExtension < medianExtension || origLengthExtension > medianExtension
                        tempDataExtension = do_timewarp(extensionData, origLengthExtension, medianExtension);
                    else
                        tempDataExtension = extensionData;
                    end
                    allEpochs{epoch} = [tempDataFlexion tempDataExtension];
                    dataFile{k}{1, trialIdx}{epoch} = [tempDataFlexion tempDataExtension];
                end
                allTrials{trialIdx} = allEpochs;
            end

        end
        
        % Speichern der neuen Datei für alle Epochen des Subjekts
        validVarName = matlab.lang.makeValidName(varName{1});
        cellArray(:,3) = dataFile;  % Aktualisiere die veränderten Daten
        eval([validVarName, ' = cellArray;']); % Erstelle die Variable mit dem Namen der Datei
        
        [~, fileName, ~] = fileparts(matFiles(i).name); % Hole den Originaldateinamen ohne .mat-Endung
        saveFilePath = fullfile(outputDir, [fileName, '.mat']); % Gleicher Name für neue Datei
        
        % Erstelle das Verzeichnis, falls es noch nicht existiert
        if ~exist(fileparts(saveFilePath), 'dir')
            mkdir(fileparts(saveFilePath));
        end
        
        % Speichere die Variable unter dem gleichen Namen wie die Datei
        save(saveFilePath, validVarName, '-v7.3');
        fprintf('Verarbeitet: %s -> %s\n', varName{1}, saveFilePath);


    end
end


function timewarpedData = do_timewarp(data, origLength, newLength)

    [n_ICs, ~] = size(data);

    origPoints = linspace(1, origLength, origLength);
    newPoints = linspace(1, origLength, newLength);

    timewarpedData = zeros(n_ICs, newLength);

    for ch = 1:n_ICs
        timewarpedData(ch, :) = interp1(origPoints, data(ch, :), newPoints, "pchip");
    end
    
end
