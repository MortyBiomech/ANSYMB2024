function morlet_wavelet()
    % Festgelegte Probanden
    subjects = {'sub_11', 'sub_15', 'sub_16', 'sub_18'};
    
    % Arbeitsverzeichnis
    base_folder = pwd;  % Aktuelles Verzeichnis, in dem sich das Skript befindet
    analysis_folder = fullfile(base_folder, 'analysis');
    
    % Parameter des Morlet Wavelets
    disp('Initialisiere Parameter für das Morlet Wavelet...');
    Fs = 500; % Sampling-Rate in Hz
    t = -2:1/Fs:2; 
    f_min = 1; % Minimale Frequenz in Hz
    f_max = 50; % Maximale Frequenz in Hz
    num_frequencies = 250; % Anzahl der Frequenzbänder
    frequencies = logspace(log10(f_min), log10(f_max), num_frequencies); % Logarithmisch verteilte Frequenzen
    width = linspace(4, 10, num_frequencies); % Breite des Morlet Wavelets (Zyklen)
    disp('Parameter initialisiert.');

    % Erstelle den Ordner 'analysis', falls nicht vorhanden
    if ~exist(analysis_folder, 'dir')
        mkdir(analysis_folder);
        disp('Analyse-Ordner wurde erstellt.');
    else
        disp('Analyse-Ordner existiert bereits.');
    end
    
    % Lade die ROI-Datei
    disp('Lade ROI-Datei...');
    roi_file = fullfile(base_folder, 'ROIs_0_FlextoFlex.mat');
    if ~isfile(roi_file)
        error('Die Datei ROIs_0_FlextoFlex.mat konnte nicht gefunden werden.');
    end
    roi_data = load(roi_file);
    disp('ROI-Datei erfolgreich geladen.');

    % Überprüfen, ob die "ROIs"-Daten vorhanden sind
    if ~isfield(roi_data, 'ROIs')
        error('Die Datei ROIs_0_FlextoFlex.mat enthält keine Daten mit dem Namen "ROIs".');
    end
    rois = roi_data.ROIs;

    % Lade alle Probandendaten vorab
    disp('Lade Probandendaten...');
    subject_data = struct();
    trial_info_data = struct();
    for subj_idx = 1:numel(subjects)
        subject = subjects{subj_idx};
        subject_folder = fullfile(base_folder, subject);
        epochs_file = fullfile(subject_folder, 'Epochs_FlextoFlex_based.mat');
        trials_info_file = fullfile(subject_folder, 'Trials_Info.mat');
        
        if ~isfile(epochs_file)
            warning(['Die Datei ', epochs_file, ' konnte nicht gefunden werden. Überspringe ', subject, '.']);
            continue;
        end

        % Lade die Daten und speichere sie in der Struktur
        temp_data = load(epochs_file, 'Epochs_FlextoFlex_based');
        subject_data.(subject) = temp_data.Epochs_FlextoFlex_based;
        temp2_data = load(trials_info_file, 'Trials_Info');
        trial_info_data.(subject) = temp2_data.Trials_Info;
        disp(['Probandendaten für ', subject, ' erfolgreich geladen.']);
    end

    % Gehirnregionen
    regions = fieldnames(rois);

    % Iteriere über die Gehirnregionen
    for region_idx = 1:numel(regions)
        region_name = regions{region_idx};
        disp(['--- Verarbeite Region: ', region_name, ' ---']);
        
        % Initialisiere die Ergebnistabelle als nx3 Cell-Array
        result_table = {};

        % Iteriere über die Probanden
        for subj_idx = 1:numel(subjects)
            subject = subjects{subj_idx};
            subject_num = subject(5:end);
            disp(['   -> Verarbeite Proband: ', subject, ' (Nr. ', subject_num, ')']);
            
            % Überprüfe, ob die Probandendaten existieren
            if ~isfield(subject_data, subject)
                warning(['Keine geladenen Daten für ', subject, ' gefunden. Überspringe.']);
                continue;
            end
            n_trials = numel(subject_data.(subject));
            disp(['   -> Anzahl der Trials: ', num2str(n_trials)]);

            % Finde relevante ICs für diesen Probanden und die aktuelle Region
            roi_table = rois.(region_name);
            roi_subjects = roi_table(:, 1);
            if isnumeric(roi_subjects{1})
                roi_subjects = cellfun(@num2str, roi_subjects, 'UniformOutput', false);
            end
            roi_subjects = strtrim(roi_subjects);
            relevant_rows = strcmp(roi_subjects, subject_num);
            relevant_ICs = cell2mat(roi_table(relevant_rows, 2));

            if isempty(relevant_ICs)
                warning(['Keine relevanten ICs für ', subject, ' (Nr. ', subject_num, ') in Region ', region_name, ' gefunden.']);
                continue;
            end
            disp(['   -> Anzahl relevanter ICs: ', num2str(numel(relevant_ICs))]);
            
            for IC_index = 1:numel(relevant_ICs)
                
                % Extrahiere Daten für alle Trials und relevante ICs
                trial_data = cell(4, n_trials);
                
                for trial_idx = 1:n_trials
                    trial = subject_data.(subject){trial_idx}; 
                    trial_sources = trial.EEG_stream.Preprocessed.Sources;
    
                    general = subject_data.(subject){trial_idx};
                    pressureInfo = general.General.Pressure;
                    description = general.General.Description;

                    trial_info_subject = trial_info_data.(subject){trial_idx};
                    trials_info_content = trial_info_subject.Events.EEG_stream.Preprocessed;

                    n_epochs = numel(trial_sources);
                    TF_data_trial = cell(1, n_epochs);
    
                    for epoch_idx = 1:n_epochs
                        signal = trial_sources{epoch_idx};
                        relevant_signal = signal(relevant_ICs(IC_index), :);
                        power_matrix = zeros(num_frequencies, size(relevant_signal, 2));
    
                        for fi = 1:num_frequencies
                            freq = frequencies(fi);
                            w = width(fi);
                            sigma_t = w / (2 * pi * freq);
                            wavelet = exp(2 * 1i * pi * freq * t) .* exp(-t.^2 / (2 * sigma_t^2));
                            wavelet = wavelet / sqrt(sum(abs(wavelet).^2));
                            convolution = conv(relevant_signal, wavelet, 'same');
                            power_matrix(fi, :) = abs(convolution).^2;
                        end
                        TF_data_trial{epoch_idx} = power_matrix;
                    end
                    trial_data{1, trial_idx} = TF_data_trial;
                    trial_data{2, trial_idx} = pressureInfo;
                    trial_data{3, trial_idx} = trials_info_content;
                    trial_data{4, trial_idx} = description;
    
                    % Fortschrittsanzeige pro Trial
                    disp(['      -> Trial ', num2str(trial_idx), ' von ', num2str(n_trials), ' abgeschlossen.']);
                end
    
                % Aktualisiere die Ergebnistabelle (nx3)
                result_table = [result_table; {subject_num, relevant_ICs(IC_index), trial_data}];
            end
        end
        % Speichere die Regionsergebnisse
        save_path = fullfile(analysis_folder, [region_name, '.mat']);
        region_variable_name = matlab.lang.makeValidName(region_name); % Gültiger MATLAB-Variablenname
        eval([region_variable_name, ' = result_table;']); % Weist die Tabelle der Variablen mit Regionsnamen zu
        save(save_path, region_variable_name, '-v7.3'); % Speichert die Variable mit ihrem Regionsnamen
        disp(['--- Region ', region_name, ' abgeschlossen. Ergebnisse gespeichert. ---']);

    end
end
