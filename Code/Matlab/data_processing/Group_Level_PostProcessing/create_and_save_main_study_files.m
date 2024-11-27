function create_and_save_main_study_files(subject_list, data_path, processed_data_path, ALLEEG)
    
    
    % Load datasets into ALLEEG
    for i = 1:length(subject_list)
        file_name = ['sub-', subject_list{i}, '_cleaned_with_ICA.set'];
        dataset_path = fullfile([processed_data_path, 'sub-', ...
            subject_list{i}, filesep]);
        EEG = pop_loadset('filename', file_name, 'filepath', dataset_path);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
    end
    
    
    %% Create STUDY main_study_all_ICs_RV-15
    study_name = 'main_study_all_ICs_RV-15.study';
    study_folder = [data_path, '7_STUDY'];
    
    % Create STUDY from loaded datasets using commands (RV <= 15%)
    commands = cell(1, length(subject_list));
    for i = 1:length(subject_list)
        commands{i} = {'index', i, 'subject', subject_list{i}, ...
            'inbrain', 'on', 'dipselect', 0.15};
    end
    
    [STUDY, ALLEEG] = std_editset([], ALLEEG, ...
        'name', study_name, ...
        'filepath', study_folder, ...
        'commands', commands, ...
        'updatedat', 'on', 'savedat', 'off');
    
    %% Save STUDY (main_study_all_ICs_RV-15)
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', study_name, 'filepath', study_folder);
    
    
    %% Save STUDY with brain ICs (RV <= 15%)
    % Identify the brain ICs manually (it was done before by manually checking
    % all of the ICs and selecting potential brain components and reject other
    % non-brain sources.
    Brain_ICs = cortical_ICs_indentifier(subject_list, ALLEEG);
    
    for i = 1:length(subject_list)
        EEG = ALLEEG(i); % Get subject dataset
        all_ICs = 1:size(EEG.icaweights, 1);
        retain_ICs = union(Brain_ICs{i, 2}, Brain_ICs{i, 3});
        EEG = pop_subcomp(EEG, setdiff(all_ICs, retain_ICs), 0); % Remove other ICs
        ALLEEG(i) = EEG; % Store back the modified dataset
    end
    
    % Create STUDY
    study_name = 'main_study_brain_ICs_RV-15.study';
    
    [STUDY, ALLEEG] = std_editset([], ALLEEG, ...
        'name', study_name, ...
        'filepath', study_folder, ...
        'commands', commands, ...
        'updatedat', 'on', 'savedat', 'off');
    
    % Save STUDY
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', study_name, 'filepath', study_folder);

end