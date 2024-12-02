function create_and_save_main_study_files(subject_list, data_path, processed_data_path, ALLEEG)
    
    
    % Load datasets into ALLEEG
    for i = 1:length(subject_list)
        file_name = ['sub-', num2str(subject_list(i)), '_cleaned_with_ICA.set'];
        dataset_path = fullfile([processed_data_path, 'sub-', ...
            num2str(subject_list(i)), filesep]);
        EEG = pop_loadset('filename', file_name, 'filepath', dataset_path);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, i);
    end
    
    
    %% Create STUDY main_study_all_ICs_RV-15
    study_name = 'main_study_all_ICs_RV-15.study';
    study_folder = [data_path, '7_STUDY'];
    
    % Create STUDY from loaded datasets using commands (RV <= 15%)
    commands = cell(1, length(subject_list));
    for i = 1:length(subject_list)
        commands{i} = {'index', i, 'subject', num2str(subject_list(i)), ...
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

    % This function was run and brain ICs for subjects 5 to 10 were found.
    % Brain_ICs = cortical_ICs_indentifier(subject_list, ALLEEG);
    load Brain_PotentialBrain_AcceptedPotentialBrain.mat ICs
    
    subjects_ids_in_STUDY = unique(STUDY.cluster.sets, 'stable');

    for i = 1:length(subjects_ids_in_STUDY)
        
        subjects_ids_index = find(STUDY.cluster.sets == ...
            subjects_ids_in_STUDY(i));
        
        all_RV_15_ICs = STUDY.cluster.comps(subjects_ids_index);
        all_potential_brain_ICs = union(ICs{subjects_ids_in_STUDY(i), 1}, ...
            ICs{subjects_ids_in_STUDY(i), 3});
        
        [~, common_elements_indxs] = intersect(all_RV_15_ICs, all_potential_brain_ICs);
        indices_to_keep = subjects_ids_index(common_elements_indxs);
        
        STUDY.cluster.sets(subjects_ids_index(~ismember(subjects_ids_index, indices_to_keep))) = [];
        STUDY.cluster.comps(subjects_ids_index(~ismember(subjects_ids_index, indices_to_keep))) = [];

        STUDY.datasetinfo(subjects_ids_in_STUDY(i)).comps...
            (~ismember(1:length(all_RV_15_ICs), common_elements_indxs)) = [];
    end

    
    % Save new STUDY with just potential Brain ICs
    file_name = 'main_study_potential_brain_ICs_RV-15.study';

    STUDY.name = file_name;
    
    % Save STUDY
    [STUDY, ALLEEG] = pop_savestudy(STUDY, ALLEEG, 'filename', file_name, 'filepath', study_folder);

end