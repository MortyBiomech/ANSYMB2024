function ROIs = ROIs_with_features(all_STUDY_names, all_STUDY_files, subject_list, ...
    epoch_type, features_from_epochs, data_path, main_project_folder)

    ROIs = struct();
    for i = 1:length(all_STUDY_names)
        best_cluster = ...
            all_STUDY_files{1, i}.etc.bemobil.clustering.cluster_ROI_index; 
        N = length(all_STUDY_files{1, i}.cluster(1, best_cluster).sets);
        ROIs.(all_STUDY_names{1, i}) = cell(N, 4);
    end
    
    
    for i = 1:length(all_STUDY_names)
        best_cluster = ...
            all_STUDY_files{1, i}.etc.bemobil.clustering.cluster_ROI_index; 
        sets = all_STUDY_files{1, i}.cluster(1, best_cluster).sets;
        comps = all_STUDY_files{1, i}.cluster(1, best_cluster).comps;
    
        ROIs.(all_STUDY_names{1, i})(1:length(sets), 1) = num2cell(subject_list(sets));
        ROIs.(all_STUDY_names{1, i})(1:length(sets), 2) = num2cell(comps);
        
    end
    
    
    %% Load epoched datasets to fill the metafile
    for i = 1:length(subject_list)
        
        disp(['Subject ', num2str(subject_list(i))])
        cd([data_path, '6_Trials_Info_and_Epoched_data\sub-', ...
            num2str(subject_list(i))]);
        
        data = load(epoch_type);
        [~, baseName, ~] = fileparts(epoch_type);
        data = data.(baseName);
        
        Trials_Info = load('Trials_Info.mat');
        [~, baseName, ~] = fileparts('Trials_Info.mat');
        Trials_Info = Trials_Info.(baseName);

        cd(fullfile(main_project_folder, 'Code\Matlab\data_processing\', ...
            'Group_Level_PostProcessing\Classifiction_analysis'))
    

        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject_list(i));
    

        if subject_list(i) < 10 % data structure for subject 10 and above was saved differently
            
            % frequencies of the Power Spectral Density
            frequencies = data{1, 1}.EEG_stream.Preprocessed.Freq_Domain.Freqs;
            
            for j = 1:length(all_STUDY_files)
                
                subject_sets_indx = ...
                    find(cell2mat(ROIs.(all_STUDY_names{1, j})(:,1)) == subject_list(i));
                subject_comps = cell2mat(ROIs.(all_STUDY_names{1, j})(subject_sets_indx, 2));
                if ~isempty(subject_comps)
                    
                    for k = 1:length(subject_sets_indx)
                        IC = subject_comps(k);
                        [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
                            ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
                            RMS_features_generator(condition_indices, data, IC, frequencies);
                    end
        
                end
        
            end

        else

            % frequencies of the Power Spectral Density
            [data, frequencies] = calculate_PSD(data);
            for j = 1:length(all_STUDY_files)
                
                subject_sets_indx = ...
                    find(cell2mat(ROIs.(all_STUDY_names{1, j})(:,1)) == subject_list(i));
                subject_comps = cell2mat(ROIs.(all_STUDY_names{1, j})(subject_sets_indx, 2));
                if ~isempty(subject_comps)
                    
                    for k = 1:length(subject_sets_indx)
                        IC = subject_comps(k);
                        [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
                            ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
                            RMS_features_generator(condition_indices, data, IC, frequencies);
                    end
        
                end
        
            end

        end

        clear data
        clear Trials_Info

    end


    %% Save ROI (check if there is a older version, save a new one)
    classisfication_path = [data_path, '8_classification\'];
    ROIs_features_path = [classisfication_path, 'ROIs_features'];
    
    % Define the base name and extension
    baseName = 'ROIs';
    suffix = features_from_epochs;
    folder = ROIs_features_path; 
    extension = '.mat'; 
    
    % Initialize version number
    version = 0;
    fileName = sprintf('%s_%d_%s%s', baseName, version, suffix, extension);
    fullFilePath = fullfile(folder, fileName);
    
    % Increment version until a unique file name is found
    while exist(fullFilePath, 'file') == 2
        version = version + 1;
        fileName = sprintf('%s_%d_%s%s', baseName, version, suffix, extension);
        fullFilePath = fullfile(folder, fileName);
    end
    
    % Save your file
    save(fullFilePath, 'ROIs'); 
    disp(['File saved as: ' fileName]);

end
make_timewarp