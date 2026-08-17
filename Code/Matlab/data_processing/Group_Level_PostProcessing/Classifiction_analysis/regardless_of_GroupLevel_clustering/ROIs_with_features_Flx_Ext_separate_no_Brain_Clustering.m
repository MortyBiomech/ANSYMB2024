function ROIs = ROIs_with_features_Flx_Ext_separate_no_Brain_Clustering(epoch_type, ...
    features_from_epochs, data_path, per_trial_or_all_epochs, pressure_score, ...
    grouplevel_postprocessing_path)

    %% load data containing the finalized Brain ICs per subject
    subject_list = 5:18;
    ICs = load(fullfile(grouplevel_postprocessing_path, ...
        'Brain_PotentialBrain_AcceptedPotentialBrain.mat'));
    name = fieldnames(ICs);
    ICs = ICs.(name{1});

    Brain_ICs = cell(size(subject_list));
    for i = 1:numel(subject_list)
        Brain_ICs{i} = [cell2mat(ICs(i, 1)); cell2mat(ICs(i, 3))];
    end
    
    
    %% Create meta STUDY files
    all_STUDY_names = cell(size(subject_list));
    for i = 1:length(all_STUDY_names)
        all_STUDY_names{i} = ['TotalBrain_sub_', num2str(subject_list(i))];
    end
    
    all_STUDY_files = cell(size(all_STUDY_names));
    
    for i = 1:length(all_STUDY_files)
        all_STUDY_files{1, i} = Brain_ICs{1, i};
    end
    
    ROIs = struct();
    for i = 1:length(all_STUDY_names)
        N = length(all_STUDY_files{1, i});
        ROIs.(all_STUDY_names{1, i}) = cell(N, 3);
    end
    
    for i = 1:length(all_STUDY_names)
        sets = length(all_STUDY_files{1, i});
        comps = all_STUDY_files{1, i};
        ROIs.(all_STUDY_names{1, i})(1:sets, 1) = num2cell(subject_list(i));
        ROIs.(all_STUDY_names{1, i})(1:sets, 2) = num2cell(comps);
    end
    
    
    %% Load epoched datasets to fill the metafile
    for i = 1:length(subject_list)
        
        %%
        disp(['Subject ', num2str(subject_list(i))])
        filepath = [data_path, '6_Trials_Info_and_Epoched_data\sub-', ...
            num2str(subject_list(i))];
        data = load(fullfile(filepath, epoch_type));
        name = fieldnames(data);
        data = data.(name{1});
        
        Trials_Info = load(fullfile(filepath, 'Trials_Info.mat'));
        name = fieldnames(Trials_Info);
        Trials_Info = Trials_Info.(name{1});

        

        if strcmp(pressure_score, 'pressure')
            % Find all conditions indices in trials 
            % P1, P3, and P6
            condition_indices = condition_indices_identifier(...
                Trials_Info, subject_list(i));
        elseif strcmp(pressure_score, 'score')
            % Find all conditions indices in trials 
            % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
            condition_indices = condition_indices_identifier_ScoreBased(...
                Trials_Info, subject_list(i));
        end
    

        % frequencies of the Power Spectral Density
        [data, frequencies] = calculate_PSD_Flx_Ext_separate(...
            data, Trials_Info, subject_list(i));
        for j = 1:length(all_STUDY_files)
            
            subject_sets_indx = ...
                find(cell2mat(ROIs.(all_STUDY_names{1, j})(:,1)) == subject_list(i));
            subject_comps = cell2mat(ROIs.(all_STUDY_names{1, j})(subject_sets_indx, 2));
            if ~isempty(subject_comps)
                
                for k = 1:length(subject_sets_indx)
                    IC = subject_comps(k);
                    % [ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3), ...
                    %  ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 4)] = ...
                    %     RMS_features_generator_Flx_Ext_separate(...
                    %         condition_indices, data, IC, frequencies, ...
                    %         per_trial_or_all_epochs, pressure_score);

                    ROIs.(all_STUDY_names{1, j})(subject_sets_indx(k), 3) = ...
                        PSD_integral_feature_generator_Flx_Ext_separate(...
                            condition_indices, data, IC, frequencies, ...
                            per_trial_or_all_epochs, pressure_score);
                end
    
            end
    
        end


        clear data
        clear Trials_Info

    end


    %% Save ROI (check if there is a older version, save a new one)
    classisfication_path = [data_path, '8_classification\'];
    ROIs_features_path = [classisfication_path, ...
        'ROIs_features\regardless_of_clustering'];
    
    % Define the base name and extension
    baseName = 'ROIs';
    suffix = features_from_epochs;
    folder = ROIs_features_path; 
    extension = ['_', per_trial_or_all_epochs, '_', ...
        pressure_score, '_Based','_Flx_Ext_separate_NoBrainClustering_PSD_integral.mat']; 
    
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
