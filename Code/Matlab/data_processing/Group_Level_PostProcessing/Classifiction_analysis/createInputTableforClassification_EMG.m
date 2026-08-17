function subjectTable = ...
    createInputTableforClassification_EMG(subject, EMG_features_path, ...
    pressure_score, epoched_data_path, classes)


    filename = ['sub-', num2str(subject), '\sub-', num2str(subject), ...
        '_structured_EMG_data.mat'];
    data = load(fullfile(EMG_features_path, filename));
    name = fieldnames(data);
    data = data.(name{1});

    Trials_Info = load([epoched_data_path, '\sub-', num2str(subject), ...
        '\Trials_Info.mat']);
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});


    %% identify the conditions
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = ...
            condition_indices_identifier(Trials_Info, subject);
    elseif strcmp(pressure_score, 'score')
        % Find all conditions indices in trials 
        % S1, S2, and S3 
        condition_indices = ...
            condition_indices_identifier_ScoreBased(Trials_Info, subject);
    end


    %% Creating the tables for classification
    if strcmp(classes(1), 'P')
        P1_count = 0; 
        P3_count = 0;
        P6_count = 0;
        P1_features = [];
        P3_features = [];
        P6_features = [];
        condition1 = condition_indices.P1;
        condition2 = condition_indices.P3;
        condition3 = condition_indices.P6;
    elseif strcmp(classes(1), 'S')
        S1_count = 0;
        S2_count = 0;
        S3_count = 0;
        S1_features = [];
        S2_features = [];
        S3_features = [];
        condition1 = condition_indices.S1;
        condition2 = condition_indices.S2;
        condition3 = condition_indices.S3;
    end



    for j = 1:length(data)

        if isempty(data{1, j}.Classification_Features.per_trial)
            
            if strcmp(classes(1), 'P')
                if ismember(j, condition1) % condition_indices.P1
                    P1_features = cat(1, P1_features, NaN(1, 8));
                    P1_count = P1_count + 1;
                elseif ismember(j, condition2) % condition_indices.P3
                    P3_features = cat(1, P3_features, NaN(1, 8));
                    P3_count = P3_count + 1;
                elseif ismember(j, condition3) % condition_indices.P6
                    P6_features = cat(1, P6_features, NaN(1, 8));
                    P6_count = P6_count + 1;
                end
            elseif strcmp(classes(1), 'S')
                if ismember(j, condition1) % condition_indices.P1
                    S1_features = cat(1, S1_features, NaN(1, 8));
                    S1_count = S1_count + 1;
                elseif ismember(j, condition2) % condition_indices.P3
                    S2_features = cat(1, S2_features, NaN(1, 8));
                    S2_count = S2_count + 1;
                elseif ismember(j, condition3) % condition_indices.P6
                    S3_features = cat(1, S3_features, NaN(1, 8));
                    S3_count = S3_count + 1;
                end
            end

            continue

        end


        if strcmp(classes(1), 'P')
            if ismember(j, condition1) % condition_indices.P1
                P1_features = cat(1, P1_features, data{1, j}.Classification_Features.per_trial);
                P1_count = P1_count + 1;
            elseif ismember(j, condition2) % condition_indices.P3
                P3_features = cat(1, P3_features, data{1, j}.Classification_Features.per_trial);
                P3_count = P3_count + 1;
            elseif ismember(j, condition3) % condition_indices.P6
                P6_features = cat(1, P6_features, data{1, j}.Classification_Features.per_trial);
                P6_count = P6_count + 1;
            end
        elseif strcmp(classes(1), 'S')
            if ismember(j, condition1) % condition_indices.P1
                S1_features = cat(1, S1_features, data{1, j}.Classification_Features.per_trial);
                S1_count = S1_count + 1;
            elseif ismember(j, condition2) % condition_indices.P3
                S2_features = cat(1, S2_features, data{1, j}.Classification_Features.per_trial);
                S2_count = S2_count + 1;
            elseif ismember(j, condition3) % condition_indices.P6
                S3_features = cat(1, S3_features, data{1, j}.Classification_Features.per_trial);
                S3_count = S3_count + 1;
            end
        end
        
    end


    if strcmp(classes(1), 'P')
        all_features = [P1_features; P3_features; P6_features];
        classLabels = [repmat({'P1'}, P1_count, 1); ...
                       repmat({'P3'}, P3_count, 1); ...
                       repmat({'P6'}, P6_count, 1)];
    elseif strcmp(classes(1), 'S')
        all_features = [S1_features; S2_features; S3_features];
        classLabels = [repmat({'S1'}, S1_count, 1); ...
                       repmat({'S2'}, S2_count, 1); ...
                       repmat({'S3'}, S3_count, 1)];
    end
    classLabels = categorical(classLabels); % Convert to categorical

    % table columns
    columnNames = {'Vastus_med_R_Flexsion_RMS', 'Vastus_med_R_Extension_RMS', ...
                   'Rectus_femoris_R_Flexsion_RMS', 'Rectus_femoris_R_Extension_RMS', ...
                   'Gastrocnemius_R_Flexsion_RMS', 'Gastrocnemius_R_Extension_RMS', ...
                   'Biceps_femoris_R_Flexsion_RMS', 'Biceps_femoris_R_Extension_RMS', ...
                   'Label'};

    % Create the table for the current subject
    subjectTable = array2table(all_features, 'VariableNames', columnNames(1:end-1));
    subjectTable.Label = classLabels; % Add the categorical labels


end