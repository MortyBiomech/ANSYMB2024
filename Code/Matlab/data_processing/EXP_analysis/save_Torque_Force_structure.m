function save_Torque_Force_structure(epoched_data_path, Exp_analysis_code_path, EXP_Analysis_path)

    %% read the total torque (generated in Opensim analysis folder) 
    subject_list = 5:18;
    
    for i = 1:length(subject_list)
    
        disp(['subject ', num2str(subject_list(i))])
        %% Load Trials_Info
        filepath = [epoched_data_path, 'sub-', num2str(subject_list(i))];
        Trials_Info = load(fullfile(filepath, 'Trials_Info.mat'));
        name = fieldnames(Trials_Info);
        Trials_Info = Trials_Info.(name{1});
    
        %% Load Epoched data
        filepath = [epoched_data_path, 'sub-', num2str(subject_list(i))];
        data = load(fullfile(filepath, 'Epochs_FlextoFlex_based.mat'));
        name = fieldnames(data);
        data = data.(name{1});
        Times = cellfun(@(x) x.EXP_stream.Times, data, 'UniformOutput', false);
    
    
        %% Fill the big structure with Knee Torque and Force sensor data
        folderName_ID = [Exp_analysis_code_path, '\Subjects_data\sub-', ...
                num2str(subject_list(i)), '\ID'];
        % taking the .sto file names
        ID_Files = dir(fullfile(folderName_ID, '*')); % Get all files and folders
        IDNames = {ID_Files(~[ID_Files.isdir]).name}'; % Extract file names and store in a cell array
    
        ID = cell(length(IDNames), 1);
        ID(:, 1) =  cellfun(@(x) read_motionFile([folderName_ID, '\', x]), IDNames, 'UniformOutput', false);
    
        numEntries = length(IDNames); % Get number of entries
        extractedNumbers = zeros(numEntries, 2); % Initialize matrix
        
        for j = 1:numEntries
            str = IDNames{j}; % Get the string
            tokens = regexp(str, 'Trial(\d+)_Epoch(\d+)', 'tokens'); % Extract numbers
            if ~isempty(tokens)
                extractedNumbers(j, :) = str2double(tokens{1}); % Convert to numbers
            end
        end
    
         
        Knee_Torque_Force_Sensor = struct('Torque_raw', [], 'Torque_TimeWarped', [], 'Pressure', [], ...
            'Score', [], 'Description', [], 'Force_sensor_raw', [], 'Force_sensor_TimeWarped', [], ...
            'Flexion_Length', [], 'Extension_Length', [], ...
            'Knee_Angle_raw', [], 'Knee_Angle_TimeWarped', []);
        Knee_Torque_Force_Sensor = repmat({Knee_Torque_Force_Sensor}, size(Trials_Info));
        trials = extractedNumbers(:, 1);
        epochs = extractedNumbers(:, 2);
        for j = 1:length(Trials_Info)
            rows = trials == j;
            [~, sorted_rows] = sort(epochs(rows), 'ascend');
            trial_IDs = ID(rows, 1)'; 
            trial_IDs = trial_IDs(sorted_rows);
            Knee_Torque_Force_Sensor{1, j}.Torque_raw = cellfun(@(x) x.data(:, [1,17]), trial_IDs, 'UniformOutput', false);
            
            Knee_Torque_Force_Sensor{1, j}.Pressure = Trials_Info{1, j}.General.Pressure;
            Knee_Torque_Force_Sensor{1, j}.Score    = Trials_Info{1, j}.General.Score;
    
        end
        Knee_Torque_Force_Sensor = cellfun(@(x, y) ...
            setfield(x, 'Knee_Angle_raw', y.EXP_stream.Encoder_angle), ...
            Knee_Torque_Force_Sensor, data, 'UniformOutput', false);
        
    
        if subject_list(i) < 10
            Knee_Torque_Force_Sensor = cellfun(@(x) setfield(x, 'Description', 'Experiment'), ...
                Knee_Torque_Force_Sensor, 'UniformOutput', false); 
        else
            descriptions = cellfun(@(x) x.General.Description, Trials_Info, ...
                'UniformOutput', false);
            Knee_Torque_Force_Sensor = cellfun(@(x, y) setfield(x, 'Description', y), ...
                Knee_Torque_Force_Sensor, descriptions, 'UniformOutput', false); 
        end
    
    
    
        %% Complete the structure
        % Add Calibrated Force sensor data (subjects 11, 12, 15, 16, 17, 18) 
        % Add time-warped knee torque and force sensor data for above-mentioned subjects

        if ismember(subject_list(i), [11, 12, 15, 16, 17, 18])
            filepath = [EXP_Analysis_path, 'sub-', num2str(subject_list(i))];
            calibrated_Force = load(fullfile(filepath, 'calibrated_Force.mat'));
            name = fieldnames(calibrated_Force);
            calibrated_Force = calibrated_Force.(name{1});
            calibrated_Force = calibrated_Force.F_cal;
        
            % assign the calibrated force data to the main structure
            Knee_Torque_Force_Sensor = cellfun(@(x, y) setfield(x, 'Force_sensor_raw', y), ...
                Knee_Torque_Force_Sensor, calibrated_Force, 'UniformOutput', false);
        end

        Knee_Torque_Force_Sensor = ...
                Torque_Force_timeWarping(Trials_Info, Knee_Torque_Force_Sensor, Times, data, subject_list(i));
        
         
    
        %% Save the Torque-Force structure
        foldertoSave = [EXP_Analysis_path, 'sub-', num2str(subject_list(i))];
        if ~isfolder(foldertoSave)
            mkdir(foldertoSave)
        end
        save(fullfile(foldertoSave, 'KneeTorque_ForceSensor_data.mat'), 'Knee_Torque_Force_Sensor', '-v7.3');
    
    
    end


end