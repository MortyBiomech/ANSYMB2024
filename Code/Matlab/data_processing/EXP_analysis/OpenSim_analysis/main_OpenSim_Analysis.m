clc
clear


%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
Exp_analysis_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\EXP_analysis\OpenSim_analysis'];


%% Add OpenSim to the MATLAB path and load needed files
import org.opensim.modeling.*

model_weight = 75.337;
model_height = 170;

% Load main IK file
IK = read_motionFile([Exp_analysis_path, '\IK.mot']);
labels = IK.labels;

% Load the XML setup file for inverse dynamics (ID)
setupFile = 'ID_setup.xml'; % Ensure this file is in the working directory

% Load the Inverse Dynamics Tool with the setup file
idTool = InverseDynamicsTool(setupFile);


%% Load meta-data and each subject's epoched data
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;

for i = 1:length(subject_list)
    
    %% load data
    disp(['sub-', num2str(subject_list(i)), ', Generating IK data ...']);
    filename = ['sub-', num2str(subject_list(i)),'\Epochs_FlextoFlex_based.mat'];
    epoched_data = load(fullfile(epoched_data_path, filename));
    name = fieldnames(epoched_data);
    main_data = epoched_data.(name{1});
    clear epoched_data


    % %% Interpolate the encoder data (have better resolution for Opensim ID)
    % for j = 1:length(main_data)
    %     for k = 1:length(main_data{1, j}.EXP_stream.Times)
    %         tt = main_data{1, j}.EXP_stream.Times{1, k};
    %         signal = main_data{1, j}.EXP_stream.Encoder_angle{1, k};
    %         main_data{1, j}.EXP_stream.Encoder_angle{1, k}  = ...
    %             interp1(tt, signal, linspace(tt(1), tt(end), 2*length(tt)), "linear");
    %         main_data{1, j}.EXP_stream.Times{1, k} = linspace(tt(1), tt(end), 2*length(tt));
    %     end
    % end

    
    %% Generate IK data for each epoch 
    folderName_IK = [Exp_analysis_path, '\Subjects_data\sub-', ...
        num2str(subject_list(i)), '\IK'];
    if ~isfolder(folderName_IK)
        mkdir(folderName_IK)
    end

    cd(folderName_IK)
    count = 0;
    for j = 1:length(main_data)
        Times = main_data{1, j}.EXP_stream.Times;
        if isempty(Times)
            continue
        end
        Knee_angles = main_data{1, j}.EXP_stream.Encoder_angle;
        epoch = num2cell(1:length(Times)); 

        
        % %%% remove the epochs with positive knee angles
        % Knee_angles_negativity = true(1, size(Knee_angles, 2));
        % for k = 1:length(Knee_angles)
        %     if any(Knee_angles{1, k} > 0)
        %         Knee_angles_negativity(k) = false;
        %         count = count + 1;
        %     end
        % end
        % Times = Times(Knee_angles_negativity);
        % Knee_angles = Knee_angles(Knee_angles_negativity);
        % epoch = epoch(Knee_angles_negativity);
        
        
        cellfun(@(x1, x2, x3) generate_mot_file(labels, x1, x2, subject_list(i), j, x3), ...
            Times, Knee_angles, epoch, 'UniformOutput', false);
    end


    %% scale generic Opensim model
    subject_weight = metadata.Weight(metadata.ID == subject_list(i));
    subject_height = metadata.Height(metadata.ID == subject_list(i));
    
    w_scale = subject_weight / model_weight;
    h_scale = subject_height / model_height;


    %% step 1
    %%% scale the weight

    cd(Exp_analysis_path)
    import org.opensim.modeling.*

    % Load the reference OpenSim model
    modelFile = 'unscaled_generic.osim';  % OpenSim model file
    model = Model(modelFile);
    
    % Initialize the model (needed before modifications)
    model.initSystem();
    
    % Iterate through each body in the model
    bodySet = model.getBodySet();
    numBodies = bodySet.getSize();

    for j = 0:numBodies-1
        body = bodySet.get(j);  % Get the body
        oldMass = body.getMass();  % Get current mass
        newMass = oldMass * w_scale;  % Add 1 kg
        body.setMass(newMass);  % Update mass
        % fprintf('Updated mass of %s: %.2f -> %.2f kg\n', body.getName(), oldMass, newMass);
    end
    
    folderName_model = [Exp_analysis_path, '\Subjects_data\sub-', ...
        num2str(subject_list(i)), '\model'];
    if ~isfolder(folderName_model)
        mkdir(folderName_model)
    end
    % Save the modified model
    newModelFile = [folderName_model, '\modified_model.osim'];
    model.print(newModelFile);
    % fprintf('Modified model saved as %s\n', newModelFile);

    cd(folderName_model)


    %% step 2
    %%% scale the Height
    import org.opensim.modeling.*
    modelFile = 'modified_model.osim';  
    model = Model(modelFile);
    
    % Initialize the model (required before modifications)
    model.initSystem();
    
    % Iterate through each body in the model
    bodySet = model.getBodySet();
    numBodies = bodySet.getSize();
    
    for j = 0:numBodies-1
        body = bodySet.get(j);  % Get the body
        % fprintf('Processing body: %s\n', body.getName());
        
        % Get attached geometry
        geomSet = body.getPropertyByName('attached_geometry'); 
        numGeoms = geomSet.size();
        
        for k = 0:numGeoms-1
            geom = Geometry.safeDownCast(geomSet.getValueAsObject(k)); % Cast to Geometry
            if isempty(geom)
                continue;
            end
            Vec3();
            % Get existing scale factor
            scaleFactor = geom.get_scale_factors();
            % old_sf = vec3ToDouble(scaleFactor);
            old_sf = [scaleFactor.get(0),scaleFactor.get(1),scaleFactor.get(2)];
            % Multiply scale factor by 1.1
            newScaleFactor = old_sf .* h_scale;
            
            geom.set_scale_factors(Vec3(newScaleFactor(1),newScaleFactor(2),newScaleFactor(3)));
            
        end
    end

    % Save the modified model
    newModelFile = [folderName_model, '\modified_model_v2.osim'];
    model.print(newModelFile);
    % fprintf('Modified model saved as %s\n', newModelFile);


    

    %% Inverse Dynamic (ID)
    import org.opensim.modeling.*
    modelFile = 'modified_model_v2.osim';  % OpenSim model file
    model = Model(modelFile);
    model.initSystem();
    
    % change the model in ID setup file
    idTool.setModel(model);

    % taking the IK file names
    IK_Files = dir(fullfile(folderName_IK, '*')); % Get all files and folders
    IKNames = {IK_Files(~[IK_Files.isdir]).name}'; % Extract file names and store in a cell array

    folderName_ID = [Exp_analysis_path, '\Subjects_data\sub-', ...
        num2str(subject_list(i)), '\ID'];
    if ~isfolder(folderName_ID)
        mkdir(folderName_ID)
    end

    cellfun(@(x1) generate_sto_files(x1, idTool, folderName_ID, folderName_IK), ...
        IKNames, 'UniformOutput', false);


end


% cond = 14;epoch = 5; sub = 18;
% Time = Epochs_FlextoFlex_based{1,cond}.EXP_stream.Times{1,epoch};
% Knee_ang = Epochs_FlextoFlex_based{1,cond}.EXP_stream.Encoder_angle{1,epoch};

% IK = read_motionFile('IK.mot');


