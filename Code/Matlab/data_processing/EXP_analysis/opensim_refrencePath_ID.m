function ID = opensim_refrencePath_ID(Exp_analysis_code_path, subject, t, s)

    cd(Exp_analysis_code_path)

    %% Add OpenSim to the MATLAB path and load needed files
    import org.opensim.modeling.*
    
    % Load main IK file
    IK = read_motionFile([Exp_analysis_code_path, '\IK.mot']);
    labels = IK.labels;
    
    % Load the XML setup file for inverse dynamics (ID)
    setupFile = 'ID_setup.xml'; % Ensure this file is in the working directory
    
    % Load the Inverse Dynamics Tool with the setup file
    idTool = InverseDynamicsTool(setupFile);



    %% Load meta-data and create each subject's reference path (IK)

    folderName_IK = [Exp_analysis_code_path, '\Subjects_data\sub-', ...
        num2str(subject), '\reference_path\IK'];
    if ~isfolder(folderName_IK)
        mkdir(folderName_IK)
    end

    cd(folderName_IK)
    generate_mot_file(labels, t, s, subject, 0, 0);
    

    folderName_model = [Exp_analysis_code_path, '\Subjects_data\sub-', ...
        num2str(subject), '\model'];
    cd(folderName_model)

    import org.opensim.modeling.*
    modelFile = 'modified_model_v2.osim';  
    model = Model(modelFile);
    model.initSystem();
    % change the model in ID setup file
    idTool.setModel(model);

    % taking the IK file names
    IK_Files = dir(fullfile(folderName_IK, '*')); % Get all files and folders
    IKNames = {IK_Files(~[IK_Files.isdir]).name}'; % Extract file names and store in a cell array

    folderName_ID = [Exp_analysis_code_path, '\Subjects_data\sub-', ...
        num2str(subject), '\reference_path\ID'];
    if ~isfolder(folderName_ID)
        mkdir(folderName_ID)
    end

    generate_sto_files(IKNames{1}, idTool, folderName_ID, folderName_IK);
    ID_Files = dir(fullfile(folderName_ID, '*')); % Get all files and folders
    IDNames = {ID_Files(~[ID_Files.isdir]).name}'; % Extract file names and store in a cell array
    
    ID = read_motionFile([folderName_ID, '\', IDNames{1}]);
    ID = ID.data(:, [1, 17]);
    



end