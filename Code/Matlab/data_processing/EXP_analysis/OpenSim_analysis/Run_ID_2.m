% Add OpenSim to the MATLAB path
import org.opensim.modeling.*

% Define file paths (Update these with your actual files)
modelFile = 'modified_model_v2.osim';  % OpenSim model file
motionFile = 'Sub18Cond24Epoch5IK.mot'; % Kinematics file
setupFile = 'ID_setup.xml'; % Ensure this file is in the working directory

% Load the Inverse Dynamics Tool with the setup file
idTool = InverseDynamicsTool(setupFile);

% Load the OpenSim model
model = Model(modelFile);
model.initSystem();
idTool.setModel(model);

% Set input motion file (Kinematics)
idTool.setCoordinatesFileName(motionFile);

% Ensure output folder exists
outputFolder = fullfile(pwd, 'ID11'); % Creates 'ID11' in your current working directory

% Ensure the output directory exists
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end


% Update OpenSim tool paths
idTool.setResultsDir(outputFolder); % Set new results directory
outputFile = fullfile(outputFolder, 'inverse_dynamics_results2.sto');
idTool.setOutputGenForceFileName(outputFile);


% Set time range for analysis (automatically detected from motion file)
motionData = Storage(motionFile);
initial_time = motionData.getFirstTime();
final_time = motionData.getLastTime();
idTool.setStartTime(initial_time);
idTool.setEndTime(final_time);

% Set low-pass filter cutoff frequency
idTool.setLowpassCutoffFrequency(3);

% Run the Inverse Dynamics tool
try
    idTool.run();

    idTool.print('ID_DebugSetup.xml');

    % Check OpenSim logs for errors
    logFile = 'out.log'; % Temporary log file
    system(['type ' logFile]); % Display OpenSim log output

catch ME
    error('Inverse Dynamics Tool encountered an error: %s', ME.message);
end

% Verify that the output file was created successfully
if isfile(outputFile)
    fprintf('Inverse Dynamics completed successfully. Results saved in: %s\n', outputFile);
else
    error('Inverse Dynamics failed: Output file was not created.');
end

% Read the output motion file
inverse_dynamic_path = outputFile;
try
    ID = read_motionFile([folderName_ID, '\', 'Sub-5_Trial100_Epoch1_.sto']);
catch ME
    error('Error reading the motion file: %s', ME.message);
end
