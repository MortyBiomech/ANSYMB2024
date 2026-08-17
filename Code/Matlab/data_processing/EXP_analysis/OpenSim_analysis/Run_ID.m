% Add OpenSim to the MATLAB path
import org.opensim.modeling.*

% Define file paths (Update these with your actual files)
modelFile = 'modified_model_v2.osim';  % OpenSim model file
motionFile = 'Sub18Cond14Epoch5IK.mot'; % Kinematics file
% Load the XML setup file
setupFile = 'ID_setup.xml'; % Ensure this file is in the working directory

% Load the Inverse Dynamics Tool with the setup file
idTool = InverseDynamicsTool(setupFile);
% Load the OpenSim model
model = Model(modelFile);
model.initSystem();

idTool.setModel(model);

% Set input motion file (Kinematics)
idTool.setCoordinatesFileName(motionFile);

% Set external loads file (optional, update if using external forces)
% idTool.setExternalLoadsFileName('external_loads.xml');

% Set output file name
if ~isfolder('ID10')
    mkdir('ID10')
else
outputFile = 'ID10/inverse_dynamics_results.sto';
end

% outputFile = fullfile(pwd, 'ID11', 'inverse_dynamics_results.sto');
idTool.setOutputGenForceFileName(outputFile);

% Set time range for analysis (optional: set based on motion data)
state = model.initSystem();
motionData = Storage(motionFile);
initial_time = motionData.getFirstTime();
final_time = motionData.getLastTime();
idTool.setStartTime(initial_time);
idTool.setEndTime(final_time);

% Set low-pass filter cutoff frequency
idTool.setLowpassCutoffFrequency(3);
% idTool.print('ID_setup2.xml')
% Run the Inverse Dynamics tool
idTool.run();

% Display completion message
fprintf('Inverse Dynamics completed. Results saved in: %s\n', outputFile);

inverse_dynamic_path = outputFile;
ID = read_motionFile(inverse_dynamic_path);