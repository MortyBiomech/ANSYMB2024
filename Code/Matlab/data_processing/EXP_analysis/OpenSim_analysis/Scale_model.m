clc; clear; close all
%%
model_weight = 75.337;
model_height = 170;
subject_weight = 80;
subject_height = 180;

w_scale = subject_weight / model_weight;
h_scale = subject_height / model_height;

%%
% Add OpenSim to the MATLAB path
import org.opensim.modeling.*

% Load the OpenSim model
modelFile = 'unscaled_generic.osim';  % OpenSim model file
model = Model(modelFile);

% Initialize the model (needed before modifications)
model.initSystem();

% Iterate through each body in the model
bodySet = model.getBodySet();
numBodies = bodySet.getSize();

for i = 0:numBodies-1
    body = bodySet.get(i);  % Get the body
    oldMass = body.getMass();  % Get current mass
    newMass = oldMass * w_scale;  % Add 1 kg
    body.setMass(newMass);  % Update mass
    % fprintf('Updated mass of %s: %.2f -> %.2f kg\n', body.getName(), oldMass, newMass);
end

% Save the modified model
newModelFile = 'modified_model.osim';
model.print(newModelFile);
% fprintf('Modified model saved as %s\n', newModelFile);



%%
% Add OpenSim to the MATLAB path
import org.opensim.modeling.*

% Load the OpenSim model
modelFile = 'modified_model.osim';  % Replace with your actual OpenSim model file
model = Model(modelFile);

% Initialize the model (required before modifications)
model.initSystem();

% Iterate through each body in the model
bodySet = model.getBodySet();
numBodies = bodySet.getSize();

for i = 0:numBodies-1
    body = bodySet.get(i);  % Get the body
    % fprintf('Processing body: %s\n', body.getName());
    
    % Get attached geometry
    geomSet = body.getPropertyByName('attached_geometry'); 
    numGeoms = geomSet.size();
    
    for j = 0:numGeoms-1
        geom = Geometry.safeDownCast(geomSet.getValueAsObject(j)); % Cast to Geometry
        if isempty(geom)
            continue;
        end
        Vec3()
        % Get existing scale factor
        scaleFactor = geom.get_scale_factors();
%         old_sf = vec3ToDouble(scaleFactor);
        old_sf = [scaleFactor.get(0),scaleFactor.get(1),scaleFactor.get(2)];
        % Multiply scale factor by 1.1
        newScaleFactor = old_sf .* h_scale;
        
        geom.set_scale_factors(Vec3(newScaleFactor(1),newScaleFactor(2),newScaleFactor(3)));
        
%         fprintf('Updated scale factor of geometry for body %s: [%.3f, %.3f, %.3f]\n', ...
%                 body.getName(), newScaleFactor.get(0), newScaleFactor.get(1), newScaleFactor.get(2));
    end
end

% Save the modified model
newModelFile = 'modified_model_v2.osim';
model.print(newModelFile);
% fprintf('Modified model saved as %s\n', newModelFile);
