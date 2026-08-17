function [] = writeGRFsToMOT(data, fname)
% Purpose:  Write ground reaction forces applied at COP to a 
%           motion file (fname) for input into the SimTrack
%           workflow.
%
% Input:   GRFTz is a structure containing the following data
%          tStart is the starting time of the data set
%          sF is the sampling frequency of the data
%          fname is the name of the file to be written.
%
% Output:   The file 'fname' is written to the current directory.
% ASeth, 09-07


% Generate column labels for forces, COPs, and vertical torques.
% Order:  rGRF(xyz), rCOP(xyz), lGRF(xyz), lCOP(xyz), rT(xyz), lT(xyz)
label = data.labels;
% label{1} = 'ground_force_1_vx';
% label{2} = 'ground_force_1_vy';
% label{3} = 'ground_force_1_vz';
% label{4} = 'ground_force_1_px';
% label{5} = 'ground_force_1_py';
% label{6} = 'ground_force_1_pz';
% label{7} = 'ground_force_2_vx';
% label{8} = 'ground_force_2_vy';
% label{9} = 'ground_force_2_vz';
% label{10} = 'ground_force_2_px';
% label{11} = 'ground_force_2_py';
% label{12} = 'ground_force_2_pz';
% label{13} = 'ground_torque_1_x';
% label{14} = 'ground_torque_1_y';
% label{15} = 'ground_torque_1_z';
% label{16} = 'ground_torque_2_x';
% label{17} = 'ground_torque_2_y';
% label{18} = 'ground_torque_2_z';
forceIndex = length(label);

    
% Initialize 'motion file data matrix' for writing data of interest.
nRows = length(data.data(:,1));
nCols = length(label);   % plus time
motData = zeros(nRows, nCols);

% Write time array to data matrix.
% time = [tStart:1/sF:(tStart + (nRows-1)/sF)]'; 
% time = GRFTz.time; 
% motData(:, 1) = time;

% Write force data to data matrix.
% NOTE:  each field of mCS.forces has xyz components.
forceData = data.data ;
time = data.data(:,1);
% If the coordinate frame does not have FY as vertical
% if isFZ
%     rot90aboutX = [1 0 0;  0 0 1; 0 -1 0];  
%     forceData = rot3DVectors(rot90aboutX, forceData);
% end
         
motData(:, 1:end) = forceData;          

% Open file for writing.
fid = fopen(fname, 'w');
if fid == -1
    error(['unable to open ', fname])
end

% Write header.
fprintf(fid, 'name %s\n', fname);
fprintf(fid, 'datacolumns %d\n', nCols);
fprintf(fid, 'datarows %d\n', nRows);
fprintf(fid, 'inDegrees=no %s\n','');
%fprintf(fid, 'range %d %d\n', time(1), time(nRows));
fprintf(fid, 'endheader\n\n');

% Write column labels.
%fprintf(fid, '%20s\t', 'time');
for i = 1:nCols
	fprintf(fid, '%20s\t', label{i});
end

% Write data.
for i = 1:nRows
    fprintf(fid, '\n'); 
	for j = 1:nCols
        fprintf(fid, '%20.8f\t', motData(i, j));
    end
end

fclose(fid);
return;

