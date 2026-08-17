function generate_sto_files(IKName, idTool, folderName_ID, folderName_IK)

    import org.opensim.modeling.*
    
    % Set input motion file (Kinematics)
    motionFile = [folderName_IK, '\',  IKName];
    idTool.setCoordinatesFileName(motionFile);
    
    
    % Update OpenSim tool paths
    idTool.setResultsDir(folderName_ID); % Set new results directory
    outputFile = [folderName_ID, '\', IKName(1:end-6), 'ID.sto'];
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
    idTool.run();
    

end