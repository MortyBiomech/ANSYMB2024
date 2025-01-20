function save_calibrated_force(epoch_type, epoched_data_path, subject, ...
    source_data_path, data_path)

    %% Load Force sensor data
    data = load(fullfile([epoched_data_path, 'sub-', num2str(subject), filesep, ...
        epoch_type]));
    name = fieldnames(data);
    data = data.(name{1});
        
    
    %% plot an initial figure to see the force data
    figure()
    hold on
    trial = 100;
    cellfun(@(s, e) plot(s, e, 'LineWidth', 2), ...
        data{1, trial}.EXP_stream.Times, data{1, trial}.EXP_stream.Forces)
    
    title('Sample force data')
    xlabel('Time (s)')
    ylabel('Sensor data')
    hold off
    
    
    %% Force sensor calibration per subject
    data0 = load_xdf([source_data_path, 'sub-', ...
        num2str(subject), '\ses-S001', '\eeg\', ...
        'sub-', num2str(11), '_ses-S001_task-Default_run-001_eeg.xdf']);

    for i = 1:length(data0)
        if strcmp(data0{1, i}.info.name, 'Encoder_Pressure_Preference_Force')
            cell_id = i;
        end
    end

    
    %% find the point to zoom on x axis
    plot(data0{1, cell_id}.time_series(5, :))
    [xSelected, ~] = ginput(1);
    X_max_onplot = xSelected; 
    close all

    %% calibration point for 2 kg weight
    figure();
    plot(data0{1, cell_id}.time_series(5, :))
    xlim([0 X_max_onplot])

    [xSelected, ~] = ginput(2);
    
    xmin = min(xSelected);
    xmax = max(xSelected);
    
    x = 1:size(data0{1, 3}.time_series, 2);
    idx = x >= xmin & x <= xmax;
    y_2_kg = data0{1, 3}.time_series(5, idx);
    mean_y_2_kg = mean(y_2_kg);
    
    close all
    

    %% calibration point for 2.5 kg weight
    figure();
    plot(data0{1, cell_id}.time_series(5, :))
    xlim([0 X_max_onplot])

    [xSelected, ~] = ginput(2);
    
    xmin = min(xSelected);
    xmax = max(xSelected);
    
    x = 1:size(data0{1, 3}.time_series, 2);
    idx = x >= xmin & x <= xmax;
    y_2o5_kg = data0{1, 3}.time_series(5, idx);
    mean_y_2o5_kg = mean(y_2o5_kg);
    
    close all
    

    %% calibration point for 4 kg weight
    figure();
    plot(data0{1, cell_id}.time_series(5, :))
    xlim([0 X_max_onplot])

    [xSelected, ~] = ginput(2);
    
    xmin = min(xSelected);
    xmax = max(xSelected);
    
    x = 1:size(data0{1, 3}.time_series, 2);
    idx = x >= xmin & x <= xmax;
    y_4_kg = data0{1, 3}.time_series(5, idx);
    mean_y_4_kg = mean(y_4_kg);
    
    close all
    

    %% fit a linear model to the points
    F_known = [2 2.5 4]*9.81; % N
    Raw_sensor = [mean_y_2_kg mean_y_2o5_kg mean_y_4_kg];
    
    p = polyfit(Raw_sensor, F_known, 1);  % “1” means a 1st-order (linear) fit
    a = p(1);
    b = p(2);
    

    %% Transforming the raw force data into Newton values
    F_cal = cell(size(data));
    for i = 1:length(data)
        if ~contains(data{1, i}.General.Description, 'Reject')
            F_cal{1, i} = cellfun(@(x) a.*x + b, ...
                data{1, i}.EXP_stream.Forces, 'UniformOutput', false);
        end
    end
    
    %% plot an initial figure to see the calibrated force data
    figure()
    hold on
    trial = 100;
    cellfun(@(s, e) plot(s, e, 'LineWidth', 2), ...
        data{1, trial}.EXP_stream.Times, F_cal{1, trial})

    title('Sample force data')
    xlabel('Time [s]')
    ylabel('Force [N]')

    ylim([0 70])


    %% Save calibrated force data for each subject
    mkdir([data_path, '9_'])

end