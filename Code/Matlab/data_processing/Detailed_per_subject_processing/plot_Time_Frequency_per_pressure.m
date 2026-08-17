function plot_Time_Frequency_per_pressure(ROIs_feature_path, ROIs_per_pressure, subject, epoched_data_path)

    %% Load ROI with subjects and ICs in each brain region
    data = load(fullfile(ROIs_feature_path, ROIs_per_pressure));
    name = fieldnames(data);
    ROIs = data.(name{1});
    
    
    %% Loop into the brain regions to extract the time-frequency content
    
    % parameters for Morlet wavelet computation
    fs = 500; % sampling frequency (Hz)
    lower_freq = 1; % lower frequency for extracting the TF content
    upper_freq = 50; % upper frequency for extractign the TF content
    VoicesPerOctave = 30; 
    
    
    
    N = fieldnames(ROIs);
    for i = 1:length(N)
    
        %% Predefine the structure 
        region_name = N{i};

        subject_rows = cell2mat(ROIs.(region_name)(:,1)) == subject;
        if sum(subject_rows == 1) == 0
            continue
        end
        region_data = cell(sum(subject_rows == 1), 4);
        disp(['subject ', num2str(subject) , ' - ', region_name])
    
        % Subject IDs
        region_data(:, 1) = ROIs.(region_name)(subject_rows, 1);
        % ICs IDs
        region_data(:, 2) = ROIs.(region_name)(subject_rows, 2);

        filename = ['sub-', num2str(region_data{1, 1}),'\Epochs_FlextoFlex_based.mat'];
        data = load(fullfile(epoched_data_path, filename));
        name = fieldnames(data);
        main_data = data.(name{1});

        filename = ['sub-', num2str(region_data{1, 1}),'\Trials_Info.mat'];
        data = load(fullfile(epoched_data_path, filename));
        name = fieldnames(data);
        Trials_Info = data.(name{1});
        
        for j = 1:length(region_data(:,1))
            
            ic_id = region_data{j, 2};
    
            if region_data{j, 1} >= 10 % different saving structure for subject 10 and above
                all_ICs = cellfun(@(x) x.EEG_stream.Preprocessed.Sources, main_data, 'UniformOutput', false);
                
                inside_structure = struct('Signal', [], 'Pressure', [], 'Description', []);
                region_data{j, 3} = repmat({inside_structure}, size(all_ICs));
                for k = 1:length(all_ICs)
                    region_data{j, 3}{1, k}.Pressure = Trials_Info{1, k}.General.Pressure;
                    region_data{j, 3}{1, k}.Description = Trials_Info{1, k}.General.Description;
                    if ~isempty(all_ICs{1, k})
                        region_data{j, 3}{1, k}.Signal = cellfun(@(x) ...
                            x(ic_id, :), all_ICs{1, k}, 'UniformOutput', false);
                    end
                end
            else
                all_ICs = cellfun(@(x) ...
                    x.EEG_stream.Preprocessed.Time_Domain.Sources.Not_Length_Normalized, ...
                    main_data, 'UniformOutput', false);
                
                inside_structure = struct('Signal', [], 'Pressure', [], 'Description', []);
                region_data{j, 3} = repmat({inside_structure}, size(all_ICs));
                for k = 1:length(all_ICs)
                    region_data{j, 3}{1, k}.Pressure = Trials_Info{1, k}.General.Pressure;
                    region_data{j, 3}{1, k}.Description = 'Experiment';
                    if ~isempty(all_ICs{1, k})
                        region_data{j, 3}{1, k}.Signal = cellfun(@(x) ...
                            x(ic_id, :), all_ICs{1, k}, 'UniformOutput', false);
                    end
                end
            end
    
    
            % Compute the Morlet Wavelet of each signal and store it in the 4th
            % column of region_data
            inside_structure = struct('TF_content_beforTimeWarp', [], ...
                'Frequencies', [],'Event_indx', [], 'TF_content_afterTimeWarp', []);
            region_data{j, 4} = repmat({inside_structure}, size(all_ICs));
            for k = 1:length(all_ICs)
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    [cwt_coeffs, freqs] = ...
                        cellfun(@(x) cwt(x, 'amor', fs, 'FrequencyLimits', [lower_freq upper_freq], 'VoicesPerOctave', VoicesPerOctave), ...
                        region_data{j, 3}{1, k}.Signal, 'UniformOutput', false); 
                    
                    % [cwt_coeffs, freqs, ~] = ...
                    %     cellfun(@(x) morlet_transform(x, fs, frequencies, width), ...
                    %     region_data{j, 3}{1, k}.Signal, 'UniformOutput', false); 
                     
                    region_data{j, 4}{1, k}.TF_content_beforTimeWarp = ...
                        cellfun(@(x) abs(x).^2, cwt_coeffs, 'UniformOutput', false);
                    region_data{j, 4}{1, k}.Frequencies = ...
                        cellfun(@(x) x, freqs, 'UniformOutput', false);
                end
            end
            
            
            % filling the event indexes to perform time-warping afterward
            for k = 1:length(region_data{j,3})
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    a = Trials_Info{1, k}.Events.EEG_stream.Preprocessed.flextoflex_start_indx;
                    b = Trials_Info{1, k}.Events.EEG_stream.Preprocessed.extension_start_indx; 
                    b = b(1:length(a));
                    c = Trials_Info{1, k}.Events.EEG_stream.Preprocessed.flextoflex_end_indx;
                    region_data{j, 4}{1, k}.Event_indx = [a(:) b(:) c(:)];
                end
            end
        
            
        end

        % memory issue!
        clear data main_data Trials_Info all_ICs
    
        %% Performing the linear time-warping
        % first step: finding the median of flexion and extension segments
        events_vector = [];
        for j = 1:length(region_data(:,1))
            for k = 1:length(region_data{j,3})
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    events_vector = cat(1, events_vector, ...
                        region_data{j, 4}{1, k}.Event_indx);
                end
            end
        end
    
        flexion_length = events_vector(:,2) - events_vector(:,1);
        extension_length = events_vector(:,3) - events_vector(:,2);
        flexion_part_median = median(flexion_length);
        extension_part_median = median(extension_length);
    
        S = 3; % median +- S*(Standard Deviation)
        upperlim_flexion_length = flexion_part_median + S*std(flexion_length);
        lowerlim_flexion_length = flexion_part_median - S*std(flexion_length);
        upperlim_extension_length = extension_part_median + S*std(extension_length);
        lowerlim_extension_length = extension_part_median - S*std(extension_length);
    
        clear events_vector flexion_length extension_length
    
        %%
        for j = 1:length(region_data(:,1))
            for k = 1:length(region_data{j,3})
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    for c = 1:length(region_data{j, 4}{1, k}.TF_content_beforTimeWarp)
                        events = region_data{j, 4}{1, k}.Event_indx(c, :);
                        constraint1 = and(events(2) - events(1) > lowerlim_flexion_length, ...
                            events(2) - events(1) < upperlim_flexion_length);
                        constraint2 = and(events(3) - events(2) > lowerlim_extension_length, ...
                            events(3) - events(2) < upperlim_extension_length);
                        if constraint1 && constraint2
                            TF_matrix = region_data{j, 4}{1, k}.TF_content_beforTimeWarp{1, c};
                            
                            flexion_segment = TF_matrix(:, 1:events(2)-events(1)+1);
                            X = 1:size(flexion_segment,2);
                            Xq = linspace(1, size(flexion_segment,2), flexion_part_median);
                            flexion_segment_TimeWarped = interp1(X', flexion_segment', Xq', "linear");
                            flexion_segment_TimeWarped = flexion_segment_TimeWarped';
    
    
                            extension_segment = TF_matrix(:, size(flexion_segment,2)+1:end);
                            X = 1:size(extension_segment,2);
                            Xq = linspace(1, size(extension_segment,2), extension_part_median);
                            extension_segment_segment_TimeWarped = interp1(X', extension_segment', Xq', "linear");
                            extension_segment_segment_TimeWarped = extension_segment_segment_TimeWarped';
    
                            region_data{j, 4}{1, k}.TF_content_afterTimeWarp{1, c} = ...
                                [flexion_segment_TimeWarped, extension_segment_segment_TimeWarped];
                        end
                    end
                end
            end
        end
    
        clear TF_matrix flexion_segment flexion_segment_TimeWarped 
        clear extension_segment extension_segment_segment_TimeWarped
    
        %%
        
        cycle_time = linspace(0, 100, flexion_part_median + extension_part_median);
    
        frequency_vector_lenghts = [];
        frequency_vector_epochs = [];
        frequency_vector_trials = [];
        frequency_vector_ICs = [];
        for j = 1:length(region_data(:,1))
            for k = 1:length(region_data{j,3})
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    
                    freq_vector_length = cellfun(@(x) size(x,1), region_data{j, 4}{1, k}.TF_content_afterTimeWarp);
                    frequency_vector_lenghts = cat(2, frequency_vector_lenghts, freq_vector_length);
                    frequency_vector_epochs = cat(2, frequency_vector_epochs, 1:length(freq_vector_length));
                    frequency_vector_trials = cat(2, frequency_vector_trials, k*ones(1, length(freq_vector_length)));
                    frequency_vector_ICs = cat(2, frequency_vector_ICs, j*ones(1, length(freq_vector_length)));
    
                end
            end
        end
    
        v = frequency_vector_lenghts;
        v( v == 0) = v( v == 0) + 1000;
        [min_freq_vector_length, min_freq_vector_length_indx] = ...
            min(v);
        j = frequency_vector_ICs(min_freq_vector_length_indx);
        k = frequency_vector_trials(min_freq_vector_length_indx);
        c = frequency_vector_epochs(min_freq_vector_length_indx);
        
        frequency_vector = region_data{j, 4}{1, k}.Frequencies{1, c};
    
        clear frequency_vector_lenghts frequency_vector_epochs
        clear frequency_vector_trials frequency_vector_ICs
    

        %% Create the 3d matrix containing all TF data for the entire Region (all ICs of one subject together)
        P1 = 0;
        P3 = 0;
        P6 = 0;

        for j = 1:length(region_data(:,1))
            
            for k = 1:length(region_data{j,3})
                if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                        strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
    
                    % Identify non-empty cells
                    nonEmptyIdx = ~cellfun(@isempty, region_data{j, 4}{1, k}.TF_content_afterTimeWarp);
                    
                    P = region_data{j, 3}{1, k}.Pressure;
                    switch P
                        case 1
                            P1 = P1 + length(nonEmptyIdx);
                        case 3
                            P3 = P3 + length(nonEmptyIdx);
                        case 6
                            P6 = P6 + length(nonEmptyIdx);
                    end
                end
            end

        end
    
    
        %%
    
        TF_TimeWarped_P1 = zeros(min_freq_vector_length, ...
            flexion_part_median + extension_part_median, P1);
        TF_TimeWarped_P3 = zeros(min_freq_vector_length, ...
            flexion_part_median + extension_part_median, P3);
        TF_TimeWarped_P6 = zeros(min_freq_vector_length, ...
            flexion_part_median + extension_part_median, P6);
        
        % Track the number of elements added
        countP1 = 0;
        countP3 = 0;
        countP6 = 0;
    
        for j = 1:length(region_data(:,1))
            for k = 1:length(region_data{j,3})
                
                if isempty(region_data{j, 3}{1, k}.Signal) || ...
                        ~strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                    continue
                end
    
                % Identify non-empty cells
                nonEmptyIdx = ~cellfun(@isempty, region_data{j, 4}{1, k}.TF_content_afterTimeWarp);
                % Extract first N rows from non-empty matrices
                extractedCells = cellfun(@(x) x(1:min_freq_vector_length, :), ...
                    region_data{j, 4}{1, k}.TF_content_afterTimeWarp(nonEmptyIdx), ...
                    'UniformOutput', false);
                % Convert to a 3D matrix (only for non-empty cells)
                resultMatrix = cat(3, extractedCells{:});
                % Find number of new slices
                numSlices = size(resultMatrix, 3);
    
                P = region_data{j, 3}{1, k}.Pressure;
                switch P
                    case 1
                        TF_TimeWarped_P1(:,:,countP1+(1:numSlices)) = single(resultMatrix);
                        countP1 = countP1 + numSlices;
                    case 3
                        TF_TimeWarped_P3(:,:,countP3+(1:numSlices)) = single(resultMatrix);
                        countP3 = countP3 + numSlices;
                    case 6
                        TF_TimeWarped_P6(:,:,countP6+(1:numSlices)) = single(resultMatrix);
                        countP6 = countP6 + numSlices;
                end
                
            end
        end
    
    
        %%
        % Define extreme colors
        synch_color = [214, 40, 40]/255; % Replace with your desired RGB value for the maximum
        desynch_color = [58, 134, 255]/255;  % Replace with your desired RGB value for the minimum
        
        % Define the number of colors for the colormap
        num_colors = 256;
        
        % Create a gradient for the negative side (red to white)
        neg_colors = [linspace(desynch_color(1), 1, num_colors/2)', ...
                      linspace(desynch_color(2), 1, num_colors/2)', ...
                      linspace(desynch_color(3), 1, num_colors/2)'];
        
        % Create a gradient for the positive side (white to blue)
        pos_colors = [linspace(1, synch_color(1), num_colors/2)', ...
                      linspace(1, synch_color(2), num_colors/2)', ...
                      linspace(1, synch_color(3), num_colors/2)'];
        
        % Combine the gradients to create the full colormap
        custom_cmap = [neg_colors; pos_colors];
    
    
        %%
        global_min = [];
        global_max = [];
        figure()
        t = tiledlayout(1, 3, "TileSpacing", "compact");
        for p = {'P1', 'P3', 'P6'}
            
            TF_file_name = ['TF_TimeWarped_', p{1}];
            TF_TimeWarped = eval(TF_file_name);
            P_ref = mean(TF_TimeWarped, 2);
            TF_TimeWarped_norm = TF_TimeWarped - P_ref;
        
        
            nexttile
            
            imagesc(cycle_time, flipud(frequency_vector), flipud(mean(TF_TimeWarped_norm, 3)))
            colormap(custom_cmap)
            axis xy
            axis square
            xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
            title(p{1}, 'FontWeight', 'normal')
            ylabel('Frequency [Hz]')
            xlabel('Cycle [%]')

            set(gca, 'FontSize', 12)
        
            % Calculate the global range for all heatmaps
            global_min = cat(2, global_min, min(min(mean(TF_TimeWarped_norm, 3))));
            global_max = cat(2, global_max, max(max(mean(TF_TimeWarped_norm, 3))));
            
        end

        % Set symmetric limits around zero for consistent colormap
        global_limit = max(abs([min(global_min), max(global_max)]));
        allAxes = findall(t, 'Type', 'Axes');
        % Set CLim for each tile
        allAxes(1).CLim = [-global_limit, global_limit];
        allAxes(2).CLim = [-global_limit, global_limit];
        allAxes(3).CLim = [-global_limit, global_limit];
        

        region_name_plot = strrep(region_name, '_', ' ');
        if length(region_data(:, 1)) > 1
            title(t, ['Subject ', num2str(region_data{j,1}), ' - ', region_name_plot], 'FontSize', 16, 'FontWeight', 'normal')
        else
            title(t, ['Subject ', num2str(region_data{j,1}), ...
                ', IC ', num2str(region_data{j,2}),' - ', region_name_plot], 'FontSize', 16, 'FontWeight', 'normal')
        end

        cb = colorbar(allAxes(1));
        cb.FontSize = 12;
        cb.Label.String = '\Delta Power';
        cb.Label.FontSize = 12;
        cb.Label.Rotation = 90;

        set(gcf, "Position", [100, 200, 1300, 400])

        drawnow


        %% Plot the TF data of each IC in the region
        if length(region_data(:, 1)) > 1

            for j = 1:length(region_data(:,1))
    
                P1 = 0;
                P3 = 0;
                P6 = 0;
                
                for k = 1:length(region_data{j,3})
                    if ~isempty(region_data{j, 3}{1, k}.Signal) && ...
                            strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
            
                        % Identify non-empty cells
                        nonEmptyIdx = ~cellfun(@isempty, region_data{j, 4}{1, k}.TF_content_afterTimeWarp);
                        
                        P = region_data{j, 3}{1, k}.Pressure;
                        switch P
                            case 1
                                P1 = P1 + length(nonEmptyIdx);
                            case 3
                                P3 = P3 + length(nonEmptyIdx);
                            case 6
                                P6 = P6 + length(nonEmptyIdx);
                        end
                    end
                end
            
            
            
            
                %%
                
                TF_TimeWarped_P1 = zeros(min_freq_vector_length, ...
                    flexion_part_median + extension_part_median, P1);
                TF_TimeWarped_P3 = zeros(min_freq_vector_length, ...
                    flexion_part_median + extension_part_median, P3);
                TF_TimeWarped_P6 = zeros(min_freq_vector_length, ...
                    flexion_part_median + extension_part_median, P6);
                
                % Track the number of elements added
                countP1 = 0;
                countP3 = 0;
                countP6 = 0;
                
                for k = 1:length(region_data{j,3})
                    
                    if isempty(region_data{j, 3}{1, k}.Signal) || ...
                            ~strcmp(region_data{j, 3}{1, k}.Description, 'Experiment')
                        continue
                    end
            
                    % Identify non-empty cells
                    nonEmptyIdx = ~cellfun(@isempty, region_data{j, 4}{1, k}.TF_content_afterTimeWarp);
                    % Extract first N rows from non-empty matrices
                    extractedCells = cellfun(@(x) x(1:min_freq_vector_length, :), ...
                        region_data{j, 4}{1, k}.TF_content_afterTimeWarp(nonEmptyIdx), ...
                        'UniformOutput', false);
                    % Convert to a 3D matrix (only for non-empty cells)
                    resultMatrix = cat(3, extractedCells{:});
                    % Find number of new slices
                    numSlices = size(resultMatrix, 3);
            
                    P = region_data{j, 3}{1, k}.Pressure;
                    switch P
                        case 1
                            TF_TimeWarped_P1(:,:,countP1+(1:numSlices)) = single(resultMatrix);
                            countP1 = countP1 + numSlices;
                        case 3
                            TF_TimeWarped_P3(:,:,countP3+(1:numSlices)) = single(resultMatrix);
                            countP3 = countP3 + numSlices;
                        case 6
                            TF_TimeWarped_P6(:,:,countP6+(1:numSlices)) = single(resultMatrix);
                            countP6 = countP6 + numSlices;
                    end
                    
                end


                %%
                global_min = [];
                global_max = [];
                figure()
                t = tiledlayout(1, 3, "TileSpacing", "compact");
                for p = {'P1', 'P3', 'P6'}
                   
                    TF_file_name = ['TF_TimeWarped_', p{1}];
                    TF_TimeWarped = eval(TF_file_name);
                    P_ref = mean(TF_TimeWarped, 2);
                    TF_TimeWarped_norm = TF_TimeWarped - P_ref;
                
                
                    nexttile
                    
                    imagesc(cycle_time, flipud(frequency_vector), flipud(mean(TF_TimeWarped_norm, 3)))
                    colormap(custom_cmap)
                    axis xy
                    axis square
                    xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
                    title(p{1}, 'FontWeight', 'normal')
                    ylabel('Frequency [Hz]')
                    xlabel('Cycle [%]')
        
                    set(gca, 'FontSize', 12)
                
                    % Calculate the global range for all heatmaps
                    global_min = cat(2, global_min, min(min(mean(TF_TimeWarped_norm, 3))));
                    global_max = cat(2, global_max, max(max(mean(TF_TimeWarped_norm, 3))));
                    
                end
        
                % Set symmetric limits around zero for consistent colormap
                global_limit = max(abs([min(global_min), max(global_max)]));
                allAxes = findall(t, 'Type', 'Axes');
                % Set CLim for each tile
                allAxes(1).CLim = [-global_limit, global_limit];
                allAxes(2).CLim = [-global_limit, global_limit];
                allAxes(3).CLim = [-global_limit, global_limit];
                
        
                region_name_plot = strrep(region_name, '_', ' ');
                title(t, ['Subject ', num2str(region_data{j,1}), ...
                    ', IC ', num2str(region_data{j,2}),' - ', region_name_plot], 'FontSize', 16, 'FontWeight', 'normal')
                
        
                cb = colorbar(allAxes(1));
                cb.FontSize = 12;
                cb.Label.String = '\Delta Power';
                cb.Label.FontSize = 12;
                cb.Label.Rotation = 90;
        
                set(gcf, "Position", [100, 200, 1300, 400])
                
                drawnow
                
            end

        end



    end


end