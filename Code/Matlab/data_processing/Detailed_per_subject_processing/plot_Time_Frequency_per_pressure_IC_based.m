function plot_Time_Frequency_per_pressure_IC_based(Subject_ICs, subject, epoched_data_path, ic, main_tile)

    %% Loop into the subject ICs to extract the time-frequency content

    
    % parameters for Morlet wavelet computation
    fs = 500; % sampling frequency (Hz)
    lower_freq = 1; % lower frequency for extracting the TF content
    upper_freq = 50; % upper frequency for extractign the TF content
    VoicesPerOctave = 30; 
    
    
    
    N = fieldnames(Subject_ICs);
    for i = 1:length(N)
    
        %% Predefine the structure 
        region_name = N{i};

        subject_rows = cell2mat(Subject_ICs.(region_name)(:,1)) == subject;
        if sum(subject_rows == 1) == 0
            continue
        end
        region_data = cell(1, 4);
        disp(['subject ', num2str(subject) , ' - ', region_name])
    
        % Subject IDs
        region_data{1, 1} = subject;
        % ICs IDs
        region_data{1, 2} = ic;

        filename = ['sub-', num2str(region_data{1, 1}),'\Epochs_FlextoFlex_based.mat'];
        data = load(fullfile(epoched_data_path, filename));
        name = fieldnames(data);
        main_data = data.(name{1});

        filename = ['sub-', num2str(region_data{1, 1}),'\Trials_Info.mat'];
        data = load(fullfile(epoched_data_path, filename));
        name = fieldnames(data);
        Trials_Info = data.(name{1});
        
        ic_id = ic;

        for j = 1:length(region_data(:,1))
            
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
    

        %% Create the 3d matrix containing all TF data for the entire ICs (all ICs of one subject together)
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



        k = 1;
        neg_colors2 = neg_colors + k*(neg_colors - neg_colors(1, :)).*(neg_colors - 1); 
        x = flipud(pos_colors);
        pos_colors2 = flipud(x + k*(x - x(1, :)).*(x - 1));
        custom_cmap2 = [neg_colors2; pos_colors2];
    
    


        %% Plot the TF data of each IC in the region

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
            % global_min = [];
            % global_max = [];

            t1 = tiledlayout(main_tile, 10, 18, "TileSpacing", "compact", "Padding", "loose");
            t1.Layout.Tile = 5;
            t1.Layout.TileSpan = [10, 18];

            pressures = {'P1', 'P3', 'P6'};
            for p = 1:3
               
                TF_file_name = ['TF_TimeWarped_', pressures{p}];
                TF_TimeWarped = eval(TF_file_name);
                P_ref = mean(TF_TimeWarped, 2);
                TF_TimeWarped_norm = TF_TimeWarped - P_ref;
            
            
                nexttile(t1, p*6-5, [5, 6])
                
                imagesc(cycle_time, flipud(frequency_vector), flipud(mean(TF_TimeWarped_norm, 3)))
                % colormap("jet") 
                % colormap(custom_cmap)
                axis xy
                % axis square
                set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex')

                xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
                for b = [4, 8, 13, 30]
                    yline(b, 'LineStyle', '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
                end
                
                title(pressures{p}, 'FontWeight', 'normal', 'Interpreter', 'latex', 'FontSize', 10)
                if strcmp(pressures{p}, 'P1')
                    ylabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 10)
                end
                % xlabel('Cycle [%]')
    
                % % Calculate the global range
                % global_min = cat(2, global_min, min(min(mean(TF_TimeWarped_norm, 3))));
                % global_max = cat(2, global_max, max(max(mean(TF_TimeWarped_norm, 3))));
                
            end
    
            % Set symmetric limits around zero for consistent colormap
            % global_limit = max(abs([min(global_min), max(global_max)]));
            allAxes = findall(t1, 'Type', 'Axes');

            CLIMIT_ax1 = max(abs([allAxes(1).CLim(1), allAxes(1).CLim(2)]));
            CLIMIT_ax2 = max(abs([allAxes(2).CLim(1), allAxes(2).CLim(2)]));
            CLIMIT_ax3 = max(abs([allAxes(3).CLim(1), allAxes(3).CLim(2)]));
            % Set CLim for each tile
            allAxes(1).CLim = [-CLIMIT_ax1, CLIMIT_ax1];
            allAxes(2).CLim = [-CLIMIT_ax2, CLIMIT_ax2];
            allAxes(3).CLim = [-CLIMIT_ax3, CLIMIT_ax3];

            % global_max = max([CLIMIT_ax1, CLIMIT_ax2, CLIMIT_ax3]);
            % allAxes(1).CLim = [-global_max, global_max];
            % allAxes(2).CLim = [-global_max, global_max];
            % allAxes(3).CLim = [-global_max, global_max];
            
    
            % title(t, ['Subject ', num2str(region_data{j,1}), ...
            %     ', IC ', num2str(region_data{j,2})], 'FontSize', 16, 'FontWeight', 'normal')
            
    
            cb1 = colorbar(allAxes(1));
            cb1.FontSize = 10;
            cb1.Label.String = 'Power Baseline Corrected [$\mu V^2$]';
            cb1.Label.Interpreter = 'latex';
            cb1.TickLabelInterpreter = "latex";
            cb1.Label.FontSize = 10;
            cb1.Label.Rotation = 90;

            cb2 = colorbar(allAxes(2));
            cb2.FontSize = 10;
            cb2.TickLabelInterpreter = "latex";
            cb2.Label.FontSize = 10;
            cb2.Label.Rotation = 90;

            cb3 = colorbar(allAxes(3));
            cb3.FontSize = 10;
            cb3.TickLabelInterpreter = "latex";
            cb3.Label.FontSize = 10;
            cb3.Label.Rotation = 90;




            % t2 = tiledlayout(main_tile, 5, 18, "TileSpacing", "compact", "Padding", "loose");
            % t2.Layout.Tile = 115;
            % t2.Layout.TileSpan = [5, 18];

            ax4 = nexttile(t1, 91, [5, 6]);
            imagesc(cycle_time, flipud(frequency_vector), ...
                flipud(mean(TF_TimeWarped_P3, 3)-mean(TF_TimeWarped_P1, 3)))
            colormap("jet") 
            % colormap(custom_cmap)
            axis xy
            % axis square
            set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex')

            xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
            for b = [4, 8, 13, 30]
                yline(b, 'LineStyle', '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
            end
            title('P3 - P1', 'FontWeight', 'normal', 'Interpreter', 'latex')
            ylabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 10)
            xlabel('Cycle [\%]', 'Interpreter', 'latex', 'FontSize', 10)



            ax5 = nexttile(t1, 97, [5, 6]);
            imagesc(cycle_time, flipud(frequency_vector), ...
                flipud(mean(TF_TimeWarped_P6, 3)-mean(TF_TimeWarped_P1, 3)))
            colormap("jet") 
            % colormap(custom_cmap)
            axis xy
            % axis square
            set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex')

            xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
            for b = [4, 8, 13, 30]
                yline(b, 'LineStyle', '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
            end
            title('P6 - P1', 'FontWeight', 'normal', 'Interpreter', 'latex')
            % ylabel('Frequency [Hz]')
            xlabel('Cycle [\%]', 'Interpreter', 'latex', 'FontSize', 10)
            


            ax6 = nexttile(t1, 103, [5, 6]);
            imagesc(cycle_time, flipud(frequency_vector), ...
                flipud(mean(TF_TimeWarped_P6, 3)-mean(TF_TimeWarped_P3, 3)))
            % colormap("jet") 
            colormap(custom_cmap2)
            axis xy
            % axis square
            set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex')

            xline(cycle_time(flexion_part_median+1), 'LineStyle', '--', 'LineWidth', 2);
            for b = [4, 8, 13, 30]
                yline(b, 'LineStyle', '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
            end
            title('P6 - P3', 'FontWeight', 'normal', 'Interpreter', 'latex')
            % ylabel('Frequency [Hz]')
            xlabel('Cycle [\%]', 'Interpreter', 'latex', 'FontSize', 10)



            % CLIMIT = [min([ax4.CLim(1), ax5.CLim(1), ax6.CLim(1)]), ...
            %     max([ax4.CLim(2), ax5.CLim(2), ax6.CLim(2)])];
            % global_lim = max(abs(CLIMIT));

            CLIMIT_ax4 = max(abs([ax4.CLim(1), ax4.CLim(2)]));
            CLIMIT_ax5 = max(abs([ax5.CLim(1), ax5.CLim(2)]));
            CLIMIT_ax6 = max(abs([ax6.CLim(1), ax6.CLim(2)]));

            ax4.CLim = [-CLIMIT_ax4, CLIMIT_ax4];
            ax5.CLim = [-CLIMIT_ax5, CLIMIT_ax5];
            ax6.CLim = [-CLIMIT_ax6, CLIMIT_ax6];

            cb4 = colorbar(ax4);
            cb4.FontSize = 10;
            cb4.TickLabelInterpreter = "latex";
            cb4.Label.FontSize = 10;
            cb4.Label.Rotation = 90;


            cb5 = colorbar(ax5);
            cb5.FontSize = 10;
            cb5.TickLabelInterpreter = "latex";
            cb5.Label.FontSize = 10;
            cb5.Label.Rotation = 90;


            cb6 = colorbar(ax6);
            cb6.FontSize = 10;
            cb6.Label.String = '$\Delta$ Power [$\mu V^2$]';
            cb6.Label.Interpreter = 'latex';
            cb6.TickLabelInterpreter = "latex";
            cb6.Label.FontSize = 10;
            cb6.Label.Rotation = 90;
    
            % set(gcf, "Position", [100, 50, 1300, 750])
            
            drawnow
            
        end

        



    end


end