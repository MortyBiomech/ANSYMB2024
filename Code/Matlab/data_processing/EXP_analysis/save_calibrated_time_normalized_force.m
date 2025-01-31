function save_calibrated_time_normalized_force(trials_info, data, subject)

    %% Create events for each FlextoFlex epoch per trial
    events = cell(1, length(trials_info));
    for i = 1:length(events)
        flextoflex_start = trials_info{1, i}.Events.EXP_stream.flextoflex_start_indx;
        flextoflex_end   = trials_info{1, i}.Events.EXP_stream.flextoflex_end_indx;
        extension_start  = trials_info{1, i}.Events.EXP_stream.extension_start_indx;
        events{1, i} = [flextoflex_start' extension_start' flextoflex_end'];
    end
    
    
    %% Create main structure for storing force data
    Force_data_structure = struct('General', [], 'Events', [], ...
        'F_cal_not_length_normalized', [], 'F_cal_length_normalized', []);
    Force_data = repmat({Force_data_structure}, 1, length(trials_info));
    for i = 1:length(Force_data)
        Force_data{1, i}.General = trials_info{1, i}.General;
        Force_data{1, i}.Events = events{1, i};
        Force_data{1, i}.F_cal_not_length_normalized = data.F_cal{1, i};
    end
    
    
    %% Find the median length per pressure as a reference
    condition_indices = condition_indices_identifier(trials_info, subject);
    
    Lepoch_p1 = [];
    Lepoch_p3 = [];
    Lepoch_p6 = [];
    
    for i = 1:length(Force_data)
        if ismember(i, condition_indices.P1)
            Lepoch_p1 = cat(1, Lepoch_p1, ...
                [repmat(i, size(Force_data{1, i}.Events, 1), 1) , ...
                 transpose(1:size(Force_data{1, i}.Events, 1)) , ...
                 Force_data{1, i}.Events(:, 3) - Force_data{1, i}.Events(:, 1) + 1]);
        elseif ismember(i, condition_indices.P3)
            Lepoch_p3 = cat(1, Lepoch_p3, ...
                [repmat(i, size(Force_data{1, i}.Events, 1), 1) , ...
                 transpose(1:size(Force_data{1, i}.Events, 1)) , ...
                 Force_data{1, i}.Events(:, 3) - Force_data{1, i}.Events(:, 1) + 1]);
        elseif ismember(i, condition_indices.P6)
            Lepoch_p6 = cat(1, Lepoch_p6, ...
                [repmat(i, size(Force_data{1, i}.Events, 1), 1) , ...
                 transpose(1:size(Force_data{1, i}.Events, 1)) , ...
                 Force_data{1, i}.Events(:, 3) - Force_data{1, i}.Events(:, 1) + 1]);
        end
    end
    

    Lepoch_p1_median = median(Lepoch_p1(:, 3));
    [Lepoch_p1_median_indx, ~] = find(Lepoch_p1(:, 3) == Lepoch_p1_median, 1, "first");
    L_flexion_p1_median = ...
        Force_data{1, Lepoch_p1(Lepoch_p1_median_indx, 1)}.Events(Lepoch_p1(Lepoch_p1_median_indx, 2), 2) - ...
        Force_data{1, Lepoch_p1(Lepoch_p1_median_indx, 1)}.Events(Lepoch_p1(Lepoch_p1_median_indx, 2), 1) + 1;
    L_extension_p1_median = ...
        Force_data{1, Lepoch_p1(Lepoch_p1_median_indx, 1)}.Events(Lepoch_p1(Lepoch_p1_median_indx, 2), 3) - ...
        Force_data{1, Lepoch_p1(Lepoch_p1_median_indx, 1)}.Events(Lepoch_p1(Lepoch_p1_median_indx, 2), 2);

    Lepoch_p3_median = median(Lepoch_p3(:, 3));
    [Lepoch_p3_median_indx, ~] = find(Lepoch_p3(:, 3) == Lepoch_p3_median, 1, "first");
    L_flexion_p3_median = ...
        Force_data{1, Lepoch_p3(Lepoch_p3_median_indx, 1)}.Events(Lepoch_p3(Lepoch_p3_median_indx, 2), 2) - ...
        Force_data{1, Lepoch_p3(Lepoch_p3_median_indx, 1)}.Events(Lepoch_p3(Lepoch_p3_median_indx, 2), 1) + 1;
    L_extension_p3_median = ...
        Force_data{1, Lepoch_p3(Lepoch_p3_median_indx, 1)}.Events(Lepoch_p3(Lepoch_p3_median_indx, 2), 3) - ...
        Force_data{1, Lepoch_p3(Lepoch_p3_median_indx, 1)}.Events(Lepoch_p3(Lepoch_p3_median_indx, 2), 2);

    Lepoch_p6_median = median(Lepoch_p6(:, 3));
    [Lepoch_p6_median_indx, ~] = find(Lepoch_p6(:, 3) == Lepoch_p6_median, 1, "first");
    L_flexion_p6_median = ...
        Force_data{1, Lepoch_p6(Lepoch_p6_median_indx, 1)}.Events(Lepoch_p6(Lepoch_p6_median_indx, 2), 2) - ...
        Force_data{1, Lepoch_p6(Lepoch_p6_median_indx, 1)}.Events(Lepoch_p6(Lepoch_p6_median_indx, 2), 1) + 1;
    L_extension_p6_median = ...
        Force_data{1, Lepoch_p6(Lepoch_p6_median_indx, 1)}.Events(Lepoch_p6(Lepoch_p6_median_indx, 2), 3) - ...
        Force_data{1, Lepoch_p6(Lepoch_p6_median_indx, 1)}.Events(Lepoch_p6(Lepoch_p6_median_indx, 2), 2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % special case for subject 12
    % due to some technical issues force data is corrupted in sessions 3
    % and 4 in this subject and will be ignored in the final illustration
    ignore_trials = [];
    if subject == 12
        for i = 1:length(Force_data)
            if ismember(Force_data{1, i}.General.Session, [3, 4])
                ignore_trials = cat(2, ignore_trials, i);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Force_data = add_length_normalized_force(Force_data, condition_indices, ignore_trials, ...
        L_flexion_p1_median, L_extension_p1_median, ...
        L_flexion_p3_median, L_extension_p3_median, ...
        L_flexion_p6_median, L_extension_p6_median);
    
    
    %% Calculate the mean and std per pressure
    F_p1 = [];
    F_p3 = [];
    F_p6 = [];
    
    for i = 1:length(Force_data)
        if ismember(i, condition_indices.P1)
            F_p1 = cat(1, F_p1, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P3)
            F_p3 = cat(1, F_p3, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P6)
            F_p6 = cat(1, F_p6, Force_data{1, i}.F_cal_length_normalized);
        end
    end
    
    time_p1 = linspace(0, 100, Lepoch_p1_median);
    time_p3 = linspace(0, 100, Lepoch_p3_median);
    time_p6 = linspace(0, 100, Lepoch_p6_median);
    
    
    %% plot the calibrated force data (mean +- std)

    mean_force_p1 = mean(F_p1, 1);
    std_force_p1 = std(F_p1, 1);
    mean_force_p3 = mean(F_p3, 1);
    std_force_p3 = std(F_p3, 1);
    mean_force_p6 = mean(F_p6, 1);
    std_force_p6 = std(F_p6, 1);
    
    All_colours = struct('p1_dark_blue', [0, 0.4470, 0.7410], 'p1_light_blue', [0.3010, 0.7450, 0.9330], ...
            'p3_dark_orange', [0.8500, 0.3250, 0.0980], 'p3_light_orange', [0.9290, 0.6940, 0.1250], ...
            'p6_dark_green', [0.4660, 0.6740, 0.1880], 'p6_light_green', [0.5960, 0.8740, 0.5410]);
    
    figure();
    hold on
    fill([time_p1, fliplr(time_p1)], ...
         [mean_force_p1 + std_force_p1, fliplr(mean_force_p1 - std_force_p1)], ...
         All_colours.p1_light_blue, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Shaded area in light blue
    plot(time_p1, mean_force_p1, 'Color', All_colours.p1_dark_blue,'LineWidth', 2)
    
    fill([time_p3, fliplr(time_p3)], ...
         [mean_force_p3 + std_force_p3, fliplr(mean_force_p3 - std_force_p3)], ...
         All_colours.p3_light_orange, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p3, mean_force_p3, 'Color', All_colours.p3_dark_orange, 'LineWidth', 2)
    
    fill([time_p6, fliplr(time_p6)], ...
         [mean_force_p6 + std_force_p6, fliplr(mean_force_p6 - std_force_p6)], ...
         All_colours.p6_light_green, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p6, mean_force_p6, 'Color', All_colours.p6_dark_green, 'LineWidth', 2)
    
    xline(time_p1(L_flexion_p1_median), 'LineWidth', 2, 'LineStyle', '-.' , 'Color', All_colours.p1_dark_blue)
    xline(time_p3(L_flexion_p3_median), 'LineWidth', 2, 'LineStyle', '-.', 'Color', All_colours.p3_dark_orange)
    xline(time_p6(L_flexion_p6_median), 'LineWidth', 2, 'LineStyle', '-.', 'Color', All_colours.p6_dark_green)

    title('Option 1 - same flexion lenght per pressure', 'Interpreter', 'latex')
    xlabel('Cycle [\%]', 'Interpreter', 'latex')
    ylabel('Force [N]', 'Interpreter', 'latex')
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')






    %% Option 2: same flexion and extension length across all pressures
    Lepoch_p1p3p6 = [Lepoch_p1; Lepoch_p3; Lepoch_p6];
    Lepoch_median = median(Lepoch_p1p3p6(:, 3));
    [Lepoch_median_indx, ~] = find(Lepoch_p1p3p6(:, 3) == Lepoch_median, 1, "first");
    L_flexion_median = ...
        Force_data{1, Lepoch_p1p3p6(Lepoch_median_indx, 1)}.Events(Lepoch_p1p3p6(Lepoch_median_indx, 2), 2) - ...
        Force_data{1, Lepoch_p1p3p6(Lepoch_median_indx, 1)}.Events(Lepoch_p1p3p6(Lepoch_median_indx, 2), 1) + 1;
    L_extension_median = ...
        Force_data{1, Lepoch_p1p3p6(Lepoch_median_indx, 1)}.Events(Lepoch_p1p3p6(Lepoch_median_indx, 2), 3) - ...
        Force_data{1, Lepoch_p1p3p6(Lepoch_median_indx, 1)}.Events(Lepoch_p1p3p6(Lepoch_median_indx, 2), 2);
    L_flexion_p1_median = L_flexion_median;
    L_flexion_p3_median = L_flexion_median;
    L_flexion_p6_median = L_flexion_median;
    L_extension_p1_median = L_extension_median;
    L_extension_p3_median = L_extension_median;
    L_extension_p6_median = L_extension_median;


    Force_data_structure = struct('General', [], 'Events', [], ...
        'F_cal_not_length_normalized', [], 'F_cal_length_normalized', []);
    Force_data = repmat({Force_data_structure}, 1, length(trials_info));
    for i = 1:length(Force_data)
        Force_data{1, i}.General = trials_info{1, i}.General;
        Force_data{1, i}.Events = events{1, i};
        Force_data{1, i}.F_cal_not_length_normalized = data.F_cal{1, i};
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % special case for subject 12
    % due to some technical issues force data is corrupted in sessions 3
    % and 4 in this subject and will be ignored in the final illustration
    ignore_trials = [];
    if subject == 12
        for i = 1:length(Force_data)
            if ismember(Force_data{1, i}.General.Session, [3, 4])
                ignore_trials = cat(2, ignore_trials, i);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Force_data = add_length_normalized_force(Force_data, condition_indices, ignore_trials, ...
        L_flexion_p1_median, L_extension_p1_median, ...
        L_flexion_p3_median, L_extension_p3_median, ...
        L_flexion_p6_median, L_extension_p6_median);



    F_p1 = [];
    F_p3 = [];
    F_p6 = [];
    
    for i = 1:length(Force_data)
        if ismember(i, condition_indices.P1)
            F_p1 = cat(1, F_p1, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P3)
            F_p3 = cat(1, F_p3, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P6)
            F_p6 = cat(1, F_p6, Force_data{1, i}.F_cal_length_normalized);
        end
    end

    mean_force_p1 = mean(F_p1, 1);
    std_force_p1 = std(F_p1, 1);
    mean_force_p3 = mean(F_p3, 1);
    std_force_p3 = std(F_p3, 1);
    mean_force_p6 = mean(F_p6, 1);
    std_force_p6 = std(F_p6, 1);

    time_p1 = linspace(0, 100, Lepoch_median);
    time_p3 = linspace(0, 100, Lepoch_median);
    time_p6 = linspace(0, 100, Lepoch_median);
    

 
    All_colours = struct('p1_dark_blue', [0, 0.4470, 0.7410], 'p1_light_blue', [0.3010, 0.7450, 0.9330], ...
            'p3_dark_orange', [0.8500, 0.3250, 0.0980], 'p3_light_orange', [0.9290, 0.6940, 0.1250], ...
            'p6_dark_green', [0.4660, 0.6740, 0.1880], 'p6_light_green', [0.5960, 0.8740, 0.5410]);
    
    figure();
    hold on
    fill([time_p1, fliplr(time_p1)], ...
         [mean_force_p1 + std_force_p1, fliplr(mean_force_p1 - std_force_p1)], ...
         All_colours.p1_light_blue, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Shaded area in light blue
    plot(time_p1, mean_force_p1, 'Color', All_colours.p1_dark_blue,'LineWidth', 2)
    
    fill([time_p3, fliplr(time_p3)], ...
         [mean_force_p3 + std_force_p3, fliplr(mean_force_p3 - std_force_p3)], ...
         All_colours.p3_light_orange, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p3, mean_force_p3, 'Color', All_colours.p3_dark_orange, 'LineWidth', 2)
    
    fill([time_p6, fliplr(time_p6)], ...
         [mean_force_p6 + std_force_p6, fliplr(mean_force_p6 - std_force_p6)], ...
         All_colours.p6_light_green, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p6, mean_force_p6, 'Color', All_colours.p6_dark_green, 'LineWidth', 2)
    
    xline(time_p1(L_flexion_p1_median), 'LineWidth', 2, 'LineStyle', '-.' , 'Color', All_colours.p1_dark_blue)
    xline(time_p3(L_flexion_p3_median), 'LineWidth', 2, 'LineStyle', '-.', 'Color', All_colours.p3_dark_orange)
    xline(time_p6(L_flexion_p6_median), 'LineWidth', 2, 'LineStyle', '-.', 'Color', All_colours.p6_dark_green)

    title('Option 2 - same flexion lenght across all pressures', 'Interpreter', 'latex')
    xlabel('Cycle [\%]', 'Interpreter', 'latex')
    ylabel('Force [N]', 'Interpreter', 'latex')
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')







    %% Option 3: total FlextoFlex cycles are considered for time normalization (no separate flexion and extension normalization)
    Force_data_structure = struct('General', [], 'Events', [], ...
        'F_cal_not_length_normalized', [], 'F_cal_length_normalized', []);
    Force_data = repmat({Force_data_structure}, 1, length(trials_info));
    for i = 1:length(Force_data)
        Force_data{1, i}.General = trials_info{1, i}.General;
        Force_data{1, i}.Events = events{1, i};
        Force_data{1, i}.F_cal_not_length_normalized = data.F_cal{1, i};
    end
    
    

    Lepoch_p1_median = median(Lepoch_p1(:, 3));
    Lepoch_p3_median = median(Lepoch_p3(:, 3));
    Lepoch_p6_median = median(Lepoch_p6(:, 3));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % special case for subject 12
    % due to some technical issues force data is corrupted in sessions 3
    % and 4 in this subject and will be ignored in the final illustration
    ignore_trials = [];
    if subject == 12
        for i = 1:length(Force_data)
            if ismember(Force_data{1, i}.General.Session, [3, 4])
                ignore_trials = cat(2, ignore_trials, i);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Force_data, p1_extension_start_percent, p3_extension_start_percent, p6_extension_start_percent] = ...
        add_length_normalized_force_option3(Force_data, condition_indices, ignore_trials, ...
        Lepoch_p1_median, Lepoch_p3_median, Lepoch_p6_median);
    
    

    F_p1 = [];
    F_p3 = [];
    F_p6 = [];
    
    for i = 1:length(Force_data)
        if ismember(i, condition_indices.P1)
            F_p1 = cat(1, F_p1, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P3)
            F_p3 = cat(1, F_p3, Force_data{1, i}.F_cal_length_normalized);
        elseif ismember(i, condition_indices.P6)
            F_p6 = cat(1, F_p6, Force_data{1, i}.F_cal_length_normalized);
        end
    end
    
    time_p1 = linspace(0, 100, Lepoch_p1_median);
    time_p3 = linspace(0, 100, Lepoch_p3_median);
    time_p6 = linspace(0, 100, Lepoch_p6_median);
    
    

    
    mean_force_p1 = mean(F_p1, 1);
    std_force_p1 = std(F_p1, 1);
    mean_force_p3 = mean(F_p3, 1);
    std_force_p3 = std(F_p3, 1);
    mean_force_p6 = mean(F_p6, 1);
    std_force_p6 = std(F_p6, 1);
    
    All_colours = struct('p1_dark_blue', [0, 0.4470, 0.7410], 'p1_light_blue', [0.3010, 0.7450, 0.9330], ...
            'p3_dark_orange', [0.8500, 0.3250, 0.0980], 'p3_light_orange', [0.9290, 0.6940, 0.1250], ...
            'p6_dark_green', [0.4660, 0.6740, 0.1880], 'p6_light_green', [0.5960, 0.8740, 0.5410]);
    
    figure();
    hold on
    fill([time_p1, fliplr(time_p1)], ...
         [mean_force_p1 + std_force_p1, fliplr(mean_force_p1 - std_force_p1)], ...
         All_colours.p1_light_blue, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Shaded area in light blue
    plot(time_p1, mean_force_p1, 'Color', All_colours.p1_dark_blue,'LineWidth', 2)
    
    fill([time_p3, fliplr(time_p3)], ...
         [mean_force_p3 + std_force_p3, fliplr(mean_force_p3 - std_force_p3)], ...
         All_colours.p3_light_orange, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p3, mean_force_p3, 'Color', All_colours.p3_dark_orange, 'LineWidth', 2)
    
    fill([time_p6, fliplr(time_p6)], ...
         [mean_force_p6 + std_force_p6, fliplr(mean_force_p6 - std_force_p6)], ...
         All_colours.p6_light_green, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(time_p6, mean_force_p6, 'Color', All_colours.p6_dark_green, 'LineWidth', 2)
    

    plot(mean(p1_extension_start_percent)*100, max(mean_force_p1), ...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    errorbar(mean(p1_extension_start_percent)*100, max(mean_force_p1), ...
        std(p1_extension_start_percent)*100, 'horizontal', 'Color', 'k', 'LineWidth', 2, 'CapSize', 8)
    plot(mean(p3_extension_start_percent)*100, max(mean_force_p3), ...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    errorbar(mean(p3_extension_start_percent)*100, max(mean_force_p3), ...
        std(p3_extension_start_percent)*100, 'horizontal', 'Color', 'k', 'LineWidth', 2, 'CapSize', 8)
    plot(mean(p6_extension_start_percent)*100, max(mean_force_p6), ...
        'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    errorbar(mean(p6_extension_start_percent)*100, max(mean_force_p6), ...
        std(p6_extension_start_percent)*100, 'horizontal', 'Color', 'k', 'LineWidth', 2, 'CapSize', 8)



    title('Option 3 - No separate flexion and extension normalization', 'Interpreter', 'latex')
    xlabel('Cycle [\%]', 'Interpreter', 'latex')
    ylabel('Force [N]', 'Interpreter', 'latex')
    set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex')


end