clc
clear

%% Add and Define Necessary Paths
% main folder containing all codes and data:
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); 
% main processing and data path:
data_processing_path = [main_project_folder, '\Code\Matlab\data_processing'];
data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
% Required path for the code
EMG_features_path = [data_processing_path, ...
    '\EMG_processing\structured_EMG_data'];
ROIs_feature_path = [data_path, '8_Classification\ROIs_features'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
grouplevel_postprocessing_path = [data_processing_path, ...
    '\Group_Level_PostProcessing'];
cleaned_EEG_data_path = [data_path, '5_single-subject-EEG-analysis\'];


%% Main Processing
subject_list = 17; %[5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18];
% Pressures_per_scores = cell(length(subject_list), 10);

for i = 1:length(subject_list)
    
    %% load data
    disp(['subject ', num2str(subject_list(i))])

    filename = ['sub-', num2str(subject_list(i)), '\sub-', num2str(subject_list(i)),'_structured_EMG_data.mat'];
    data = load(fullfile(EMG_features_path, filename));
    name = fieldnames(data);
    data = data.(name{1});

    % for j = 1:length(data)
    % 
    %     if ~strcmp(data{1, j}.Description, 'Experiment')
    %         continue
    %     end
    % 
    %     score = data{1, j}.Score;
    %     pressure = data{1, j}.Pressure;
    % 
    %     if score ~= 0
    %         Pressures_per_scores{i, score} = ...
    %             cat(1, Pressures_per_scores{i, score}, pressure);
    %     end
    % end


    % %% Plot the figure to show the number of pressure conditions per score
    % subject_score_cell = Pressures_per_scores(i, :);
    % 
    % pressure_count = zeros(length(subject_score_cell), 3);
    % for j = 1:length(subject_score_cell)
    %     pressure_count(j, 1) = sum(subject_score_cell{j} == 1);  % Count of 1's
    %     pressure_count(j, 2) = sum(subject_score_cell{j} == 3);  % Count of 3's
    %     pressure_count(j, 3) = sum(subject_score_cell{j} == 6);  % Count of 6's
    % end
    % 
    % % Prepare X (scores), Y (conditions), and Size data for scatter plot
    % [X, Y] = meshgrid(1:10, [1, 3, 6]); % 10 classes (x-axis), 3 conditions (y-axis)
    % X = X(:);
    % Y = Y(:);
    % sizes = pressure_count';         % Transpose to align with Y (conditions)
    % sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    % sizes(sizes == 0) = NaN;
    % 
    % % Create the bubble chart
    % figure;
    % scatter(X, Y, sizes, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor','none');
    % 
    % % Customize the plot
    % xlabel('Scores');
    % title(['Subject ', num2str(subject_list(i)),' - Count of Pressure Conditions per Scores']);
    % yticks([1, 3, 6]);
    % xticks(1:10)
    % yticklabels({'P1', 'P3', 'P6'});
    % set(gca, 'FontSize', 14)
    % grid on
    % 
    % % Add text labels for each bubble to show counts
    % for k = 1:length(X)
    %     if sizes(k) > 0
    %         text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
    %             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w', 'FontWeight', 'bold');
    %     end
    % end
    % 
    % % Adjust axis limits for better visibility
    % xlim([-0.5, 10.5]);
    % ylim([-2, 8]);
    % set(gcf, "Position", [300, 200, 800, 400])


    %% Plot the Muscle forces, create the structures based on pressures
    Vastus_med_R = struct('P1', [], 'P3', [], 'P6', [], ...
        'S1', [], 'S2', [], 'S3', []);
    Rectus_femoris_R = struct('P1', [], 'P3', [], 'P6', [], ...
        'S1', [], 'S2', [], 'S3', []);
    Gastrocnemius_R = struct('P1', [], 'P3', [], 'P6', [], ...
        'S1', [], 'S2', [], 'S3', []);
    Biceps_femoris_R = struct('P1', [], 'P3', [], 'P6', [], ...
        'S1', [], 'S2', [], 'S3', []);

    for j = 1:length(data)
        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end

        non_outliers = data{1, j}.Outlier == 0;
        non_empty_signal_timewarped = cell2mat(cellfun(@(x) ~isempty(x),data{1, j}.Signal_TimeWarped, 'UniformOutput', false));

        temp = data{1, j}.Signal_TimeWarped(1, and(non_outliers(1, :), non_empty_signal_timewarped));
        signals1 = cellfun(@(x) x(1, :), temp, 'UniformOutput', false);
        signals2 = cellfun(@(x) x(2, :), temp, 'UniformOutput', false);
        signals3 = cellfun(@(x) x(3, :), temp, 'UniformOutput', false);
        signals4 = cellfun(@(x) x(4, :), temp, 'UniformOutput', false);

        P = data{1, j}.Pressure;
        switch P
            case 1
                Vastus_med_R.P1 = cat(1, Vastus_med_R.P1, vertcat(signals1{:}));
                Rectus_femoris_R.P1 = cat(1, Rectus_femoris_R.P1, vertcat(signals2{:}));
                Gastrocnemius_R.P1 = cat(1, Gastrocnemius_R.P1, vertcat(signals3{:}));
                Biceps_femoris_R.P1 = cat(1, Biceps_femoris_R.P1, vertcat(signals4{:}));
            case 3
                Vastus_med_R.P3 = cat(1, Vastus_med_R.P3, vertcat(signals1{:}));
                Rectus_femoris_R.P3 = cat(1, Rectus_femoris_R.P3, vertcat(signals2{:}));
                Gastrocnemius_R.P3 = cat(1, Gastrocnemius_R.P3, vertcat(signals3{:}));
                Biceps_femoris_R.P3 = cat(1, Biceps_femoris_R.P3, vertcat(signals4{:}));
            case 6
                Vastus_med_R.P6 = cat(1, Vastus_med_R.P6, vertcat(signals1{:}));
                Rectus_femoris_R.P6 = cat(1, Rectus_femoris_R.P6, vertcat(signals2{:}));
                Gastrocnemius_R.P6 = cat(1, Gastrocnemius_R.P6, vertcat(signals3{:}));
                Biceps_femoris_R.P6 = cat(1, Biceps_femoris_R.P6, vertcat(signals4{:}));
        end


    end

    %% remove the epochs with more than 3times of std of maximums
    max_values = max(Vastus_med_R.P1, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Vastus_med_R.P1 = Vastus_med_R.P1(~remove_idx, :);

    max_values = max(Vastus_med_R.P3, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Vastus_med_R.P3 = Vastus_med_R.P3(~remove_idx, :);

    max_values = max(Vastus_med_R.P6, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Vastus_med_R.P6 = Vastus_med_R.P6(~remove_idx, :);


    max_values = max(Rectus_femoris_R.P1, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Rectus_femoris_R.P1 = Rectus_femoris_R.P1(~remove_idx, :);

    max_values = max(Rectus_femoris_R.P3, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Rectus_femoris_R.P3 = Rectus_femoris_R.P3(~remove_idx, :);

    max_values = max(Rectus_femoris_R.P6, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Rectus_femoris_R.P6 = Rectus_femoris_R.P6(~remove_idx, :);


    max_values = max(Gastrocnemius_R.P1, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Gastrocnemius_R.P1 = Gastrocnemius_R.P1(~remove_idx, :);

    max_values = max(Gastrocnemius_R.P3, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Gastrocnemius_R.P3 = Gastrocnemius_R.P3(~remove_idx, :);

    max_values = max(Gastrocnemius_R.P6, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Gastrocnemius_R.P6 = Gastrocnemius_R.P6(~remove_idx, :);


    max_values = max(Biceps_femoris_R.P1, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Biceps_femoris_R.P1 = Biceps_femoris_R.P1(~remove_idx, :);
    
    max_values = max(Biceps_femoris_R.P3, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Biceps_femoris_R.P3 = Biceps_femoris_R.P3(~remove_idx, :);

    max_values = max(Biceps_femoris_R.P6, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Biceps_femoris_R.P6 = Biceps_femoris_R.P6(~remove_idx, :);


    %% Muscles No_PAM conditions
    
    if subject_list(i) >= 10

        Vastus_med_R_ref = [];
        Rectus_femoris_R_ref = [];
        Gastrocnemius_R_ref = [];
        Biceps_femoris_R_ref = [];
    
        for j = 1:length(data)
            if strcmp(data{1, j}.Description, 'No_PAM')
                
                non_outliers = data{1, j}.Outlier == 0;
                non_empty_signal_timewarped = cell2mat(cellfun(@(x) ~isempty(x),data{1, j}.Signal_TimeWarped, 'UniformOutput', false));
        
                temp = data{1, j}.Signal_TimeWarped(1, and(non_outliers(1, :), non_empty_signal_timewarped));
                signals1 = cellfun(@(x) x(1, :), temp, 'UniformOutput', false);
                signals2 = cellfun(@(x) x(2, :), temp, 'UniformOutput', false);
                signals3 = cellfun(@(x) x(3, :), temp, 'UniformOutput', false);
                signals4 = cellfun(@(x) x(4, :), temp, 'UniformOutput', false);
        
                
                Vastus_med_R_ref = cat(1, Vastus_med_R_ref, vertcat(signals1{:}));
                Rectus_femoris_R_ref = cat(1, Rectus_femoris_R_ref, vertcat(signals2{:}));
                Gastrocnemius_R_ref = cat(1, Gastrocnemius_R_ref, vertcat(signals3{:}));
                Biceps_femoris_R_ref = cat(1, Biceps_femoris_R_ref, vertcat(signals4{:}));
                   
            end
    
        end

    end

    %% remove the epochs with more than 3times of std of maximums
    max_values = max(Vastus_med_R_ref, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Vastus_med_R_ref = Vastus_med_R_ref(~remove_idx, :);

    max_values = max(Rectus_femoris_R_ref, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Rectus_femoris_R_ref = Rectus_femoris_R_ref(~remove_idx, :);

    max_values = max(Gastrocnemius_R_ref, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Gastrocnemius_R_ref = Gastrocnemius_R_ref(~remove_idx, :);

    max_values = max(Biceps_femoris_R_ref, [], 2);
    remove_idx = max_values > (mean(max_values) + 2*std(max_values));
    Biceps_femoris_R_ref = Biceps_femoris_R_ref(~remove_idx, :);




    %% Plot the Muscle forces: 2*2 figure to show all muscles and pressures
    X_axis_cycle = linspace(0, 100, sum(data{1, j}.Flexion_Extension_Lengths)+1);
    end_of_flexion_indx = data{1, j}.Flexion_Extension_Lengths(1,1) + 1;

    figure(); hold on;
    plot_muscle_activity_per_Pressure(X_axis_cycle, end_of_flexion_indx, ...
        Vastus_med_R, Rectus_femoris_R, ...
        Gastrocnemius_R, Biceps_femoris_R, subject_list(i))
    legend({'P1', '', 'P3', '', 'P6', ''})

    if subject_list(i) >= 10
        plot_muscle_activity_NoPAM(X_axis_cycle, ...
            Vastus_med_R_ref, Rectus_femoris_R_ref, ...
            Gastrocnemius_R_ref, Biceps_femoris_R_ref)
        legend({'P1', '', 'P3', '', 'P6', '', '', 'No PAM', ''})
    end


    

end








%% Plot the heatmat of EEG features per person
ROIs_per_pressure = 'ROIs_2_FlextoFlex_per_trial.mat';

for i = 1:numel(subject_list)

    [RMS_Freq_Region, ...
        heatmap_diff_1_2, heatmap_diff_1_3, heatmap_diff_2_3, ...
        significance_1_2, significance_1_3, significance_2_3] = ...
        EEG_features_heatmap_tables_per_pressure(ROIs_feature_path, ROIs_per_pressure, subject_list(i));
    
    plot_EEG_features_heatmap_per_pressure(RMS_Freq_Region, ...
        heatmap_diff_1_2, heatmap_diff_1_3, heatmap_diff_2_3, ...
        significance_1_2, significance_1_3, significance_2_3, subject_list(i))

end



%% Plot the Time-Frequency figures for each subject (pressure based)
ROIs_per_pressure = 'ROIs_2_FlextoFlex_per_trial.mat';

for i = 1:numel(subject_list)

    plot_Time_Frequency_per_pressure(ROIs_feature_path, ...
        ROIs_per_pressure, subject_list(i), epoched_data_path)

end





 
