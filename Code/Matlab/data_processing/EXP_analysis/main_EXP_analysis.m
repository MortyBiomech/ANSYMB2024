clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
source_data_path = [data_path, '0_source_data\'];
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
EXP_Analysis_path = [data_path, '9_EXP_Analysis\'];
Exp_analysis_code_path = [main_project_folder, ...
    '\Code\Matlab\data_processing\EXP_analysis\OpenSim_analysis'];


%% Transform and Save the Calibrated Force (not length normalized)
% epoch_type = 'Epochs_FlextoFlex_based.mat';
% subject = 11; % [11, 12, 15, 16, 17, 18]; These 6 datasets have force sensor data
% 
% calibrated_Force = save_calibrated_force(epoch_type, epoched_data_path, subject, ...
%     source_data_path, data_path);

% Note: This code was executed once and the force sensor data was saved.
% Just load them when you need. 


%% Save the Torque_Force Time-Warped data per subject
% save_Torque_Force_structure(epoched_data_path, Exp_analysis_code_path, EXP_Analysis_path)

% Note: This function was executed once and the structures were saved.
% Just load them when you need them. 
% You can find them in:
% -> EXP_Analysis_path = [data_path, '9_EXP_Analysis\']



%% Plot Knee Torque data per subject
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;
pressure_score = 'pressure';

class1_mean_color = [0 0.4470 0.7410];
class1_shade_color = [96, 168, 214]/255;
class2_mean_color = [0.4660 0.6740 0.1880];
class2_shade_color = [170, 204, 126]/255;
class3_mean_color = [0.6350 0.0780 0.1840];
class3_shade_color = [221, 168, 177]/255;

for i = 1:numel(subject_list)

    subject_weight = metadata.Weight(metadata.ID == subject_list(i));
   
    %% Load Trials_Info and KneeTorque_ForceSensor data
    folderName = [epoched_data_path, 'sub-', num2str(subject_list(i))];
    Trials_Info = load(fullfile(folderName, 'Trials_Info.mat'));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    folderName = [EXP_Analysis_path, 'sub-', num2str(subject_list(i))];
    Torque_Force = load(fullfile(folderName, 'KneeTorque_ForceSensor_data.mat'));
    name = fieldnames(Torque_Force);
    Torque_Force = Torque_Force.(name{1});


    %% Plot the Knee torque data per pressure condition
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject_list(i));
        condition.class1 = condition_indices.P1;
        condition.class2 = condition_indices.P3;
        condition.class3 = condition_indices.P6;
    elseif strcmp(pressure_score, 'score')
        % Find all conditions indices in trials 
        % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
        condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i));
        condition.class1 = condition_indices.S1;
        condition.class2 = condition_indices.S2;
        condition.class3 = condition_indices.S3;
    end

    
    KneeTorque_class1 = []; % P1 or S1
    KneeTorque_class2 = []; % P3 or S2
    KneeTorque_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if ismember(j, condition.class1)
            KneeTorque_class1 = cat(1, KneeTorque_class1, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class2)
            KneeTorque_class2 = cat(1, KneeTorque_class2, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class3)
            KneeTorque_class3 = cat(1, KneeTorque_class3, Torque_Force{1, j}.Torque_TimeWarped);
        end
    end
    KneeTorque_class1 = KneeTorque_class1/subject_weight;
    KneeTorque_class2 = KneeTorque_class2/subject_weight;
    KneeTorque_class3 = KneeTorque_class3/subject_weight;

    
    %% plot the mean and std area
    X = linspace(0, 100, size(KneeTorque_class1, 2));
    figure();
    hold on

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class1)+std(KneeTorque_class1,0,1), ...
         fliplr(mean(KneeTorque_class1)-std(KneeTorque_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class2)+std(KneeTorque_class2,0,1), ...
         fliplr(mean(KneeTorque_class2)-std(KneeTorque_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class3)+std(KneeTorque_class3,0,1), ...
         fliplr(mean(KneeTorque_class3)-std(KneeTorque_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'', '', '', 'P1', 'P3', 'P6',''}, 'Box', 'off')
    elseif strcmp(pressure_score, 'score')
        legend({'', '', '', 'S1', 'S2', 'S3',''}, 'Box', 'off')
    end


    xlabel('Cycle(%)')
    ylabel('Total Knee Torque (Nm/kg)')
    title(['Subject ', num2str(subject_list(i))])
    set(gca, 'FontSize', 12)
    
end



%% Add force sensor data for subjects: 11, 12, 15, 16, 17, 18
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;
pressure_score = 'pressure';

class1_mean_color = [0 0.4470 0.7410];
class1_shade_color = [96, 168, 214]/255;
class2_mean_color = [0.4660 0.6740 0.1880];
class2_shade_color = [170, 204, 126]/255;
class3_mean_color = [0.6350 0.0780 0.1840];
class3_shade_color = [221, 168, 177]/255;

for subject = [11, 12, 15, 16, 17, 18]

    subject_weight = metadata.Weight(metadata.ID == subject);
   
    %% Load Trials_Info and KneeTorque_ForceSensor data
    folderName = [epoched_data_path, 'sub-', num2str(subject)];
    Trials_Info = load(fullfile(folderName, 'Trials_Info.mat'));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    folderName = [EXP_Analysis_path, 'sub-', num2str(subject)];
    Torque_Force = load(fullfile(folderName, 'KneeTorque_ForceSensor_data.mat'));
    name = fieldnames(Torque_Force);
    Torque_Force = Torque_Force.(name{1});


    %% collect the Force Sensor data per condition
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject);
        condition.class1 = condition_indices.P1;
        condition.class2 = condition_indices.P3;
        condition.class3 = condition_indices.P6;
    elseif strcmp(pressure_score, 'score')
        % Find all conditions indices in trials 
        % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
        condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject);
        condition.class1 = condition_indices.S1;
        condition.class2 = condition_indices.S2;
        condition.class3 = condition_indices.S3;
    end

    
    ForceSensor_class1 = []; % P1 or S1
    ForceSensor_class2 = []; % P3 or S2
    ForceSensor_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            ForceSensor_class1 = cat(1, ForceSensor_class1, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class2)
            ForceSensor_class2 = cat(1, ForceSensor_class2, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class3)
            ForceSensor_class3 = cat(1, ForceSensor_class3, Torque_Force{1, j}.Force_sensor_TimeWarped);
        end
    end

    lever_arm = 0.2;
    ForceSensor_class1 = -ForceSensor_class1*lever_arm/subject_weight;
    ForceSensor_class2 = -ForceSensor_class2*lever_arm/subject_weight;
    ForceSensor_class3 = -ForceSensor_class3*lever_arm/subject_weight;
    % ForceSensor_class1 = -ForceSensor_class1;
    % ForceSensor_class2 = -ForceSensor_class2;
    % ForceSensor_class3 = -ForceSensor_class3;

    KneeTorque_class1 = []; % P1 or S1
    KneeTorque_class2 = []; % P3 or S2
    KneeTorque_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            KneeTorque_class1 = cat(1, KneeTorque_class1, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class2)
            KneeTorque_class2 = cat(1, KneeTorque_class2, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class3)
            KneeTorque_class3 = cat(1, KneeTorque_class3, Torque_Force{1, j}.Torque_TimeWarped);
        end
    end
    KneeTorque_class1 = KneeTorque_class1/subject_weight;
    KneeTorque_class2 = KneeTorque_class2/subject_weight;
    KneeTorque_class3 = KneeTorque_class3/subject_weight;


    %% plot the mean and std area
    X = linspace(0, 100, size(ForceSensor_class1, 2));
    fig = figure;
    fig.WindowState = 'maximized'; % Maximizes the figure window
    tiledlayout(2,6)
    ax1 = nexttile(1, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(KneeTorque_class1)+std(KneeTorque_class1,0,1), ...
         fliplr(mean(KneeTorque_class1)-std(KneeTorque_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class2)+std(KneeTorque_class2,0,1), ...
         fliplr(mean(KneeTorque_class2)-std(KneeTorque_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class3)+std(KneeTorque_class3,0,1), ...
         fliplr(mean(KneeTorque_class3)-std(KneeTorque_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on
    
    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total Knee Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)



    ax2 = nexttile(4, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(ForceSensor_class1)+std(ForceSensor_class1,0,1), ...
         fliplr(mean(ForceSensor_class1)-std(ForceSensor_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class2)+std(ForceSensor_class2,0,1), ...
         fliplr(mean(ForceSensor_class2)-std(ForceSensor_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class3)+std(ForceSensor_class3,0,1), ...
         fliplr(mean(ForceSensor_class3)-std(ForceSensor_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(ForceSensor_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(ForceSensor_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(ForceSensor_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'', '', '', 'P1', 'P3', 'P6',''}, 'Box', 'off', 'Location', 'southeast')
    elseif strcmp(pressure_score, 'score')
        legend({'', '', '', 'S1', 'S2', 'S3',''}, 'Box', 'off')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - PAM Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)

    % Find the global y-limits
    ylims = [min([ax1.YLim(1), ax2.YLim(1)]), max([ax1.YLim(2), ax2.YLim(2)])];
    ylim(ax1, ylims)
    ylim(ax2, ylims)


    ax3 = nexttile(7, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class1, 1) - mean(ForceSensor_class1, 1), 'LineWidth', 2, 'Color', class1_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P1 Total Torque', 'P1 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    elseif strcmp(pressure_score, 'score')
        legend({'P1 Total Torque', 'P1 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax4 = nexttile(9, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class2, 1) - mean(ForceSensor_class2, 1), 'LineWidth', 2, 'Color', class2_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P3 Total Torque', 'P3 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    elseif strcmp(pressure_score, 'score')
        legend({'S1 Total Torque', 'S2 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax5 = nexttile(11, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    plot(X, mean(KneeTorque_class3,1) - mean(ForceSensor_class3, 1), 'LineWidth', 2, 'Color', class3_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P6 Total Torque', 'P6 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    elseif strcmp(pressure_score, 'score')
        legend({'S3 Total Torque', 'S3 Biological Torque', ''}, 'Box', 'off', 'Location', 'southeast')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)
    set(gcf, 'Position', [100 50 1100 800])

    % Find the global y-limits
    ylims = [min([ax3.YLim(1), ax4.YLim(1), ax5.YLim(1)]), max([ax3.YLim(2), ax4.YLim(2), ax5.YLim(2)])];
    ylim(ax3, ylims)
    ylim(ax4, ylims)
    ylim(ax5, ylims)

end




%% PAM torques, total torque for pressure conditions + reference path
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;
pressure_score = 'pressure';

class1_mean_color = [0 0.4470 0.7410];
class1_shade_color = [96, 168, 214]/255;
class2_mean_color = [0.4660 0.6740 0.1880];
class2_shade_color = [170, 204, 126]/255;
class3_mean_color = [0.6350 0.0780 0.1840];
class3_shade_color = [221, 168, 177]/255;

for subject = [11, 12, 15, 16, 17, 18]

    subject_weight = metadata.Weight(metadata.ID == subject);
   
    %% Load Trials_Info and KneeTorque_ForceSensor data
    folderName = [epoched_data_path, 'sub-', num2str(subject)];
    Trials_Info = load(fullfile(folderName, 'Trials_Info.mat'));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    folderName = [EXP_Analysis_path, 'sub-', num2str(subject)];
    Torque_Force = load(fullfile(folderName, 'KneeTorque_ForceSensor_data.mat'));
    name = fieldnames(Torque_Force);
    Torque_Force = Torque_Force.(name{1});


    %% collect the Force Sensor data per condition
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject);
        condition.class1 = condition_indices.P1;
        condition.class2 = condition_indices.P3;
        condition.class3 = condition_indices.P6;
    elseif strcmp(pressure_score, 'score')
        % Find all conditions indices in trials 
        % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
        condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject);
        condition.class1 = condition_indices.S1;
        condition.class2 = condition_indices.S2;
        condition.class3 = condition_indices.S3;
    end

    
    ForceSensor_class1 = []; % P1 or S1
    ForceSensor_class2 = []; % P3 or S2
    ForceSensor_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            ForceSensor_class1 = cat(1, ForceSensor_class1, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class2)
            ForceSensor_class2 = cat(1, ForceSensor_class2, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class3)
            ForceSensor_class3 = cat(1, ForceSensor_class3, Torque_Force{1, j}.Force_sensor_TimeWarped);
        end
    end

    lever_arm = 0.2;
    ForceSensor_class1 = -ForceSensor_class1*lever_arm/subject_weight;
    ForceSensor_class2 = -ForceSensor_class2*lever_arm/subject_weight;
    ForceSensor_class3 = -ForceSensor_class3*lever_arm/subject_weight;
    % ForceSensor_class1 = -ForceSensor_class1;
    % ForceSensor_class2 = -ForceSensor_class2;
    % ForceSensor_class3 = -ForceSensor_class3;

    KneeTorque_class1 = []; % P1 or S1
    KneeTorque_class2 = []; % P3 or S2
    KneeTorque_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            KneeTorque_class1 = cat(1, KneeTorque_class1, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class2)
            KneeTorque_class2 = cat(1, KneeTorque_class2, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class3)
            KneeTorque_class3 = cat(1, KneeTorque_class3, Torque_Force{1, j}.Torque_TimeWarped);
        end
    end
    KneeTorque_class1 = KneeTorque_class1/subject_weight;
    KneeTorque_class2 = KneeTorque_class2/subject_weight;
    KneeTorque_class3 = KneeTorque_class3/subject_weight;


    %% reference path torque
    % ref angle:
    upperLim = metadata.Upper_Lim(metadata.ID == subject_list(i));
    lowerLim = metadata.Lower_Lim(metadata.ID == subject_list(i));
    Freq = metadata.Average_Freq(metadata.ID == subject_list(i));
   
    t = linspace(0, 1/Freq, 1000);
    s = abs((lowerLim - upperLim) / 2)*sin(2*pi*Freq*t + pi/2) + ...
        (((lowerLim - upperLim) / 2) + upperLim);

    ID_ref = opensim_refrencePath_ID(Exp_analysis_code_path, subject, t, s);
    [~, extension_start_on_ref] = min(abs(ID_ref(:, 1) - ID(end, 1)/2));
    flexion_ratio = Torque_Force{1,1}.Flexion_Length/...
        (Torque_Force{1,1}.Flexion_Length + Torque_Force{1,1}.Extension_Length);
    extension_ratio = Torque_Force{1,1}.Extension_Length/...
        (Torque_Force{1,1}.Flexion_Length + Torque_Force{1,1}.Extension_Length);
    
    flexion_part_ref = interp1(ID(1:extension_start_on_ref,1), ID(1:extension_start_on_ref,2), ...
        linspace(ID(1,1), ID(extension_start_on_ref,1), floor(flexion_ratio*size(ID_ref, 1))), "linear");
    extension_part_ref = interp1(ID(extension_start_on_ref+1:end,1), ID(extension_start_on_ref+1:end,2), ...
        linspace(ID(extension_start_on_ref+1,1), ID(end,1), ceil(extension_ratio*size(ID_ref, 1))), "linear");
    
    ID_ref_timeWarped = [flexion_part_ref, extension_part_ref]/subject_weight;



    %% plot the mean and std area
    X = linspace(0, 100, size(ForceSensor_class1, 2));
    fig = figure;
    fig.WindowState = 'maximized'; % Maximizes the figure window
    tiledlayout(2,6)
    ax1 = nexttile(1, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(KneeTorque_class1)+std(KneeTorque_class1,0,1), ...
         fliplr(mean(KneeTorque_class1)-std(KneeTorque_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class2)+std(KneeTorque_class2,0,1), ...
         fliplr(mean(KneeTorque_class2)-std(KneeTorque_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class3)+std(KneeTorque_class3,0,1), ...
         fliplr(mean(KneeTorque_class3)-std(KneeTorque_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    X_ref = linspace(0, 100, size(ID_ref_timeWarped, 2));
    plot(X_ref, ID_ref_timeWarped, 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');
   
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'', '', '', 'P1', 'P3', 'P6', 'Reference Path', ''}, 'Box', 'off', 'Location', 'best')
    elseif strcmp(pressure_score, 'score')
        legend({'', '', '', 'S1', 'S2', 'S3',''}, 'Box', 'off')
    end
    
    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total Knee Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)



    ax2 = nexttile(4, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(ForceSensor_class1)+std(ForceSensor_class1,0,1), ...
         fliplr(mean(ForceSensor_class1)-std(ForceSensor_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class2)+std(ForceSensor_class2,0,1), ...
         fliplr(mean(ForceSensor_class2)-std(ForceSensor_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class3)+std(ForceSensor_class3,0,1), ...
         fliplr(mean(ForceSensor_class3)-std(ForceSensor_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(ForceSensor_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(ForceSensor_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(ForceSensor_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    X_ref = linspace(0, 100, size(ID_ref_timeWarped, 2));
    plot(X_ref, ID_ref_timeWarped, 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');
   
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - PAM Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)

    % % Find the global y-limits
    % ylims = [min([ax1.YLim(1), ax2.YLim(1)]), max([ax1.YLim(2), ax2.YLim(2)])];
    % ylim(ax1, ylims)
    % ylim(ax2, ylims)


    ax3 = nexttile(7, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class1, 1) - mean(ForceSensor_class1, 1), 'LineWidth', 2, 'Color', class1_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P1 Total Torque', 'P1 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    elseif strcmp(pressure_score, 'score')
        legend({'P1 Total Torque', 'P1 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax4 = nexttile(9, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class2, 1) - mean(ForceSensor_class2, 1), 'LineWidth', 2, 'Color', class2_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P3 Total Torque', 'P3 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    elseif strcmp(pressure_score, 'score')
        legend({'S1 Total Torque', 'S2 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax5 = nexttile(11, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    plot(X, mean(KneeTorque_class3,1) - mean(ForceSensor_class3, 1), 'LineWidth', 2, 'Color', class3_mean_color, 'LineStyle', '--');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P6 Total Torque', 'P6 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    elseif strcmp(pressure_score, 'score')
        legend({'S3 Total Torque', 'S3 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)
    set(gcf, 'Position', [100 50 1100 800])

    % Find the global y-limits
    ylims = [min([ax3.YLim(1), ax4.YLim(1), ax5.YLim(1)]), max([ax3.YLim(2), ax4.YLim(2), ax5.YLim(2)])];
    ylim(ax3, ylims)
    ylim(ax4, ylims)
    ylim(ax5, ylims)

end







%% PAM torques, total torque for pressure conditions + NO_PAM torque
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;
pressure_score = 'pressure';

class1_mean_color = [0 0.4470 0.7410];
class1_shade_color = [96, 168, 214]/255;
class2_mean_color = [0.4660 0.6740 0.1880];
class2_shade_color = [170, 204, 126]/255;
class3_mean_color = [0.6350 0.0780 0.1840];
class3_shade_color = [221, 168, 177]/255;

for subject = [11, 12, 15, 16, 17, 18]

    subject_weight = metadata.Weight(metadata.ID == subject);
   
    %% Load Trials_Info and KneeTorque_ForceSensor data
    folderName = [epoched_data_path, 'sub-', num2str(subject)];
    Trials_Info = load(fullfile(folderName, 'Trials_Info.mat'));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    folderName = [EXP_Analysis_path, 'sub-', num2str(subject)];
    Torque_Force = load(fullfile(folderName, 'KneeTorque_ForceSensor_data.mat'));
    name = fieldnames(Torque_Force);
    Torque_Force = Torque_Force.(name{1});


    %% collect the Force Sensor data per condition
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject);
        condition.class1 = condition_indices.P1;
        condition.class2 = condition_indices.P3;
        condition.class3 = condition_indices.P6;
    elseif strcmp(pressure_score, 'score')
        % Find all conditions indices in trials 
        % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
        condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject);
        condition.class1 = condition_indices.S1;
        condition.class2 = condition_indices.S2;
        condition.class3 = condition_indices.S3;
    end

    
    ForceSensor_class1 = []; % P1 or S1
    ForceSensor_class2 = []; % P3 or S2
    ForceSensor_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            ForceSensor_class1 = cat(1, ForceSensor_class1, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class2)
            ForceSensor_class2 = cat(1, ForceSensor_class2, Torque_Force{1, j}.Force_sensor_TimeWarped);
        elseif ismember(j, condition.class3)
            ForceSensor_class3 = cat(1, ForceSensor_class3, Torque_Force{1, j}.Force_sensor_TimeWarped);
        end
    end

    lever_arm = 0.2;
    ForceSensor_class1 = -ForceSensor_class1*lever_arm/subject_weight;
    max_force = max(ForceSensor_class1, [], 2);
    ForceSensor_class1 = ForceSensor_class1 - max_force;
    ForceSensor_class2 = -ForceSensor_class2*lever_arm/subject_weight;
    max_force = max(ForceSensor_class2, [], 2);
    ForceSensor_class2 = ForceSensor_class2 - max_force;
    ForceSensor_class3 = -ForceSensor_class3*lever_arm/subject_weight;
    max_force = max(ForceSensor_class3, [], 2);
    ForceSensor_class3 = ForceSensor_class3 - max_force;
    % ForceSensor_class1 = -ForceSensor_class1;
    % ForceSensor_class2 = -ForceSensor_class2;
    % ForceSensor_class3 = -ForceSensor_class3;

    KneeTorque_class1 = []; % P1 or S1
    KneeTorque_class2 = []; % P3 or S2
    KneeTorque_class3 = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'Experiment')
            continue
        end

        if ismember(j, condition.class1)
            KneeTorque_class1 = cat(1, KneeTorque_class1, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class2)
            KneeTorque_class2 = cat(1, KneeTorque_class2, Torque_Force{1, j}.Torque_TimeWarped);
        elseif ismember(j, condition.class3)
            KneeTorque_class3 = cat(1, KneeTorque_class3, Torque_Force{1, j}.Torque_TimeWarped);
        end
    end
    KneeTorque_class1 = KneeTorque_class1/subject_weight;
    KneeTorque_class2 = KneeTorque_class2/subject_weight;
    KneeTorque_class3 = KneeTorque_class3/subject_weight;


    %% No_PAM condition Knee Torque
    KneeTorque_No_PAM = []; % P6 or S3
    for j = 1:length(Torque_Force)
        if subject == 12 && j > 75
            break
        end
        if ~strcmp(Torque_Force{1, j}.Description, 'No_PAM')
            continue
        end
        
        KneeTorque_No_PAM = cat(1, KneeTorque_No_PAM, ...
            Torque_Force{1, j}.Torque_TimeWarped);

    end
    KneeTorque_No_PAM = KneeTorque_No_PAM/subject_weight;


    %% plot the mean and std area
    X = linspace(0, 100, size(ForceSensor_class1, 2));
    fig = figure;
    fig.WindowState = 'maximized'; % Maximizes the figure window
    tiledlayout(2,6)
    ax1 = nexttile(1, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(KneeTorque_class1)+std(KneeTorque_class1,0,1), ...
         fliplr(mean(KneeTorque_class1)-std(KneeTorque_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class2)+std(KneeTorque_class2,0,1), ...
         fliplr(mean(KneeTorque_class2)-std(KneeTorque_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(KneeTorque_class3)+std(KneeTorque_class3,0,1), ...
         fliplr(mean(KneeTorque_class3)-std(KneeTorque_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    plot(X, mean(KneeTorque_No_PAM, 1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');
   
    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'', '', '', 'P1', 'P3', 'P6', 'No-PAM', ''}, 'Box', 'off', 'Location', 'best')
    elseif strcmp(pressure_score, 'score')
        legend({'', '', '', 'S1', 'S2', 'S3',''}, 'Box', 'off')
    end
    
    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total Knee Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)



    ax2 = nexttile(4, [1 3]);
    hold on
    fill([X fliplr(X)], ...
        [mean(ForceSensor_class1)+std(ForceSensor_class1,0,1), ...
         fliplr(mean(ForceSensor_class1)-std(ForceSensor_class1,0,1))], ...
         class1_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class2)+std(ForceSensor_class2,0,1), ...
         fliplr(mean(ForceSensor_class2)-std(ForceSensor_class2,0,1))], ...
         class2_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    fill([X fliplr(X)], ...
        [mean(ForceSensor_class3)+std(ForceSensor_class3,0,1), ...
         fliplr(mean(ForceSensor_class3)-std(ForceSensor_class3,0,1))], ...
         class3_shade_color, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

    plot(X, mean(ForceSensor_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(ForceSensor_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(ForceSensor_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    
    plot(X, mean(KneeTorque_No_PAM, 1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - PAM Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)
    y2lim = get(ax2, 'YLim');
    set(ax2, 'YLim', [y2lim(1), 0])

    % % Find the global y-limits
    % ylims = [min([ax1.YLim(1), ax2.YLim(1)]), max([ax1.YLim(2), ax2.YLim(2)])];
    % ylim(ax1, ylims)
    % ylim(ax2, ylims)


    ax3 = nexttile(7, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class1, 1), 'LineWidth', 3, 'Color', class1_mean_color);
    plot(X, mean(KneeTorque_class1, 1) - mean(ForceSensor_class1, 1), 'LineWidth', 2, 'Color', class1_mean_color, 'LineStyle', '--');
    plot(X, mean(KneeTorque_No_PAM, 1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P1 Total Torque', 'P1 Biological Torque', 'No-PAM Biological', ''}, 'Box', 'off', 'Location', 'northwest')
    elseif strcmp(pressure_score, 'score')
        legend({'P1 Total Torque', 'P1 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax4 = nexttile(9, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class2, 1), 'LineWidth', 3, 'Color', class2_mean_color);
    plot(X, mean(KneeTorque_class2, 1) - mean(ForceSensor_class2, 1), 'LineWidth', 2, 'Color', class2_mean_color, 'LineStyle', '--');
    plot(X, mean(KneeTorque_No_PAM, 1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P3 Total Torque', 'P3 Biological Torque', 'No-PAM Biological', ''}, 'Box', 'off', 'Location', 'northwest')
    elseif strcmp(pressure_score, 'score')
        legend({'S1 Total Torque', 'S2 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax5 = nexttile(11, [1 2]);
    hold on
    plot(X, mean(KneeTorque_class3, 1), 'LineWidth', 3, 'Color', class3_mean_color);
    plot(X, mean(KneeTorque_class3,1) - mean(ForceSensor_class3, 1), 'LineWidth', 2, 'Color', class3_mean_color, 'LineStyle', '--');
    plot(X, mean(KneeTorque_No_PAM, 1), 'LineWidth', 3, 'Color', 'k', 'LineStyle', '-');

    xline(X(Torque_Force{1, 1}.Flexion_Length), 'LineStyle', '--', 'LineWidth', 2)
    grid on

    if strcmp(pressure_score, 'pressure') 
        legend({'P6 Total Torque', 'P6 Biological Torque', 'No-PAM Biological', ''}, 'Box', 'off', 'Location', 'northwest')
    elseif strcmp(pressure_score, 'score')
        legend({'S3 Total Torque', 'S3 Biological Torque', ''}, 'Box', 'off', 'Location', 'best')
    end

    xlabel('Cycle(%)')
    ylabel('Torque (Nm/kg)')
    title(['Subject ', num2str(subject), ' - Total vs. Biological Torque'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)
    set(gcf, 'Position', [100 50 1100 800])

    % Find the global y-limits
    ylims = [min([ax3.YLim(1), ax4.YLim(1), ax5.YLim(1)]), max([ax3.YLim(2), ax4.YLim(2), ax5.YLim(2)])];
    ylim(ax3, ylims)
    ylim(ax4, ylims)
    ylim(ax5, ylims)


end