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


%% Main loop

subjects = [11, 12, 15, 16, 17, 18];
for i = subjects

    %% load the EXP structured data
    fileName = ['sub-', num2str(i), '\KneeTorque_ForceSensor_data.mat'];
    data = load(fullfile(EXP_Analysis_path, fileName));
    name = fieldnames(data);
    data = data.(name{1});


    %% Separate the data based on pressure condition
    PAM_Force_P1 = [];
    Knee_Angle_P1 = [];
    PAM_Force_P3 = [];
    Knee_Angle_P3 = [];
    PAM_Force_P6 = [];
    Knee_Angle_P6 = [];

    for j = 1:length(data)
        if i == 12 && j > 75
            continue
        end

        min_force = min(data{1, j}.Force_sensor_TimeWarped, [], 2);
        data{1, j}.Force_sensor_TimeWarped = ...
            data{1, j}.Force_sensor_TimeWarped - min_force;

        wanted_indx = true(size(min_force'));

        if i == 16
            max_force = max(data{1, j}.Force_sensor_TimeWarped, [], 2);
            max_force_biggerthan10 = max_force > 10;
            wanted_indx = max_force_biggerthan10;
        end
        
        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end

        

        P = data{1, j}.Pressure;
        switch P
            case 1
                PAM_Force_P1 = cat(1, PAM_Force_P1, data{1, j}.Force_sensor_TimeWarped(wanted_indx, :));
                Knee_Angle_P1 = cat(1, Knee_Angle_P1, data{1, j}.Knee_Angle_TimeWarped(wanted_indx, :));
            case 3
                PAM_Force_P3 = cat(1, PAM_Force_P3, data{1, j}.Force_sensor_TimeWarped(wanted_indx, :));
                Knee_Angle_P3 = cat(1, Knee_Angle_P3, data{1, j}.Knee_Angle_TimeWarped(wanted_indx, :));
            case 6
                PAM_Force_P6 = cat(1, PAM_Force_P6, data{1, j}.Force_sensor_TimeWarped(wanted_indx, :));
                Knee_Angle_P6 = cat(1, Knee_Angle_P6, data{1, j}.Knee_Angle_TimeWarped(wanted_indx, :));
        end
        
    end


    %% Plot the result
    flexion_length = data{1, 1}.Flexion_Length;

    P1_color_flx = [0 0.4470 0.7410];
    P1_color_ext = [96, 168, 214]/255;
    P3_color_flx = [0.4660 0.6740 0.1880];
    P3_color_ext = [170, 204, 126]/255;
    P6_color_flx = [0.6350 0.0780 0.1840];
    P6_color_ext = [221, 168, 177]/255;

    figure();
    t = tiledlayout(1,3);
    ax1 = nexttile(1);
    hold on
    for j = 1:size(PAM_Force_P1, 1)
        plot(Knee_Angle_P1(j, flexion_length+1:end), PAM_Force_P1(j, flexion_length+1:end), 'Color', P1_color_ext, 'LineWidth', 0.5);
    end
    for j = 1:size(PAM_Force_P1, 1)
        plot(Knee_Angle_P1(j, 1:flexion_length), PAM_Force_P1(j, 1:flexion_length), 'Color', P1_color_flx, 'LineWidth', 0.5);
    end
    xlabel('Knee Angle [degree]')    
    ylabel('PAM Force [N]')
    title('P1', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)



    ax2 = nexttile(2);
    hold on
    for j = 1:size(PAM_Force_P3, 1)
        plot(Knee_Angle_P3(j, flexion_length+1:end), PAM_Force_P3(j, flexion_length+1:end), 'Color', P3_color_ext, 'LineWidth', 0.5);
    end
    for j = 1:size(PAM_Force_P3, 1)
        plot(Knee_Angle_P3(j, 1:flexion_length), PAM_Force_P3(j, 1:flexion_length), 'Color', P3_color_flx, 'LineWidth', 0.5);
    end
    xlabel('Knee Angle [degree]')
    title('P3', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    
    ax3 = nexttile(3);
    hold on
    for j = 1:size(PAM_Force_P6, 1)
        plot(Knee_Angle_P6(j, flexion_length+1:end), PAM_Force_P6(j, flexion_length+1:end), 'Color', P6_color_ext, 'LineWidth', 0.5);
    end
    for j = 1:size(PAM_Force_P6, 1)
        plot(Knee_Angle_P6(j, 1:flexion_length), PAM_Force_P6(j, 1:flexion_length), 'Color', P6_color_flx, 'LineWidth', 0.5);
    end
    xlabel('Knee Angle [degree]')
    title('P6', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)



    title(t, ['Subject ', num2str(i), ' - PAM Force vs. Knee Angle'],...
        'FontSize', 14)
    
    YLimit = [min([ax1.YLim(1), ax2.YLim(1), ax3.YLim(1)]), ...
        max([ax1.YLim(2), ax2.YLim(2), ax3.YLim(2)])];
    set(ax1, 'YLim', YLimit);
    set(ax2, 'YLim', YLimit);
    set(ax3, 'YLim', YLimit);

    XLimit = [min([ax1.XLim(1), ax2.XLim(1), ax3.XLim(1)]), ...
        max([ax1.XLim(2), ax2.XLim(2), ax3.XLim(2)])];
    set(ax1, 'XLim', XLimit);
    set(ax2, 'XLim', XLimit);
    set(ax3, 'XLim', XLimit);

    set(gcf, 'Position', [100, 100, 1300, 700])
   



end