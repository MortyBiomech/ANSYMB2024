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
metadata = readtable(fullfile(data_path, '0_source_data', 'Subjects.xlsx'));
subject_list = 5:18;
pressure_score = 'pressure';

class1_mean_color = [0 0.4470 0.7410];
class1_shade_color = [96, 168, 214]/255;
class2_mean_color = [0.4660 0.6740 0.1880];
class2_shade_color = [170, 204, 126]/255;
class3_mean_color = [0.6350 0.0780 0.1840];
class3_shade_color = [221, 168, 177]/255;


for i = length(subject_list):-1:1

    %% Load the data and Trials_Info
    folderName = [epoched_data_path, 'sub-', num2str(subject_list(i))];
    data = load(fullfile(folderName, 'Epochs_FlextoFlex_based.mat'));
    name = fieldnames(data);
    data = data.(name{1});

    Trials_Info = load(fullfile(folderName, 'Trials_Info.mat'));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});

    
    %% Identify the conditions (pressure or score based)
    if strcmp(pressure_score, 'pressure')
        % Find all conditions indices in trials 
        % P1, P3, and P6
        condition_indices = condition_indices_identifier(Trials_Info, subject_list(i));
        condition.class1 = condition_indices.P1;
        condition.class2 = condition_indices.P3;
        condition.class3 = condition_indices.P6;
    % elseif strcmp(pressure_score, 'score')
    %     % Find all conditions indices in trials 
    %     % S1, S2, and S3 condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i))
    %     condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject_list(i));
    %     condition.class1 = condition_indices.S1;
    %     condition.class2 = condition_indices.S2;
    %     condition.class3 = condition_indices.S3;
    end

    epoch_lengths = [];
    for j = 1:length(data)
        if isempty(data{1, j}.EXP_stream.Encoder_angle)
            continue
        end
        if ~strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
            continue
        end

        epochs = data{1, j}.EXP_stream.Encoder_angle;
        epoch_lengths = cat(2, epoch_lengths, ...
            cell2mat(cellfun(@(x) length(x), epochs, 'UniformOutput', false)));
    end

    median_epoch_length = floor(median(epoch_lengths));
    upper_lim_length = median_epoch_length + floor(3*std(epoch_lengths));
    lower_lim_length = median_epoch_length - floor(3*std(epoch_lengths));

    
    %% time-normalization (without time warping)
    P1_timeNorm_encoder = [];
    P3_timeNorm_encoder = [];
    P6_timeNorm_encoder = [];
    P1_diff = [];
    P3_diff = [];
    P6_diff = [];
    for j = 1:length(data)
        if isempty(data{1, j}.EXP_stream.Encoder_angle)
            continue
        end
        if ~strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
            continue
        end

        epochs = data{1, j}.EXP_stream.Encoder_angle;
        epoch_lengths = cell2mat(cellfun(@(x) length(x), epochs, 'UniformOutput', false));
    
        constraint = and(epoch_lengths > lower_lim_length, ...
            epoch_lengths < upper_lim_length);

        times = data{1, j}.EXP_stream.Times;
        temp_angle = cellfun(@(x1, x2) interp1(x1, x2, linspace(x1(1), x1(end), 2*median_epoch_length), "linear"), ...
            times(constraint), epochs(constraint), 'UniformOutput', false);
        

        ref_angle = data{1, j}.EXP_stream.Ref_angle;
        temp_difference = cellfun(@(x1, x2) abs(x1 - x2), ...
            epochs, ref_angle, 'UniformOutput', false);
        temp_difference_tNorm = cellfun(@(x1, x2) interp1(x1, x2, linspace(x1(1), x1(end), 2*median_epoch_length), "linear"), ...
            times(constraint), temp_difference(constraint), 'UniformOutput', false);
   
        % for k = 1:length(temp_difference)
        %     if any(temp_difference{1, k} < 0)
        %         j
        %     end
        % end

        if Trials_Info{1, j}.General.Pressure == 1
            P1_timeNorm_encoder = cat(1, P1_timeNorm_encoder, ...
                cell2mat(temp_angle'));
            P1_diff = cat(1, P1_diff, ...
                cell2mat(temp_difference_tNorm'));
        elseif Trials_Info{1, j}.General.Pressure == 3
            P3_timeNorm_encoder = cat(1, P3_timeNorm_encoder, ...
                cell2mat(temp_angle'));
            P3_diff = cat(1, P3_diff, ...
                cell2mat(temp_difference_tNorm'));
        elseif Trials_Info{1, j}.General.Pressure == 6
            P6_timeNorm_encoder = cat(1, P6_timeNorm_encoder, ...
                cell2mat(temp_angle'));
            P6_diff = cat(1, P6_diff, ...
                cell2mat(temp_difference_tNorm'));
        end

    end


    %% plot the encoder angle data
    X = linspace(0, 100, size(P6_timeNorm_encoder, 2));
    % ref angle:
    upperLim = metadata.Upper_Lim(metadata.ID == subject_list(i));
    lowerLim = metadata.Lower_Lim(metadata.ID == subject_list(i));
    Freq = metadata.Average_Freq(metadata.ID == subject_list(i));
   
    t = linspace(0, 1/Freq, size(P6_timeNorm_encoder, 2));
    s = abs((lowerLim - upperLim) / 2)*sin(2*pi*Freq*t + pi/2) + ...
        (((lowerLim - upperLim) / 2) + upperLim);

    fig = figure();
    fig.WindowState = 'maximized'; % Maximizes the figure window
    f = tiledlayout(2, 3);
    ax1 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P1_timeNorm_encoder, 1)+std(P1_timeNorm_encoder, 1, 1), ...
         fliplr(mean(P1_timeNorm_encoder, 1)-std(P1_timeNorm_encoder, 1, 1))], ...
         class1_shade_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P1_timeNorm_encoder, 1), 'LineWidth', 2, 'Color', class1_mean_color)
    plot(X, s, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2)
    title('Knee Angle - P1', 'FontWeight', 'normal')
    ylabel('Degree')
    set(gca, 'FontSize', 14)


    ax2 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P3_timeNorm_encoder, 1)+std(P3_timeNorm_encoder, 0, 1), ...
         fliplr(mean(P3_timeNorm_encoder, 1)-std(P3_timeNorm_encoder, 0, 1))], ...
         class2_shade_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P3_timeNorm_encoder, 1), 'LineWidth', 2, 'Color', class2_mean_color)
    plot(X, s, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2)
    title('Knee Angle - P3', 'FontWeight', 'normal')
    set(gca, 'FontSize', 14)


    ax3 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P6_timeNorm_encoder, 1)+std(P6_timeNorm_encoder, 0, 1), ...
         fliplr(mean(P6_timeNorm_encoder, 1)-std(P6_timeNorm_encoder, 0, 1))], ...
         class3_shade_color, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P6_timeNorm_encoder, 1), 'LineWidth', 2, 'Color', class3_mean_color)
    plot(X, s, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2)
    title('Knee Angle - P6', 'FontWeight', 'normal')
    set(gca, 'FontSize', 14)


    ylims = [min([ax1.YLim(1), ax2.YLim(1), ax3.YLim(1)]), max([ax1.YLim(2), ax2.YLim(2), ax3.YLim(2)])];
    ylim(ax1, ylims)
    ylim(ax2, ylims)
    ylim(ax3, ylims)


    ax4 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P1_diff, 1)+std(P1_diff, 0, 1), ...
         fliplr(mean(P1_diff, 1)-std(P1_diff, 0, 1))], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P1_diff, 1), 'LineWidth', 2, 'Color', 'k')
    title('Tracking Error P1', 'FontWeight', 'normal')
    ylabel('Degree')    
    xlabel('Cycle [%]')
    set(gca, 'FontSize', 14)


    ax5 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P3_diff, 1)+std(P3_diff, 0, 1), ...
         fliplr(mean(P3_diff, 1)-std(P3_diff, 0, 1))], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P3_diff, 1), 'LineWidth', 2, 'Color', 'k')
    title('Tracking Error P3', 'FontWeight', 'normal')
    xlabel('Cycle [%]')
    set(gca, 'FontSize', 14)


    ax6 = nexttile;
    hold on
    fill([X fliplr(X)], ...
        [mean(P6_diff, 1)+std(P6_diff, 0, 1), ...
         fliplr(mean(P6_diff, 1)-std(P6_diff, 0, 1))], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    plot(X, mean(P6_diff, 1), 'LineWidth', 2, 'Color', 'k')
    title('Tracking Error P6', 'FontWeight', 'normal')
    xlabel('Cycle [%]')

    ylims = [min([ax4.YLim(1), ax5.YLim(1), ax6.YLim(1)]), max([ax4.YLim(2), ax5.YLim(2), ax6.YLim(2)])];
    ylim(ax4, ylims)
    ylim(ax5, ylims)
    ylim(ax6, ylims)
    set(gca, 'FontSize', 14)


    title(f, ['Subject ', num2str(subject_list(i))], 'FontSize', 18)
    hold off;
    drawnow



    %% plot the mean (std) of kinematic error across cycles per epoch
    P1_diff_mean_perCycle = mean(sum(P1_diff, 2)/size(P1_diff, 2));
    P1_diff_std_perCycle  = std(sum(P1_diff, 2)/size(P1_diff, 2));

    P3_diff_mean_perCycle = mean(sum(P3_diff, 2)/size(P3_diff, 2));
    P3_diff_std_perCycle  = std(sum(P3_diff, 2)/size(P3_diff, 2));

    P6_diff_mean_perCycle = mean(sum(P6_diff, 2)/size(P6_diff, 2));
    P6_diff_std_perCycle  = std(sum(P6_diff, 2)/size(P6_diff, 2));

    % check statistical significance - Kruskal-Wallis Test (Non-parametric)
    P1 = sum(P1_diff, 2)/size(P1_diff, 2); [h_p1, p_p1] = lillietest(P1);
    P3 = sum(P3_diff, 2)/size(P3_diff, 2); [h_p3, p_p3] = lillietest(P3);
    P6 = sum(P6_diff, 2)/size(P6_diff, 2); [h_p6, p_p6] = lillietest(P6);

    data_anova = [P1(:); P3(:); P6(:)];
    group = [repmat({'P1'}, length(P1), 1);
             repmat({'P3'}, length(P3), 1);
             repmat({'P6'}, length(P6), 1)];
    
    [p, tbl, stats] = kruskalwallis(data_anova, group, 'off');
    c = multcompare(stats, 'Display', 'off');



    figure(); hold on;
    means = [P1_diff_mean_perCycle, P3_diff_mean_perCycle, P6_diff_mean_perCycle];
    stds  = [P1_diff_std_perCycle, P3_diff_std_perCycle, P6_diff_std_perCycle];
    b = bar(means);
    b.EdgeColor = "none";
    b.FaceColor = "flat";
    b.FaceAlpha = 0.9;
    b.CData(1, :) = class1_mean_color; 
    b.CData(2, :) = class2_mean_color; 
    b.CData(3, :) = class3_mean_color; 

    % Calculate bar centers for accurate error bar placement
    x = b.XEndPoints;
    
    % Add error bars
    errorbar(x, means, stds, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 18);
    set(gca, 'XTick', x, 'XTickLabel', {'P1', 'P3', 'P6'})
    ylabel('Degree')
    title(['Subject ', num2str(subject_list(i)), ...
        ' - Tracking Error per Cycle per Epoch'], 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)

    

    % Add statistical significance
    y_max = max(means + stds) + 1;
    y_line = 0;
    for i_cmp = 1:size(c,1)
        % Get the group indices and p-value from multcompare output
        group1 = c(i_cmp, 1);
        group2 = c(i_cmp, 2);
        p_val   = c(i_cmp, 6);
        
        % Choose significance label based on p-value
        if p_val < 0.001
            sig_label = '***';
        elseif p_val < 0.01
            sig_label = '**';
        elseif p_val < 0.05
            sig_label = '*';
        else
            continue; % Skip this iteration if not significant (NO horizontal line or label)
        end
        
        % Determine the x positions for the two groups from the bar plot
        x1 = x(group1);
        x2 = x(group2);
        
        % Determine the y positions: set the line at a value just above the taller error bar
        y1 = means(group1) + stds(group1);
        y2 = means(group2) + stds(group2);
        y_line = max([y1, y2, y_line]) + 1;  % adjust offset as needed
        
        % Draw the bracket (a small horizontal line with vertical ticks at the ends)
        line([x1, x1, x2, x2], [y_line-0.2, y_line, y_line, y_line-0.2], 'Color', 'k', 'LineWidth', 1);
        
        % Place the significance label above the bracket
        text(mean([x1, x2]), y_line + 0.1, sig_label, 'HorizontalAlignment', 'center', 'FontSize', 14);
        
        % Increase y_max to avoid overlap with the next comparison
        y_max = y_max + 1.5; 
    end
    hold off;
    drawnow


end




%% Plot Knee Angle vs. PAM Pressure


