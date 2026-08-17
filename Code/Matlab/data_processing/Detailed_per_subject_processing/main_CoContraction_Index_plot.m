clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];
EMG_processing_path = [main_project_folder, '\Code\Matlab\data_processing\EMG_processing'];


%% Main Loop
subject_list = 5:18;

class1_mean_color = [0 0.4470 0.7410];
class2_mean_color = [0.4660 0.6740 0.1880];
class3_mean_color = [0.6350 0.0780 0.1840];

for i = 1:length(subject_list)

    %% load the structured EMG data
    folderName = ['structured_EMG_data\sub-', num2str(subject_list(i))];
    fileName = ['sub-', num2str(subject_list(i)), ...
        '_structured_EMG_data.mat'];
    data = load(fullfile(EMG_processing_path, folderName, fileName));
    name = fieldnames(data);
    data = data.(name{1});


    %% Normalized to the max EMG for each muscle
    signals = cellfun(@(x) x.Signal_TimeWarped, data, 'UniformOutput', false);
    for j = 1:length(signals)
        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end
        max_signals = cellfun(@(x) max(x, [], 2), signals{1, j}, 'UniformOutput', false);
        data{1, j}.Signal_TimeWarped = cellfun(@(x,y) x./y, ...
            signals{1, j}, max_signals, 'UniformOutput', false);
    end


    %% Sum the extensors and flexors
    extensors = cell(size(data));
    flexors = cell(size(data));
    for j = 1:length(data)
        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end
        epochs = data{1, j}.Signal_TimeWarped;
        emptycells = cell2mat(cellfun(@(x) ~isempty(x), epochs, 'UniformOutput', false));
        extensors{1, j} = cellfun(@(x) x(1, :) + x(2, :), epochs(emptycells), 'UniformOutput', false);
        flexors{1, j} = cellfun(@(x) x(3, :) + x(4, :), epochs(emptycells), 'UniformOutput', false);
    end


    %% Calculate the co-contraction index
    CCI = cell(size(data));
    for j = 1:length(data)
        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end
        F = flexors{1, j};
        E = extensors{1, j};
        CCI{1, j} = cellfun(@(x, y) 2.*(min(x, y)./(x + y)), F, E, 'UniformOutput', false);
    end


    %% Group the CCI data based on pressure conditioins
    CCI_P1 = [];
    CCI_P3 = [];
    CCI_P6 = [];
    for j = 1:length(data)

        if ~strcmp(data{1, j}.Description, 'Experiment')
            continue
        end

        P = data{1, j}.Pressure;
        switch P
            case 1
                CCI_P1 = cat(1, CCI_P1, cell2mat(CCI{1, j}'));
            case 3
                CCI_P3 = cat(1, CCI_P3, cell2mat(CCI{1, j}'));
            case 6
                CCI_P6 = cat(1, CCI_P6, cell2mat(CCI{1, j}'));
        end
    
    end

    max_cci_p1 = max(CCI_P1, [], 2);
    indx1 = max_cci_p1 > 1;
    CCI_P1 = CCI_P1(~indx1, :);

    max_cci_p3 = max(CCI_P3, [], 2);
    indx3 = max_cci_p3 > 1;
    CCI_P3 = CCI_P3(~indx3, :);

    max_cci_p6 = max(CCI_P6, [], 2);
    indx6 = max_cci_p6 > 1;
    CCI_P6 = CCI_P6(~indx6, :);

    
    %% Plot Co-Contraction Index
    X = linspace(0, 100, size(CCI_P1, 2));
    figure();
    t = tiledlayout(1, 3);

    ax1 = nexttile(1); hold on;
    fill([X fliplr(X)], ...
        100*[mean(CCI_P1, 1) + std(CCI_P1, 0, 1), ...
        fliplr(mean(CCI_P1, 1) - std(CCI_P1, 0, 1))], class1_mean_color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(X, 100*mean(CCI_P1, 1), 'Color', class1_mean_color, 'LineWidth', 2);
    xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    xlabel('Cycle [%]')
    ylabel('Co-Contraction Index [%]')
    title('P1', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax2 = nexttile(2); hold on;
    fill([X fliplr(X)], ...
        100*[mean(CCI_P3, 1) + std(CCI_P3, 0, 1), ...
        fliplr(mean(CCI_P3, 1) - std(CCI_P3, 0, 1))], class2_mean_color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(X, 100*mean(CCI_P3, 1), 'Color', class2_mean_color, 'LineWidth', 2);
    xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    xlabel('Cycle [%]')
    title('P3', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)


    ax3 = nexttile(3); hold on;
    fill([X fliplr(X)], ...
        100*[mean(CCI_P6, 1) + std(CCI_P6, 0, 1), ...
        fliplr(mean(CCI_P6, 1) - std(CCI_P6, 0, 1))], class3_mean_color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(X, 100*mean(CCI_P6, 1), 'Color', class3_mean_color, 'LineWidth', 2);
    xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    xlabel('Cycle [%]')
    title('P6', 'FontWeight', 'normal')
    set(gca, 'FontSize', 12)

    title(t, ['Subject ', num2str(subject_list(i)), ...
        ' - Knee Flexor/Extensor Muscles CCI'], ...
        'FontSize', 14)

    % YLIM = [min([ax1.YLim(1), ax2.YLim(1), ax3.YLim(1)]), ...
    %     max([ax1.YLim(2), ax2.YLim(2), ax3.YLim(2)])];
    % ylim(ax1, YLIM)
    % ylim(ax2, YLIM)
    % ylim(ax3, YLIM)
    ylim(ax1, [0, 100])
    ylim(ax2, [0, 100])
    ylim(ax3, [0, 100])

    set(gcf, 'Position', [100, 200, 1300, 420])
    drawnow
    % get(gcf, 'Position')



    % %% Plot Co-Contraction Index
    % X = linspace(0, 100, size(CCI_P1, 2));
    % figure();
    % t = tiledlayout(1, 3);
    % 
    % ax1 = nexttile(1); hold on;
    % 
    % plot(X, 100*CCI_P1', 'Color', [0, 0.4470, 0.7410, 0.04], 'LineWidth', 2);
    % xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    % xlabel('Cycle [%]')
    % ylabel('Co-Contraction Index [%]')
    % title('P1', 'FontWeight', 'normal')
    % set(gca, 'FontSize', 12)
    % 
    % 
    % ax2 = nexttile(2); hold on;
    % plot(X, 100*CCI_P3', 'Color', [0.4660, 0.6740, 0.1880, 0.04], 'LineWidth', 2);
    % xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    % xlabel('Cycle [%]')
    % title('P3', 'FontWeight', 'normal')
    % set(gca, 'FontSize', 12)
    % 
    % 
    % ax3 = nexttile(3); hold on;
    % plot(X, 100*CCI_P6', 'Color', [0.6350, 0.0780, 0.1840, 0.04], 'LineWidth', 2);
    % xline(X(data{1, 1}.Flexion_Extension_Lengths(1)), 'LineStyle', '--', 'LineWidth', 2)
    % xlabel('Cycle [%]')
    % title('P6', 'FontWeight', 'normal')
    % set(gca, 'FontSize', 12)
    % 
    % title(t, ['Subject ', num2str(subject_list(i)), ...
    %     ' - Knee Flexor/Extensor Muscles CCI'], ...
    %     'FontSize', 14)
    % 
    % drawnow

end