clc
clear


%% Add paths
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\Code'))
addpath(genpath('D:\Morteza\MyProjects\ANSYMB2024\data\7_STUDY\Epoched_data'))
data_path = 'D:\Morteza\MyProjects\ANSYMB2024\data\';
ersp_data_path = [data_path, '\7_STUDY\Epoched_data\Final_figures', ...
    '\ERSP\Three Pressure Conditions\p 0.01 ersp results\'];
Code_path = ['D:\Morteza\MyProjects\ANSYMB2024\Code\Matlab\data_processing\'];
Subjects_ICs_in_cluster_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\', ...
    'extracting Subjects and ICs in the brain clusters\']
Detailed_Analysis_on_TF_path = [Code_path, ...
    'Group_Level_PostProcessing\Final_paper_plot_generation\', ...
    'Detailed_Analysis_on_TF_regions\'];


%% Load the ERSP Data structures - Right Prim Motor
cluster_name = 'Right Prim Motor';
ersp_structure = load([cluster_name, ' ersp_results.mat']);
ersp_structure = ersp_structure.ersp_results;

ersp_all = ersp_structure.data.allersp;

times = ersp_structure.data.alltimes;
freqs = ersp_structure.data.allfreqs;

pcond = ersp_structure.data.all_pcond{:};

alpha_band = [8, 14];  % Hz
beta_band  = [14, 30]; % Hz

alpha_indx = [find(freqs > alpha_band(1), 1), ...
              find(freqs > alpha_band(2), 1)-1];
beta_indx  = [find(freqs > beta_band(1), 1), ...
              find(freqs > beta_band(2), 1)-1];

%% Main plots
P1_color = [1, 115, 178]/255;
P3_color = [222, 143, 5]/255;
P6_color = [148, 73, 92]/255; %[148, 73, 92]/255;
colors = [P1_color; P3_color; P6_color];


%% Option 1: Mean Power inside the whole cluster regardless of Time & Freq

sig_points = find(pcond == 1);

reshaped_ersp_subj = cellfun(@(x) reshape(x, [], size(ersp_all{1,1}, 3)), ...
    ersp_all, 'UniformOutput', false);
ersp_subj_region = cellfun(@(x) x(sig_points, :), reshaped_ersp_subj, ...
    'UniformOutput', false);


%% -- plot Option 1 ---

% nexttile(1); hold on; % cla; 
% ax11 = subplot(4, 2, 1); hold on; 

fig1 = figure('name', ['Power at significant region - Left Prim Motor'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape');
hold on;
% Jitter setup for spread
spread = 0.15;
rng(1); % Sets the random seed to 1
% x positions (jittered for clarity)
loc = 1:3;
x_jitter = loc + (rand( size(ersp_all{1,1}, 3) ,numel(loc))-0.5)*spread;

subjects_unique = 1:size(ersp_all{1,1}, 3);
n_cond = length(ersp_all);

subj_power_region = [];
for s = 1:numel(subjects_unique)
    ersp_subj_region_P1 = ersp_subj_region{1, 1}(:, s);
    ersp_subj_region_P3 = ersp_subj_region{2, 1}(:, s);
    ersp_subj_region_P6 = ersp_subj_region{3, 1}(:, s);

    subj_power_region = [subj_power_region; ...
        [mean(ersp_subj_region_P1), ...
         mean(ersp_subj_region_P3), ...
         mean(ersp_subj_region_P6)] ];
end




for s = 1:numel(subjects_unique)
    x = x_jitter(s, :);
    subj_power_s = subj_power_region(s, :);
    % 3. Paired line for within-subject
    plot(x, subj_power_s, '-', 'Color', [0.8 0.8 0.8 0.8], 'LineWidth', 0.5);
end


% 1. Box plot layer
for i = 1:n_cond
    
    data = subj_power_region(:, i);

    % boxplot(data, ...
    %     'Whisker', 2.5, ...
    %     'Positions', (i)-spread, ...
    %     'BoxStyle', 'outline', ...
    %     'Colors', colors(i, :), ...
    %     '') % Changes the outlier rule to 2.5×IQR

    % Boxplot: use boxchart (R2019b+), else use boxplot
    boxchart(ones(length(data),1)*(i)-spread, data, ...
        'BoxFaceColor', colors(i,:), ...
        'BoxFaceAlpha', 0.6, ...
        'BoxWidth', 0.18, ...
        'MarkerStyle', 'none', 'LineWidth', 2, ...
        'BoxEdgeColor', colors(i,:)*0.7, ...
        'WhiskerLineColor', colors(i,:)*0.7);
end


% 2. Individual data points (jittered scatter)
for s = 1:numel(subjects_unique)
    subj_power_s = subj_power_region(s, :);
    x = x_jitter(s, :);

    % scatter each subject's points
    for j = 1:numel(loc)
        scatter(x(j), subj_power_s(j), 60, 'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', colors(j,:), ...
            'MarkerFaceAlpha', 0.4, ...
            'MarkerEdgeAlpha', 0.4, ...
            'LineWidth',0.7);
    end
end



% 4. Overlay group means (LMM, or just mean if you prefer)
mu = mean(subj_power_region, 1);
x = [1, 2, 3] - repmat(spread, 1, 3);
plot(x, mu, 'LineWidth', 2, 'Color', [0.3 0.3 0.3])


%% Pairwise comparison tests - Option 1
% Test each condition's normality separately
[h1, p1] = swtest(subj_power_region(:, 1));  % or use lillietest
[h2, p2] = swtest(subj_power_region(:, 2));
[h3, p3] = swtest(subj_power_region(:, 3));

% based on the results of normality test we should use non-parametric tests
% Kruskal-Wallis with Post-hoc Tests
group = {'Low Pressure', 'Medium Pressure', 'High Pressure'};
[p, tbl, stats] = kruskalwallis(subj_power_region, group, 'off');

% if p < 0.05 then:
% Post-hoc pairwise comparisons using multcompare
[c, m, h, gnames] = multcompare(stats, 'CType', 'dunn-sidak', ...
    'Display', 'off');
p_values = c(:, end);


%% add the statistical tests on the plot

pressure_levels = {'Low', 'Medium', 'High'};
set(gca, 'XTick', 1:n_cond, 'XTickLabel', pressure_levels);
xlh11 = xlabel('Pressure', 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold'); 


ylh11 = ylabel('Baseline-Corrected Power (dB)', ...
    'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold'); 

xlim([0.5 3.5])
title(sprintf('Significant Region'), 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold');
% grid on;
% set(gca, 'YTick', 1:10)
set(gca,'FontSize', 24); box off;

multtestStr = sprintf(['Kruskal-Wallis with\nDunn-Sidak Correction:\n\n', ...
    'Low vs. Medium:\\ast(p < 0.01)\n', ...
    'Low vs. High:\\ast(p < 0.01)\n', ...
    'Medium vs. High: ns (p = 1.0)']);
hTxt = text(0.6, -1, multtestStr, 'HorizontalAlignment', 'left', ...
    'FontSize', 18, 'FontWeight', 'normal', 'BackgroundColor', 'none', 'EdgeColor', 'none');

hold off;



%% Option 2: Power average per specral band across the significant region

sig_points = find(pcond == 1);
[row, col] = ind2sub(size(pcond), sig_points);

freqs_sig = freqs(row);

% alpha
alpha_sig = and(freqs_sig >= alpha_band(1) , freqs_sig <= alpha_band(2));
alpha_sig_indx = sig_points(alpha_sig);

reshaped_ersp_subj = cellfun(@(x) reshape(x, [], size(ersp_all{1,1}, 3)), ...
    ersp_all, 'UniformOutput', false);
ersp_subj_alpha_region = cellfun(@(x) x(alpha_sig_indx, :), reshaped_ersp_subj, ...
    'UniformOutput', false);
% beta
beta_sig = and(freqs_sig >= beta_band(1) , freqs_sig <= beta_band(2));
beta_sig_indx = sig_points(beta_sig);

reshaped_ersp_subj = cellfun(@(x) reshape(x, [], size(ersp_all{1,1}, 3)), ...
    ersp_all, 'UniformOutput', false);
ersp_subj_beta_region = cellfun(@(x) x(beta_sig_indx, :), reshaped_ersp_subj, ...
    'UniformOutput', false);


%% -- plot Option 2 ---


fig2 = figure('name', ['Power at significant region per specral band - Right Prim Motor'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape');
ax_option2_1 = subplot(1, 2, 1);
hold on;
% Jitter setup for spread
spread = 0.15;
rng(1); % Sets the random seed to 1
% x positions (jittered for clarity)
loc = 1:3;
x_jitter = loc + (rand( size(ersp_all{1,1}, 3) ,numel(loc))-0.5)*spread;

subjects_unique = 1:size(ersp_all{1,1}, 3);
n_cond = length(ersp_all);

subj_power_alpha_region = [];
for s = 1:numel(subjects_unique)
    ersp_subj_alpha_region_P1 = ersp_subj_alpha_region{1, 1}(:, s);
    ersp_subj_alpha_region_P3 = ersp_subj_alpha_region{2, 1}(:, s);
    ersp_subj_alpha_region_P6 = ersp_subj_alpha_region{3, 1}(:, s);

    subj_power_alpha_region = [subj_power_alpha_region; ...
        [mean(ersp_subj_alpha_region_P1), ...
         mean(ersp_subj_alpha_region_P3), ...
         mean(ersp_subj_alpha_region_P6)] ];
end




for s = 1:numel(subjects_unique)
    x = x_jitter(s, :);
    subj_power_s = subj_power_alpha_region(s, :);
    % 3. Paired line for within-subject
    plot(x, subj_power_s, '-', 'Color', [0.8 0.8 0.8 0.8], 'LineWidth', 0.5);
end


% 1. Box plot layer
for i = 1:n_cond
    
    data = subj_power_alpha_region(:, i);

    % boxplot(data, ...
    %     'Whisker', 2.5, ...
    %     'Positions', (i)-spread, ...
    %     'BoxStyle', 'outline', ...
    %     'Colors', colors(i, :), ...
    %     '') % Changes the outlier rule to 2.5×IQR

    % Boxplot: use boxchart (R2019b+), else use boxplot
    boxchart(ones(length(data),1)*(i)-spread, data, ...
        'BoxFaceColor', colors(i,:), ...
        'BoxFaceAlpha', 0.6, ...
        'BoxWidth', 0.30, ...
        'MarkerStyle', 'none', 'LineWidth', 2, ...
        'BoxEdgeColor', colors(i,:)*0.7, ...
        'WhiskerLineColor', colors(i,:)*0.7);
end


% 2. Individual data points (jittered scatter)
for s = 1:numel(subjects_unique)
    subj_power_s = subj_power_alpha_region(s, :);
    x = x_jitter(s, :);

    % scatter each subject's points
    for j = 1:numel(loc)
        scatter(x(j), subj_power_s(j), 100, 'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', colors(j,:), ...
            'MarkerFaceAlpha', 0.4, ...
            'MarkerEdgeAlpha', 0.4, ...
            'LineWidth',0.7);
    end
end


% 4. Overlay group means (LMM, or just mean if you prefer)
mu = mean(subj_power_alpha_region, 1);
x = [1, 2, 3] - repmat(spread, 1, 3);
plot(x, mu, 'LineWidth', 2, 'Color', [0.3 0.3 0.3])



%% Pairwise comparison tests
% Test each condition's normality separately
[h1, p1] = swtest(subj_power_alpha_region(:, 1));  % or use lillietest
[h2, p2] = swtest(subj_power_alpha_region(:, 2));
[h3, p3] = swtest(subj_power_alpha_region(:, 3));

% based on the results of normality test we should use non-parametric tests
% Kruskal-Wallis with Post-hoc Tests
group = {'Low Pressure', 'Medium Pressure', 'High Pressure'};
[p, tbl, stats] = kruskalwallis(subj_power_alpha_region, group, 'off');

% if p < 0.05 then:
% Post-hoc pairwise comparisons using multcompare
[c, m, h, gnames] = multcompare(stats, 'CType', 'dunn-sidak', ...
    'Display', 'off');
p_values = c(:, end);


%% add the statistical tests on the plot

pressure_levels = {'Low', 'Medium', 'High'};
set(gca, 'XTick', 1:n_cond, 'XTickLabel', pressure_levels);
xlh11 = xlabel('Pressure', 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold'); 


ylh11 = ylabel('Baseline-Corrected Power (dB)', ...
    'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold'); 

% ylim([-3.5 3.5]);
xlim([0.5 3.5])
title(sprintf('Mean \\alpha Power over\nSignificant Region'), 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold');
% grid on;
% set(gca, 'YTick', 1:10)
set(gca,'FontSize', 24); box off;

% multtestStr = sprintf(['Kruskal-Wallis with\nDunn-Sidak Correction:\n\n', ...
%     'Low vs. Medium:\\ast(p < 0.01)\n', ...
%     'Low vs. High:\\ast(p < 0.001)\n', ...
%     'Medium vs. High: ns (p = 0.16)']);
multtestStr = sprintf(['Low vs. Medium:\\ast(p < 0.01)\n', ...
    'Low vs. High:\\ast(p < 0.01)\n', ...
    'Medium vs. High: ns (p = 0.99)']);
hTxt = text(0.6, -1.7, multtestStr, 'HorizontalAlignment', 'left', ...
    'FontSize', 16, 'FontWeight', 'normal', 'BackgroundColor', 'none', ...
    'EdgeColor', 'none', 'FontName', 'Arial');

hold off;



%% -- plot Option 2 --- beta power over significant region
ax_option2_2 = subplot(1, 2, 2);
hold on;
% Jitter setup for spread
spread = 0.15;
rng(1); % Sets the random seed to 1
% x positions (jittered for clarity)
loc = 1:3;
x_jitter = loc + (rand( size(ersp_all{1,1}, 3) ,numel(loc))-0.5)*spread;

subjects_unique = 1:size(ersp_all{1,1}, 3);
n_cond = length(ersp_all);

subj_power_beta_region = [];
for s = 1:numel(subjects_unique)
    ersp_subj_beta_region_P1 = ersp_subj_beta_region{1, 1}(:, s);
    ersp_subj_beta_region_P3 = ersp_subj_beta_region{2, 1}(:, s);
    ersp_subj_beta_region_P6 = ersp_subj_beta_region{3, 1}(:, s);

    subj_power_beta_region = [subj_power_beta_region; ...
        [mean(ersp_subj_beta_region_P1), ...
         mean(ersp_subj_beta_region_P3), ...
         mean(ersp_subj_beta_region_P6)] ];
end




for s = 1:numel(subjects_unique)
    x = x_jitter(s, :);
    subj_power_s = subj_power_beta_region(s, :);
    % 3. Paired line for within-subject
    plot(x, subj_power_s, '-', 'Color', [0.8 0.8 0.8 0.8], 'LineWidth', 0.5);
end


% 1. Box plot layer
for i = 1:n_cond
    
    data = subj_power_beta_region(:, i);

    % boxplot(data, ...
    %     'Whisker', 2.5, ...
    %     'Positions', (i)-spread, ...
    %     'BoxStyle', 'outline', ...
    %     'Colors', colors(i, :), ...
    %     '') % Changes the outlier rule to 2.5×IQR

    % Boxplot: use boxchart (R2019b+), else use boxplot
    boxchart(ones(length(data),1)*(i)-spread, data, ...
        'BoxFaceColor', colors(i,:), ...
        'BoxFaceAlpha', 0.6, ...
        'BoxWidth', 0.30, ...
        'MarkerStyle', 'none', 'LineWidth', 2, ...
        'BoxEdgeColor', colors(i,:)*0.7, ...
        'WhiskerLineColor', colors(i,:)*0.7);
end


% 2. Individual data points (jittered scatter)
for s = 1:numel(subjects_unique)
    subj_power_s = subj_power_beta_region(s, :);
    x = x_jitter(s, :);

    % scatter each subject's points
    for j = 1:numel(loc)
        scatter(x(j), subj_power_s(j), 100, 'MarkerFaceColor', colors(j,:), ...
            'MarkerEdgeColor', colors(j,:), ...
            'MarkerFaceAlpha', 0.4, ...
            'MarkerEdgeAlpha', 0.4, ...
            'LineWidth',0.7);
    end
end


% 4. Overlay group means (LMM, or just mean if you prefer)
mu = mean(subj_power_beta_region, 1);
x = [1, 2, 3] - repmat(spread, 1, 3);
plot(x, mu, 'LineWidth', 2, 'Color', [0.3 0.3 0.3])



%% Pairwise comparison tests
% Test each condition's normality separately
[h1, p1] = swtest(subj_power_beta_region(:, 1)); 
[h2, p2] = swtest(subj_power_beta_region(:, 2));
[h3, p3] = swtest(subj_power_beta_region(:, 3));

% h1 is 0 but still h2 and h3 are 1
% based on the results of normality test we should use non-parametric tests
% Kruskal-Wallis with Post-hoc Tests
group = {'Low Pressure', 'Medium Pressure', 'High Pressure'};
[p, tbl, stats] = kruskalwallis(subj_power_beta_region, group, 'off');

% if p < 0.05 then:
% Post-hoc pairwise comparisons using multcompare
[c, m, h, gnames] = multcompare(stats, 'CType', 'dunn-sidak', ...
    'Display', 'off');
p_values = c(:, end);


%% add the statistical tests on the plot


pressure_levels = {'Low', 'Medium', 'High'};
set(gca, 'XTick', 1:n_cond, 'XTickLabel', pressure_levels);
xlh11 = xlabel('Pressure', 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold'); 


% ylh11 = ylabel('Baseline-Corrected Power (dB)', ...
%     'FontName', 'Arial', 'FontSize', 16, ...
%     'FontWeight', 'bold'); 

ylim([-2 1]); % same as the alpha power at the 1st subplot
xlim([0.5 3.5])
title(sprintf('Mean \\beta Power over\nSignificant Region'), 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold');
% grid on;
% set(gca, 'YTick', 1:10)
set(gca,'FontSize', 24); box off;

% multtestStr = sprintf(['Kruskal-Wallis with\nDunn-Sidak Correction:\n\n', ...
%     'Low vs. Medium:\\ast(p < 0.01)\n', ...
%     'Low vs. High:\\ast(p < 0.001)\n', ...
%     'Medium vs. High: ns (p = 0.16)']);
multtestStr = sprintf(['Low vs. Medium:\\ast(p < 0.01)\n', ...
    'Low vs. High:\\ast(p < 0.05)\n', ...
    'Medium vs. High: ns (p = 0.96)']);
hTxt = text(0.6, -1.7, multtestStr, 'HorizontalAlignment', 'left', ...
    'FontSize', 16, 'FontWeight', 'normal', 'BackgroundColor', 'none', ...
    'EdgeColor', 'none', 'FontName', 'Arial');

hold off;




%% Option 3: Power average per specral band across cycle percentage (over time)

theta_band = [4, 8]; % Hz
alpha_band = [8, 14]; % Hz
beta_band  = [14, 30]; % Hz
low_gamma_band = [30, 60]; % Hz
high_gamma_band = [60, 120]; % Hz


theta_indx = [find(freqs > theta_band(1), 1), ...
              find(freqs > theta_band(2), 1)-1];
alpha_indx = [find(freqs > alpha_band(1), 1), ...
              find(freqs > alpha_band(2), 1)-1];
beta_indx  = [find(freqs > beta_band(1), 1), ...
              find(freqs > beta_band(2), 1)-1];
low_gamma_indx  = [find(freqs > low_gamma_band(1), 1), ...
                   find(freqs > low_gamma_band(2), 1)-1];
high_gamma_indx  = [find(freqs > high_gamma_band(1), 1), ...
                    find(freqs > high_gamma_band(2), 1)-1];


ersp_subj_theta_time =  cellfun(@(x) x(theta_indx(1):theta_indx(2), :, :), ersp_all, ...
    'UniformOutput', false);
ersp_theta_time = cellfun(@(x) mean(x, 1), ersp_subj_theta_time, ...
    'UniformOutput', false);

ersp_subj_alpha_time =  cellfun(@(x) x(alpha_indx(1):alpha_indx(2), :, :), ersp_all, ...
    'UniformOutput', false);
ersp_alpha_time = cellfun(@(x) mean(x, 1), ersp_subj_alpha_time, ...
    'UniformOutput', false);

ersp_subj_beta_time =  cellfun(@(x) x(beta_indx(1):beta_indx(2), :, :), ersp_all, ...
    'UniformOutput', false);
ersp_beta_time = cellfun(@(x) mean(x, 1), ersp_subj_beta_time, ...
    'UniformOutput', false);

ersp_subj_low_gamma_time =  cellfun(@(x) x(low_gamma_indx(1):low_gamma_indx(2), :, :), ersp_all, ...
    'UniformOutput', false);
ersp_low_gamma_time = cellfun(@(x) mean(x, 1), ersp_subj_low_gamma_time, ...
    'UniformOutput', false);

ersp_subj_high_gamma_time =  cellfun(@(x) x(high_gamma_indx(1):high_gamma_indx(2), :, :), ersp_all, ...
    'UniformOutput', false);
ersp_high_gamma_time = cellfun(@(x) mean(x, 1), ersp_subj_high_gamma_time, ...
    'UniformOutput', false);

N = size(ersp_all{1,1}, 3);

ersp_theta_time_P1_mean = mean(ersp_theta_time{1, 1}, 3);
ersp_theta_time_P1_std = std(ersp_theta_time{1, 1}, 0, 3);
ersp_theta_time_P1_sem = std(ersp_theta_time{1, 1}, 0, 3)/sqrt(N);
ersp_theta_time_P3_mean = mean(ersp_theta_time{2, 1}, 3);
ersp_theta_time_P3_std = std(ersp_theta_time{2, 1}, 0, 3);
ersp_theta_time_P3_sem = std(ersp_theta_time{2, 1}, 0, 3)/sqrt(N);
ersp_theta_time_P6_mean = mean(ersp_theta_time{3, 1}, 3);
ersp_theta_time_P6_std = std(ersp_theta_time{3, 1}, 0, 3);
ersp_theta_time_P6_sem = std(ersp_theta_time{3, 1}, 0, 3)/sqrt(N);

ersp_alpha_time_P1_mean = mean(ersp_alpha_time{1, 1}, 3);
ersp_alpha_time_P1_std = std(ersp_alpha_time{1, 1}, 0, 3);
ersp_alpha_time_P1_sem = std(ersp_alpha_time{1, 1}, 0, 3)/sqrt(N);
ersp_alpha_time_P3_mean = mean(ersp_alpha_time{2, 1}, 3);
ersp_alpha_time_P3_std = std(ersp_alpha_time{2, 1}, 0, 3);
ersp_alpha_time_P3_sem = std(ersp_alpha_time{2, 1}, 0, 3)/sqrt(N);
ersp_alpha_time_P6_mean = mean(ersp_alpha_time{3, 1}, 3);
ersp_alpha_time_P6_std = std(ersp_alpha_time{3, 1}, 0, 3);
ersp_alpha_time_P6_sem = std(ersp_alpha_time{3, 1}, 0, 3)/sqrt(N);

ersp_beta_time_P1_mean = mean(ersp_beta_time{1, 1}, 3);
ersp_beta_time_P1_std = std(ersp_beta_time{1, 1}, 0, 3);
ersp_beta_time_P1_sem = std(ersp_beta_time{1, 1}, 0, 3)/sqrt(N);
ersp_beta_time_P3_mean = mean(ersp_beta_time{2, 1}, 3);
ersp_beta_time_P3_std = std(ersp_beta_time{2, 1}, 0, 3);
ersp_beta_time_P3_sem = std(ersp_beta_time{2, 1}, 0, 3)/sqrt(N);
ersp_beta_time_P6_mean = mean(ersp_beta_time{3, 1}, 3);
ersp_beta_time_P6_std = std(ersp_beta_time{3, 1}, 0, 3);
ersp_beta_time_P6_sem = std(ersp_beta_time{3, 1}, 0, 3)/sqrt(N);

ersp_low_gamma_time_P1_mean = mean(ersp_low_gamma_time{1, 1}, 3);
ersp_low_gamma_time_P1_std = std(ersp_low_gamma_time{1, 1}, 0, 3);
ersp_low_gamma_time_P1_sem = std(ersp_low_gamma_time{1, 1}, 0, 3)/sqrt(N);
ersp_low_gamma_time_P3_mean = mean(ersp_low_gamma_time{2, 1}, 3);
ersp_low_gamma_time_P3_std = std(ersp_low_gamma_time{2, 1}, 0, 3);
ersp_low_gamma_time_P3_sem = std(ersp_low_gamma_time{2, 1}, 0, 3)/sqrt(N);
ersp_low_gamma_time_P6_mean = mean(ersp_low_gamma_time{3, 1}, 3);
ersp_low_gamma_time_P6_std = std(ersp_low_gamma_time{3, 1}, 0, 3);
ersp_low_gamma_time_P6_sem = std(ersp_low_gamma_time{3, 1}, 0, 3)/sqrt(N);

ersp_high_gamma_time_P1_mean = mean(ersp_high_gamma_time{1, 1}, 3);
ersp_high_gamma_time_P1_std = std(ersp_high_gamma_time{1, 1}, 0, 3);
ersp_high_gamma_time_P1_sem = std(ersp_high_gamma_time{1, 1}, 0, 3)/sqrt(N);
ersp_high_gamma_time_P3_mean = mean(ersp_high_gamma_time{2, 1}, 3);
ersp_high_gamma_time_P3_std = std(ersp_high_gamma_time{2, 1}, 0, 3);
ersp_high_gamma_time_P3_sem = std(ersp_high_gamma_time{2, 1}, 0, 3)/sqrt(N);
ersp_high_gamma_time_P6_mean = mean(ersp_high_gamma_time{3, 1}, 3);
ersp_high_gamma_time_P6_std = std(ersp_high_gamma_time{3, 1}, 0, 3);
ersp_high_gamma_time_P6_sem = std(ersp_high_gamma_time{3, 1}, 0, 3)/sqrt(N);


%% -- plot Option 3 --- power average over all specral band across time
% Get monitor information
monitors = get(0, 'MonitorPositions');

fig3 = figure('name', ['Power average over specral bands across time - Right Prim Motor'], ...
    'InvertHardcopy', 'off', 'PaperType', 'a2', ...
    'PaperOrientation', 'landscape', ...
    'Resize', 'off');


% For second monitor (row 2), add drawnow before setting position
drawnow;  % Let MATLAB finish drawing on primary monitor first
pause(0.1);  % Short pause helps

% Now set position on second monitor
% set(gcf, 'Position', monitors(1, :));  % Use entire second monitor

set(fig3, 'Position', [monitors(1,1)-150, monitors(1,2)+700, 2000, 600]);  % Use entire second monitor


%% theta band
ax_option3_1 = subplot(1, 5, 1); hold on;

fill([times fliplr(times)], ...
    [ersp_theta_time_P1_mean + ersp_theta_time_P1_sem, ...
     fliplr(ersp_theta_time_P1_mean - ersp_theta_time_P1_sem)], ...
     colors(1, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_theta_time_P1_mean, ...
    'Color', colors(1, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_theta_time_P3_mean + ersp_theta_time_P3_sem, ...
     fliplr(ersp_theta_time_P3_mean - ersp_theta_time_P3_sem)], ...
     colors(2, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_theta_time_P3_mean, ...
    'Color', colors(2, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_theta_time_P6_mean + ersp_theta_time_P6_sem, ...
     fliplr(ersp_theta_time_P6_mean - ersp_theta_time_P6_sem)], ...
     colors(3, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_theta_time_P6_mean, ...
    'Color', colors(3, :), 'LineWidth', 4)

xlim([times(1) times(end)])
set(ax_option3_1, 'XTick', [times(1), times(66), times(end)], ...
    'XTickLabel', {'0', '50', '100'}, 'XTickLabelRotation', 45);

xlh = xlabel('Cycle (%)', 'FontName', 'Arial', 'FontWeight', 'bold');

tth1 = title('$\theta$ band', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

ylh = ylabel(sprintf('Baseline-Corrected\nPower (dB)'), ...
    'FontName', 'Arial', 'FontWeight', 'bold');
ylh.Position(1) = ylh.Position(1) - 400;

set(ax_option3_1, 'FontSize', 20, 'Box', 'on')




%% alpha band
ax_option3_2 = subplot(1, 5, 2); hold on;

fill([times fliplr(times)], ...
    [ersp_alpha_time_P1_mean + ersp_alpha_time_P1_sem, ...
     fliplr(ersp_alpha_time_P1_mean - ersp_alpha_time_P1_sem)], ...
     colors(1, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_alpha_time_P1_mean, ...
    'Color', colors(1, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_alpha_time_P3_mean + ersp_alpha_time_P3_sem, ...
     fliplr(ersp_alpha_time_P3_mean - ersp_alpha_time_P3_sem)], ...
     colors(2, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_alpha_time_P3_mean, ...
    'Color', colors(2, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_alpha_time_P6_mean + ersp_alpha_time_P6_sem, ...
     fliplr(ersp_alpha_time_P6_mean - ersp_alpha_time_P6_sem)], ...
     colors(3, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_alpha_time_P6_mean, ...
    'Color', colors(3, :), 'LineWidth', 4)

xlim([times(1) times(end)])
set(ax_option3_2, 'XTick', [times(1), times(66), times(end)], ...
    'XTickLabel', {'0', '50', '100'}, 'XTickLabelRotation', 45);

xlh = xlabel('Cycle (%)', 'FontName', 'Arial', 'FontWeight', 'bold');
% ylh = ylabel('Baseline-Corrected Power (dB)', ...
%     'FontName', 'Arial', 'FontWeight', 'bold');
tth2 = title('$\alpha$ band', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

set(ax_option3_2, 'FontSize', 20, 'Box', 'on')



%% beta band
ax_option3_3 = subplot(1, 5, 3); hold on;

fill([times fliplr(times)], ...
    [ersp_beta_time_P1_mean + ersp_beta_time_P1_sem, ...
     fliplr(ersp_beta_time_P1_mean - ersp_beta_time_P1_sem)], ...
     colors(1, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_beta_time_P1_mean, ...
    'Color', colors(1, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_beta_time_P3_mean + ersp_beta_time_P3_sem, ...
     fliplr(ersp_beta_time_P3_mean - ersp_beta_time_P3_sem)], ...
     colors(2, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_beta_time_P3_mean, ...
    'Color', colors(2, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_beta_time_P6_mean + ersp_beta_time_P6_sem, ...
     fliplr(ersp_beta_time_P6_mean - ersp_beta_time_P6_sem)], ...
     colors(3, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_beta_time_P6_mean, ...
    'Color', colors(3, :), 'LineWidth', 4)

xlim([times(1) times(end)])
set(ax_option3_3, 'XTick', [times(1), times(66), times(end)], ...
    'XTickLabel', {'0', '50', '100'}, 'XTickLabelRotation', 45);

xlh = xlabel('Cycle (%)', 'FontName', 'Arial', 'FontWeight', 'bold');
% ylh = ylabel('Baseline-Corrected Power (dB)', ...
%     'FontName', 'Arial', 'FontWeight', 'bold');
tth3 = title('$\beta$ band', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

set(ax_option3_3, 'FontSize', 20, 'Box', 'on')



%% low_gamma band
ax_option3_4 = subplot(1, 5, 4); hold on;

fill([times fliplr(times)], ...
    [ersp_low_gamma_time_P1_mean + ersp_low_gamma_time_P1_sem, ...
     fliplr(ersp_low_gamma_time_P1_mean - ersp_low_gamma_time_P1_sem)], ...
     colors(1, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_low_gamma_time_P1_mean, ...
    'Color', colors(1, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_low_gamma_time_P3_mean + ersp_low_gamma_time_P3_sem, ...
     fliplr(ersp_low_gamma_time_P3_mean - ersp_low_gamma_time_P3_sem)], ...
     colors(2, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_low_gamma_time_P3_mean, ...
    'Color', colors(2, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_low_gamma_time_P6_mean + ersp_low_gamma_time_P6_sem, ...
     fliplr(ersp_low_gamma_time_P6_mean - ersp_low_gamma_time_P6_sem)], ...
     colors(3, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_low_gamma_time_P6_mean, ...
    'Color', colors(3, :), 'LineWidth', 4)

xlim([times(1) times(end)])
set(ax_option3_4, 'XTick', [times(1), times(66), times(end)], ...
    'XTickLabel', {'0', '50', '100'}, 'XTickLabelRotation', 45);

xlh = xlabel('Cycle (%)', 'FontName', 'Arial', 'FontWeight', 'bold');
% ylh = ylabel('Baseline-Corrected Power (dB)', ...
%     'FontName', 'Arial', 'FontWeight', 'bold');
tth4 = title('Low $\gamma$ band', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

set(ax_option3_4, 'FontSize', 20, 'Box', 'on')



%% high_gamma band
ax_option3_5 = subplot(1, 5, 5); hold on;

fill([times fliplr(times)], ...
    [ersp_high_gamma_time_P1_mean + ersp_high_gamma_time_P1_sem, ...
     fliplr(ersp_high_gamma_time_P1_mean - ersp_high_gamma_time_P1_sem)], ...
     colors(1, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_high_gamma_time_P1_mean, ...
    'Color', colors(1, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_high_gamma_time_P3_mean + ersp_high_gamma_time_P3_sem, ...
     fliplr(ersp_high_gamma_time_P3_mean - ersp_high_gamma_time_P3_sem)], ...
     colors(2, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_high_gamma_time_P3_mean, ...
    'Color', colors(2, :), 'LineWidth', 4)

fill([times fliplr(times)], ...
    [ersp_high_gamma_time_P6_mean + ersp_high_gamma_time_P6_sem, ...
     fliplr(ersp_high_gamma_time_P6_mean - ersp_high_gamma_time_P6_sem)], ...
     colors(3, :), 'EdgeColor','none', 'FaceAlpha', 0.3);
plot(times, ersp_high_gamma_time_P6_mean, ...
    'Color', colors(3, :), 'LineWidth', 4)

xlim([times(1) times(end)])
set(ax_option3_5, 'XTick', [times(1), times(66), times(end)], ...
    'XTickLabel', {'0', '50', '100'}, 'XTickLabelRotation', 45);

xlh = xlabel('Cycle (%)', 'FontName', 'Arial', 'FontWeight', 'bold');
% ylh = ylabel('Baseline-Corrected Power (dB)', ...
%     'FontName', 'Arial', 'FontWeight', 'bold');
tth5 = title('High $\gamma$ band', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

set(ax_option3_5, 'FontSize', 20, 'Box', 'on')


%% set a unified YLim at the end
ylimit = [-0.6 0.6];
set(ax_option3_1, 'YLim', ylimit)
set(ax_option3_2, 'YLim', ylimit)
set(ax_option3_3, 'YLim', ylimit)
set(ax_option3_4, 'YLim', ylimit)
set(ax_option3_5, 'YLim', ylimit)


%% change the size
q = 0.8;

pos = get(ax_option3_1, 'Position');
set(ax_option3_1, 'Position', [pos(1), pos(2)+0.05, pos(3), pos(4)*q])

pos = get(ax_option3_2, 'Position');
set(ax_option3_2, 'Position', [pos(1), pos(2)+0.05, pos(3), pos(4)*q])

pos = get(ax_option3_3, 'Position');
set(ax_option3_3, 'Position', [pos(1), pos(2)+0.05, pos(3), pos(4)*q])

pos = get(ax_option3_4, 'Position');
set(ax_option3_4, 'Position', [pos(1), pos(2)+0.05, pos(3), pos(4)*q])

pos = get(ax_option3_5, 'Position');
set(ax_option3_5, 'Position', [pos(1), pos(2)+0.05, pos(3), pos(4)*q])


%% change the title position
D = 0.2;
tth1.Position(2) = tth1.Position(2) + D; 
tth2.Position(2) = tth2.Position(2) + D; 
tth3.Position(2) = tth3.Position(2) + D; 
tth4.Position(2) = tth4.Position(2) + D; 
tth5.Position(2) = tth5.Position(2) + D; 


%% Add event labels
all_axes = {ax_option3_1, ax_option3_2, ax_option3_3, ax_option3_4, ax_option3_5};
for i = 1:5
axes(all_axes{1, i});

evPlotLines_correct = [times(1) times(66) times(end)];
eventLabels_new = {sprintf('FlxS'), sprintf('FlxE\nExtS'), sprintf('ExtE')};
% add event lines from time warp
for L = 1:length(evPlotLines_correct)
    if L == 1 || L == length(evPlotLines_correct)
        %v = vline(evPlotLines(L),'-k',eventLabels{1,L},[0.05 1.05]); set(v,'LineWidth',1); %solid line
        v = vline(evPlotLines_correct(L),'-k', eventLabels_new{1,L}); set(v,'LineWidth',1); %solid line
    else
        %v = vline(evPlotLines(L),':k',eventLabels{1,L},[0.05 1.05]); set(v,'LineWidth',1.2);
        v = vline(evPlotLines_correct(L),':k',eventLabels_new{1,L}); set(v,'LineWidth',1.2);
    end
end

H = findobj(gcf);
tb = findobj(H,'Type','text');

K = 0.65
for textbox = 1:3 % 1:size(tb,1)
    if     mod(textbox, 3) == 1
        pos = tb(textbox).Position;
        tb(textbox).Position = [pos(1)+30 K 0];
        set(tb(textbox),'Rotation',90) % rotate 90 degrees
        set(tb(textbox),'FontSize',10, 'FontWeight', 'normal') 
    elseif mod(textbox, 3) == 2
        pos = tb(textbox).Position;
        tb(textbox).Position = [pos(1)-10 K 0];
        set(tb(textbox),'Rotation',90) % rotate 90 degrees
        set(tb(textbox),'FontSize',10, 'FontWeight', 'normal') 
    elseif mod(textbox, 3) == 0
        pos = tb(textbox).Position;
        tb(textbox).Position = [pos(1)+15 K 0];
        set(tb(textbox),'Rotation',90) % rotate 90 degrees
        set(tb(textbox),'FontSize',10, 'FontWeight', 'normal') 
    end
    % pos = tb(textbox).Position;
    % tb(textbox).Position = [pos(1) 150 0];
    % set(tb(textbox),'Rotation',90) % rotate 90 degrees
    % set(tb(textbox),'FontSize',8) 
end

end


% add legend
axes(ax_option3_5)
lgd = legend({'', ' Low', '', ' Medium', '', ' High'});
lgd.Location = "southeast";
lgd.FontSize = 16;
lgd.Box = "off";




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

studyNames = {'Left Dorsal ACC', 'Left Parieto Occipital', ...
    'Left PreMot SuppMot', 'Left Prim Motor', 'Prime Visual', ...
    'Right Parieto Occipital', 'Right PreMot SuppMot', 'Right Prim Motor'};

for cluster = 1:length(studyNames)

    cluster_name = studyNames{cluster};
    ersp_structure = load([cluster_name, ' ersp_results.mat']);
    ersp_structure = ersp_structure.ersp_results;
    
    ersp_all = ersp_structure.data.allersp;
    
    times = ersp_structure.data.alltimes;
    freqs = ersp_structure.data.allfreqs;
    
    %% Look at the individual Time-Frequency plots in the brain cluster
    
    Subs_ICs = load([Subjects_ICs_in_cluster_path, ...
        'Subjects_ICs_in_clusters.mat']);
    Subs_ICs = Subs_ICs.SUBJECTS_ICS;
    Subs_ICs(:, 1) = cellfun(@(x) strrep(x, '_', ' '), Subs_ICs(:, 1), ...
        'UniformOutput', false);
    
    
    idx = strcmp(Subs_ICs(:, 1), cluster_name);
    all_subjects = Subs_ICs(idx, 2);
    [all_subjects_sorted, sort_indx] = sort(all_subjects{1}.Subjects);
    all_ics = Subs_ICs(idx, 2);
    all_ics_sorted = all_ics{1}.ICs(sort_indx);
    
    monitors = get(0, 'MonitorPositions');
    fig = figure('name',['Individual Time-Frequency Plots - ' ...
        cluster_name], ...
        'InvertHardcopy', 'off', 'PaperType', 'a2', ...
        'PaperOrientation', 'landscape', ...
        'Resize', 'off');
    
    % For second monitor (row 2), add drawnow before setting position
    drawnow;  % Let MATLAB finish drawing on primary monitor first
    pause(0.1);  % Short pause helps
    
    % Now set position on second monitor
    % set(gcf, 'Position', monitors(1, :));  % Use entire second monitor
    
    numSubjects = size(ersp_all{1, 1}, 3);
    set(fig, 'Position', ...
        [monitors(1,1)+50, monitors(1,2)+400, 0.9*175*numSubjects, 960]);  
    
    ax = gobjects(1, 3*numSubjects);  % Preallocate axes array
    ttlh = gobjects(1, numSubjects);  % Preallocate title array
    
    % data = [];
    % for p = 1:3
    %     data = [data, reshape(mean(ersp_all{p, 1}, 3).',1,[])];
    % end
    % IQR = iqr(data); % interquartile range
    % Q1 = quantile(data,0.25);
    % myMin = round(Q1-1.5*IQR,1);
    % erspdata_clim = [myMin myMin*(-1)];
    
    pressure_conditions = {'Low', 'Medium', 'High'};
    for p = 1:3
    
        for sub = 1:numSubjects
            ax((p-1)*numSubjects + sub) = ...
                subplot(3, numSubjects, (p-1)*numSubjects + sub); 
            hold on;
            
    
            data = [];
            data = [data, reshape(ersp_all{p, 1}(:, :, sub).', 1,[])];
            IQR = iqr(data); % interquartile range
            Q1 = quantile(data,0.25);
            myMin = round(Q1-1.5*IQR,1);
            erspdata_clim = [myMin myMin*(-1)];
    
            contourf(times, freqs, ersp_all{p,1}(:, :, sub), 50, ...
                'linecolor','none');
            xline(times(66), 'LineWidth', 1, 'LineStyle', '--')
    
            set(gca, 'clim', erspdata_clim, 'xlim', [times(1) times(end)], ...
                'ydir', 'norm', 'ylim', [freqs(1) freqs(end)], ...
                'yscale','log')
    
            
            set(gca,'XTick',[times(1) times(66) times(end)],...
                        'XTickLabel',{'0', '50', '100'}, ...
                        'ytick', [4 8 14 30 60 120], 'fontsize', 12);
            xtickangle(45)
    
            pos = get(gca, 'Position');
            % new_pos = [pos(1), pos(2), pos(3), pos(4)*0.9];
            if p == 1 
                new_pos = [pos(1), pos(2)-0.05, pos(3), pos(4)*0.9];
            elseif p == 2
                new_pos = [pos(1), pos(2), pos(3), pos(4)*0.9];
            else
                new_pos = [pos(1), pos(2)+0.05, pos(3), pos(4)*0.9];
            end
            set(gca, 'Position', new_pos);
    
            % ylabel
            if sub == 1
                ylh = ylabel(sprintf(['\\it', pressure_conditions{p}, ' Pressure\n\n', ...
                    'Frequency (Hz)']), 'fontsize', 14, 'fontweight', 'bold', ...
                    'FontName','Arial');
                ylh.Position(1) = -1000; % it was -450 I changed it!
            else
                set(gca,'YTickLabel',[]);
            end
    
            % xlabel
            if p == 3
                xlh = xlabel('Cycle (%)', 'Fontsize', 14, 'fontweight', 'bold');
                xlh.Position(2) = 1.3;
            end
    
            % % title
            % if p == 1
            %     ttlh(sub) = title(['S', num2str(all_subjects_sorted(sub) + 4), ...
            %         ', IC', num2str(all_ics_sorted(sub))]);
            % end
                    
            % drawnow
            % pause(0.1)
        
        end
    end
    
    set(gcf,'Colormap', calldefinedcolormap(), 'Color',[1 1 1]);
    
    
    % set the same color limit for each subject across all conditions
    climits = zeros(1, numSubjects);
    for sub = 1:numSubjects
        data = [];
        for p = 1:3
            data = [data, reshape(ersp_all{p, 1}(:, :, sub).', 1,[])];
        end
        IQR = iqr(data); % interquartile range
        Q1 = quantile(data,0.25);
        Q3 = quantile(data,0.75);
        myMin = round(Q1-1.5*IQR,1);
        myMax = round(Q3+1.5*IQR,1);
        myLim = max(abs([myMin, myMax]));
        climits(1, sub) = myLim;
    
        for p = 1:3
            set(ax((p-1)*numSubjects + sub), 'CLim', [myMin myMin*(-1)]);
        end
    
    
        p = 1;
        originalPos = get(ax(sub), 'Position');
        c = colorbar(ax(sub), "northoutside");
        c.Position(2) = originalPos(2) + originalPos(4)+0.007;
        ylabel(c, '(dB)');
        c.Ticks = [-climits(sub), 0, climits(sub)];
        % caxis = c.Axis;
        % c.TickLabels
        c.Position(4) = c.Position(4)*0.7;
        c.LineWidth = 0.5;
        ax(sub).Position(4) = originalPos(4);
    
        % title
        axes(ax(sub));
        ttlh(sub) = title(['S', num2str(all_subjects_sorted(sub) + 4), ...
            ', IC', num2str(all_ics_sorted(sub))]);
         
        ttlh(sub).Position(2) = 550;
    
    end
    
    
    figname = cluster_name;
    savePath = Detailed_Analysis_on_TF_path;
    savethisfig(gcf, strcat(figname,' all subjects.png'), ...
        [savePath, '\png'],'png')
    savethisfig(gcf, strcat(figname,' all subjects.svg'), ...
        [savePath, '\svg'],'svg')
    
    cd(savePath)
    
    % savethisfig(gcf, strcat(figname,'.fig'), ...
    %     [savePath,'\ERSP\',myplotParams.figname,'\Cond_vs_Baseline\fig'],'fig')
    % savethisfig(gcf, strcat(figname,'.svg'), ...
    %     [savePath,'\ERSP\',myplotParams.figname,'\Cond_vs_Baseline\svg'],'svg')
    
    close all

end