clc
clear

%% Add and Define Necessary Paths
main_project_folder = 'C:\Morteza\MyProjects\ANSYMB2024';
addpath(genpath(main_project_folder)); % main folder containing all codes and data

data_path = 'C:\Morteza\MyProjects\ANSYMB2024\data\';
epoched_data_path = [data_path, '6_Trials_Info_and_Epoched_data\'];


%% load and plot the scores per subject
subject_list = 5:18;


scores = cell(size(subject_list));
pressures = cell(size(subject_list));

R1 = cell(size(subject_list));
R2 = cell(size(subject_list));

for i = 1:length(subject_list)
    filepath = ['sub-', num2str(subject_list(i)), '\Trials_Info.mat'];
    Trials_Info = load(fullfile(epoched_data_path, filepath));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});


    subject = subject_list(i);
    S = [];
    for j = 1:length(Trials_Info)
        if subject >= 10
           if strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
                S = cat(2, S, Trials_Info{1, j}.General.Score);
           end
        else
            S = cat(2, S, Trials_Info{1, j}.General.Score);
        end
    end

    R1{1, i} = unifrnd(-1, 1, size(S(S ~= 0)));
    R2{1, i} = unifrnd(-1, 1, size(S(S ~= 0)));

    scores{1, i} = S(S ~= 0);
    
end


%% Plot the scores
figure()
hold on

for i = 1:length(subject_list)
    plot(scores{1, i} + 0.3*R1{1, i}, i*ones(size(scores{1, i})) + 0.2*R2{1, i}, ...
        'Marker', '*', 'LineStyle', 'none', 'MarkerSize', 4, 'MarkerEdgeColor', 'k')
end
set(gca, 'XTick', 1:10, 'XTickLabel', 1:10);
xlim([0 11])
set(gca, 'YTick', 1:14, 'YTickLabel', 5:18);
ylim([0 15])
xlabel('Subject Scores')
ylabel('Subjects')
set(gca, 'FontSize', 16)
set(gcf, "Position", [100 100 700 500])


%% Plot clustered scores with KNN

figure()
hold on

num_clusters = 3; % Define number of clusters
colors = {'r', 'g', 'b'}; % Define cluster colors

for i = 1:length(subject_list)
    max_score = max(scores{1, i});
    min_score = min(scores{1, i});
    
    initial_centers = [min_score, (min_score + max_score)/2, max_score]';
    
    if subject_list(i) == 6
        initial_centers = [1, 2, 7]';
    elseif subject_list(i) == 7
        initial_centers = [1, 1.5, 5]';
    elseif subject_list(i) == 8
        initial_centers = [1, 5, 8]';
    elseif subject_list(i) == 9
        initial_centers = [1, 5, 6]';
    elseif subject_list(i) == 14
        initial_centers = [1, 3, 8]';
    elseif subject_list(i) == 16
        initial_centers = [1, 3, 8]';
    elseif subject_list(i) == 17
        initial_centers = [1, 3, 8]';    
    end
    
    % Apply k-means clustering per subject
    % [idx, C] = kmeans(scores{1, i}' + R1{1, i}', num_clusters, 'Start', initial_centers);
    [idx, C] = kmeans(scores{1, i}', num_clusters, 'Start', initial_centers);

    
    % Sort cluster centers and get sorted indices
    [~, sorted_indices] = sort(C);
    
    % Assign clusters based on sorted order
    new_idx = zeros(size(idx));
    for j = 1:num_clusters
        new_idx(idx == sorted_indices(j)) = j;
    end
    
    % Plot each cluster with aligned colors
    for j = 1:num_clusters
        cluster_points = scores{1, i}(new_idx == j);
        plot(cluster_points + 0.3 * R1{1, i}(new_idx == j), ...
            i * ones(size(cluster_points)) + 0.2 * R2{1, i}(new_idx == j), ...
            '*', 'MarkerSize', 4, 'MarkerEdgeColor', colors{j});
    end
end
set(gca, 'XTick', 1:10, 'XTickLabel', 1:10);
xlim([0 11])
set(gca, 'YTick', 1:14, 'YTickLabel', 5:18);
ylim([0 15])
xlabel('Subject Scores')
ylabel('Subjects')
title('Clustered using KNN')
set(gca, 'FontSize', 16)
set(gcf, "Position", [100 100 700 500])


%% Clustering the score-pressure matrix (pressure-informed clustering)

scores_pressure = cell(size(subject_list));

for i = 1:length(subject_list)
    filepath = ['sub-', num2str(subject_list(i)), '\Trials_Info.mat'];
    Trials_Info = load(fullfile(epoched_data_path, filepath));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});


    subject = subject_list(i);
    S = [];
    for j = 1:length(Trials_Info)
        
        p = Trials_Info{1, j}.General.Pressure;
        if p == 3
            p = 1.2;
        elseif p == 6
            p = 1.4;
        end

        if subject >= 10
           if strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
                S = cat(2, S, [Trials_Info{1, j}.General.Score; p]);
           end
        else
            S = cat(2, S, [Trials_Info{1, j}.General.Score; p]);
        end

    end

    

    scores_pressure{1, i} = S(:, S(1,:) ~= 0);
end


%% Plot clusters (pressure-informed clustering)


num_clusters = 3; % Define number of clusters
colors = {'r', 'g', 'b'}; % Define cluster colors

for i = 1:length(subject_list)
    
    max_score = max(scores_pressure{1, i}(1, :));
    min_score = min(scores_pressure{1, i}(1, :));
    
    initial_centers = [min_score, (min_score + max_score)/2, max_score; 1, 1.2, 1.4]';
    
    if subject_list(i) == 6
        initial_centers = [1, 2, 7; 1, 1.2, 1.4]';
    elseif subject_list(i) == 7
        initial_centers = [1, 1.5, 5; 1, 1.2, 1.4]';
    elseif subject_list(i) == 8
        initial_centers = [1, 5, 8; 1, 1.2, 1.4]';
    elseif subject_list(i) == 9
        initial_centers = [1, 5, 6; 1, 1.2, 1.4]';
    elseif subject_list(i) == 14
        initial_centers = [1, 3, 8; 1, 1.2, 1.4]';
    elseif subject_list(i) == 16
        initial_centers = [1, 3, 8; 1, 1.2, 1.4]';
    elseif subject_list(i) == 17
        initial_centers = [1, 3, 8; 1, 1.2, 1.4]';    
    end
    
    % Apply k-means clustering per subject
    % [idx, C] = kmeans(scores{1, i}' + R1{1, i}', num_clusters, 'Start', initial_centers);
    [idx, C] = kmeans(scores_pressure{1, i}', num_clusters, 'Start', initial_centers);

    
    % Sort cluster centers and get sorted indices
    [~, sorted_indices] = sort(C(:, 1));
    
    % Assign clusters based on sorted order
    new_idx = zeros(size(idx));
    for j = 1:num_clusters
        new_idx(idx == sorted_indices(j)) = j;
    end
    

    clusters = struct('S1', [], 'S2', [], 'S3', []);
    clusters.S1 = zeros(3, 10);
    clusters.S2 = zeros(3, 10);
    clusters.S3 = zeros(3, 10);
    for j = 1:10
        % S1
        c1 = new_idx == 1;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 1;
        clusters.S1(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S1(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.4;
        clusters.S1(3, j) = sum(and(and(c1', c2), c3));


        % S2
        c1 = new_idx == 2;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 1;
        clusters.S2(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S2(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.4;
        clusters.S2(3, j) = sum(and(and(c1', c2), c3));


        % S3
        c1 = new_idx == 3;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 1;
        clusters.S3(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S3(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.4;
        clusters.S3(3, j) = sum(and(and(c1', c2), c3));
    end


    % Create the bubble chart
    figure;
    hold on; 

    % X (scores), Y (conditions), and Size data for scatter plot + color (based on clustering)
    [X, Y] = meshgrid(1:10, [1, 3, 6]); % 10 classes (x-axis), 3 conditions (y-axis)
    X = X(:);
    Y = Y(:);

    sizes = clusters.S1;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end


    sizes = clusters.S2;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end


    sizes = clusters.S3;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end

    % Customize the plot
    xlabel('Scores');
    title(['Subject ', num2str(subject_list(i)),' - Count of Pressure Conditions per Scores']);
    yticks([1, 3, 6]);
    xticks(1:10)
    yticklabels({'P1', 'P3', 'P6'});
    set(gca, 'FontSize', 14)
    grid on

    

    % Adjust axis limits for better visibility
    xlim([-0.5, 10.5]);
    ylim([-2, 8]);
    set(gcf, "Position", [300, 200, 800, 400])


    
end





%% Clustering the score-pressure matrix (pressure-informed clustering)
% despite the previous section, here pressures are not assigned with
% constant distances (1, 1.2, 1.4). We calculate the mean score at each
% pressure and then compute the ratios of the scores average. 

scores_pressure = cell(size(subject_list));

for i = 1:length(subject_list)
    filepath = ['sub-', num2str(subject_list(i)), '\Trials_Info.mat'];
    Trials_Info = load(fullfile(epoched_data_path, filepath));
    name = fieldnames(Trials_Info);
    Trials_Info = Trials_Info.(name{1});


    subject = subject_list(i);
    S = [];
    % P = [];
    for j = 1:length(Trials_Info)

        pressure = Trials_Info{1, j}.General.Pressure;
        if pressure == 1
            p = 0.2;
        elseif pressure == 3
            p = 0.6;
        elseif pressure == 6
            p = 1.2;
        end
           

        
        if subject >= 10
           if strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
                S = cat(2, S, [Trials_Info{1, j}.General.Score; p]);
                % P = cat(2, P, Trials_Info{1, j}.General.Pressure);
           end
        else
            S = cat(2, S, [Trials_Info{1, j}.General.Score; p]);
            % P = cat(2, P, Trials_Info{1, j}.General.Pressure);
        end

    end

    % S = S(S ~= 0);
    % P = P(S ~= 0);
    % ratio3_1 = mean(S(P == 3))/mean(S(P == 1));
    % ratio6_1 = mean(S(P == 6))/mean(S(P == 1));
    

    scores_pressure{1, i} = S(:, S(1,:) ~= 0);
    % pressures{1, i} = P(S ~= 0);
end



%% Plot the results
num_clusters = 3; % Define number of clusters
colors = {'r', 'g', 'b'}; % Define cluster colors

for i = 1:length(subject_list)
    
    max_score = max(scores_pressure{1, i}(1, :));
    min_score = min(scores_pressure{1, i}(1, :));
    
    initial_centers = [min_score, (min_score + max_score)/2, max_score; 0.2, 0.6, 1.2]';
    
    if subject_list(i) == 6
        initial_centers = [1, 2, 7; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 7
        initial_centers = [1, 1.5, 5; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 8
        initial_centers = [1, 5, 8; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 9
        initial_centers = [1, 5, 6; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 14
        initial_centers = [1, 3, 8; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 16
        initial_centers = [1, 3, 8; 0.2, 0.6, 1.2]';
    elseif subject_list(i) == 17
        initial_centers = [1, 3, 8; 0.2, 0.6, 1.2]';    
    end
    
    % Apply k-means clustering per subject
    % [idx, C] = kmeans(scores{1, i}' + R1{1, i}', num_clusters, 'Start', initial_centers);
    [idx, C] = kmeans(scores_pressure{1, i}', num_clusters, 'Start', initial_centers);

    
    % Sort cluster centers and get sorted indices
    [~, sorted_indices] = sort(C(:, 1));
    
    % Assign clusters based on sorted order
    new_idx = zeros(size(idx));
    for j = 1:num_clusters
        new_idx(idx == sorted_indices(j)) = j;
    end
    

    clusters = struct('S1', [], 'S2', [], 'S3', []);
    clusters.S1 = zeros(3, 10);
    clusters.S2 = zeros(3, 10);
    clusters.S3 = zeros(3, 10);
    for j = 1:10
        % S1
        c1 = new_idx == 1;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 0.2;
        clusters.S1(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 0.6;
        clusters.S1(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S1(3, j) = sum(and(and(c1', c2), c3));


        % S2
        c1 = new_idx == 2;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 0.2;
        clusters.S2(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 0.6;
        clusters.S2(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S2(3, j) = sum(and(and(c1', c2), c3));


        % S3
        c1 = new_idx == 3;
        c2 = scores_pressure{1, i}(1, :) == j;
        c3 = scores_pressure{1, i}(2, :) == 0.2;
        clusters.S3(1, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 0.6;
        clusters.S3(2, j) = sum(and(and(c1', c2), c3));

        c3 = scores_pressure{1, i}(2, :) == 1.2;
        clusters.S3(3, j) = sum(and(and(c1', c2), c3));
    end


    % Create the bubble chart
    figure;
    hold on; 

    % X (scores), Y (conditions), and Size data for scatter plot + color (based on clustering)
    [X, Y] = meshgrid(1:10, [1, 3, 6]); % 10 classes (x-axis), 3 conditions (y-axis)
    X = X(:);
    Y = Y(:);

    sizes = clusters.S1;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end


    sizes = clusters.S2;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end


    sizes = clusters.S3;         % Transpose to align with Y (conditions)
    sizes = sizes(:) * 100; % Flatten and scale for visibility in scatter
    sizes(sizes == 0) = NaN;
    scatter(X, Y, sizes, 'filled', 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
    for k = 1:length(X)
        if sizes(k) > 0
            text(X(k), Y(k), sprintf('%d', sizes(k)/100), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontWeight', 'bold');
        end
    end

    % Customize the plot
    xlabel('Scores');
    title(['Subject ', num2str(subject_list(i)),' - Count of Pressure Conditions per Scores']);
    yticks([1, 3, 6]);
    xticks(1:10)
    yticklabels({'P1', 'P3', 'P6'});
    set(gca, 'FontSize', 14)
    grid on

    

    % Adjust axis limits for better visibility
    xlim([-0.5, 10.5]);
    ylim([-2, 8]);
    set(gcf, "Position", [300, 200, 800, 400])


    
end
















%% plot clustered scores using Gaussian Mixture Model GMM

% figure()
% hold on
% 
% num_clusters = 3; % Define number of clusters
% colors = {'r', 'g', 'b'}; % Define cluster colors
% 
% for i = 1:length(subject_list)
%     valid_scores = scores{1, i}(scores{1, i} ~= 0);  % Extract non-zero scores
% 
%     % Fit Gaussian Mixture Model (GMM) with 3 clusters
%     gm = fitgmdist(valid_scores' + R1{1, i}', num_clusters);%, 'RegularizationValue', 1e-6);%, 'SharedCovariance', true);
% 
%     % Assign cluster indices using posterior probabilities
%     idx = cluster(gm, valid_scores' + R1{1, i}');
% 
%     % Sort cluster centers and get sorted indices
%     [~, sorted_indices] = sort(gm.mu); 
% 
%     % Assign clusters based on sorted order
%     new_idx = zeros(size(idx));
%     for j = 1:num_clusters
%         new_idx(idx == sorted_indices(j)) = j;
%     end
% 
%     % Plot each cluster with aligned colors
%     for j = 1:num_clusters
%         cluster_points = valid_scores(new_idx == j);
%         plot(cluster_points + 0.3 * R1{1, i}(new_idx == j), ...
%             i * ones(size(cluster_points)) + 0.2 * R2{1, i}(new_idx == j), ...
%             '*', 'MarkerSize', 4, 'MarkerEdgeColor', colors{j});
%     end
% end
% set(gca, 'XTick', 1:10, 'XTickLabel', 1:10);
% xlim([0 11])
% set(gca, 'YTick', 1:14, 'YTickLabel', 5:18);
% ylim([0 15])
% xlabel('Subject Scores')
% ylabel('Subjects')
% set(gca, 'FontSize', 16)
% set(gcf, "Position", [100 100 700 500])



