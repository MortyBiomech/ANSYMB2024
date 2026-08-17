function [new_idx, trials] = Score_Clustering(Trials_Info, subject)

    
    scores = [];
    trials = [];
    for j = 1:length(Trials_Info)
        if subject >= 10
           if strcmp(Trials_Info{1, j}.General.Description, 'Experiment')
                scores = cat(2, scores, Trials_Info{1, j}.General.Score);
                trials = cat(2, trials, j);
           end
        else
            scores = cat(2, scores, Trials_Info{1, j}.General.Score);
            trials = cat(2, trials, j);
        end
    end

    scores = scores(scores ~= 0);
    trials = trials(scores ~= 0);

    
    %% clustering scores with KNN
   
    num_clusters = 3; % Define number of clusters
    
    max_score = max(scores);
    min_score = min(scores);
    initial_centers = [min_score, (min_score + max_score)/2, max_score]';

    if subject == 6
        initial_centers = [1, 2, 7]';
    elseif subject == 7
        initial_centers = [1, 1.5, 5]';
    elseif subject == 8
        initial_centers = [1, 5, 8]';
    elseif subject == 9
        initial_centers = [1, 5, 6]';
    elseif subject == 14
        initial_centers = [1, 3, 8]';
    elseif subject == 16
        initial_centers = [1, 3, 8]';
    elseif subject == 17
        initial_centers = [1, 3, 8]';    
    end


    % Apply k-means clustering per subject
    [idx, C] = kmeans(scores', num_clusters, 'Start', initial_centers);
    
    % Sort cluster centers and get sorted indices
    [~, sorted_indices] = sort(C);
    
    % Assign clusters based on sorted order
    new_idx = zeros(size(idx));
    for j = 1:num_clusters
        new_idx(idx == sorted_indices(j)) = j;
    end
        
    new_idx = new_idx';
    

end