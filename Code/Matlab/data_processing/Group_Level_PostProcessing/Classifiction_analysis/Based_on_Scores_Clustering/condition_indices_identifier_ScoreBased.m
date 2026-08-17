function condition_indices = condition_indices_identifier_ScoreBased(Trials_Info, subject)

    [cluster_idx, trials] = Score_Clustering(Trials_Info, subject);

    S1 = [];
    S2 = [];
    S3 = [];

    for i = 1:length(trials)
        S = cluster_idx(i);
        switch S
            case 1
                S1 = cat(2, S1, trials(i));
            case 2
                S2 = cat(2, S2, trials(i));
            case 3
                S3 = cat(2, S3, trials(i));
        end
    end

    condition_indices.S1 = S1;
    condition_indices.S2 = S2;
    condition_indices.S3 = S3;

end