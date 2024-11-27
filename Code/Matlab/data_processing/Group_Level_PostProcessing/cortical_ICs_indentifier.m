% Manually identified cortical ICs for each subject. The output is used for 
% subsequent clustering, focusing exclusively on cortical ICs and excluding
% ICs categorized as "muscle" or "other."

function Output = cortical_ICs_indentifier(subject_list, ALLEEG)

    % 1st column: ICs which were categorized to be Brain by ICLabel and had
    % higher than 50% probability.
    % 2nd column: ICs which need further investigation to check if they are
    % potential Brain ICs or not.
    ICs = cell(length(ALLEEG), 3);
    P_max = cell(length(ALLEEG), 1);

    for i = 1:length(ALLEEG)
        
        P = ALLEEG(i).etc.ic_classification.ICLabel.classifications;
        [maxValues, maxIndices] = max(P, [], 2);
        P_max{i, 1} = [maxValues, maxIndices];

        % Maximum Probability is for Brain and P(Brain) is higher than 50%
        ICs{i, 1} = find(and(P_max{i, 1}(:, 2) == 1, P_max{i, 1}(:, 1) > 0.5));

        % These ICs need further investigation:
        % Maximum Probability is for Brain but P(Brain) is less than 50%
        ICs{i, 2} = find(and(P_max{i, 1}(:, 2) == 1, P_max{i, 1}(:, 1) <= 0.5));
        % Maximum Probability is for "Muscle" but P(Muscle) is less than 50%
        ICs{i, 2} = cat(1, ICs{i, 2}, find(and(P_max{i, 1}(:, 2) == 2, P_max{i, 1}(:, 1) <= 0.5)));
        % Maximum Probability is for "Other" but P(Other) is less than 50%
        ICs{i, 2} = cat(1, ICs{i, 2}, find(and(P_max{i, 1}(:, 2) == 7, P_max{i, 1}(:, 1) <= 0.5)));

    end


    for i = 1:length(ALLEEG)
    
        for j = 1:length(ICs{i, 2})
            [~, ~, ~] = pop_prop_extended(ALLEEG(i), 0, ICs{i, 2}(j, 1), ...
                NaN, {'freqrange', [1 60]});
        end

    end
    
    %% GUI to select accepted potential Brain ICs
    ic_selection_gui(ALLEEG, ICs, subject_list)
    save(fullfile(pwd, 'Brain_PotentialBrain_AcceptedPotentialBrain.mat'), 'ICs')
    







    % output = cell(length(subject_list), 3);
    % 
    % for i = 1:length(subject_list)
    %     output{i, 1} = subject_list{1, i};
    % end
    % 
    % %% subject 5
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{1, 2} = [3, 4, 5, 7, 11, 13, 19, 20, 21, 26, 29, 30, 34, 36, 37, 39, 40, 44, 46];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{1, 3} = [8, 41];
    % 
    % %% subject 6
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{2, 2} = [11, 20, 24, 28, 34, 39, 48, 49];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{2, 3} = 27;
    % 
    % %% subject 7
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{3, 2} = [2, 4, 9, 12, 14, 15, 44, 57];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{3, 3} = [22, 33, 35, 37, 45];
    % 
    % %% subject 8
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{4, 2} = [1, 2, 3, 4, 5, 7, 12, 14, 15, 18, 20, 23, 29, 32, 34, 37, 39, 44, 46, 52];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{4, 3} = [6, 9, 10, 16, 19, 30, 31, 33, 35];
    % 
    % %% subject 9
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{5, 2} = [2, 3, 4, 5, 8, 13, 16, 25, 51];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{5, 3} = [7, 21, 25];
    % 
    % %% subject 10
    % % brain IC by iclabel (after checking manually, some ICs are removed based on PSD content.)
    % output{6, 2} = [1, 3, 4, 5, 6, 7, 8, 10, 11, 13, 21, 31, 33, 35, 39, 49, 50];
    % % non-brain IC by iclabel but it can be brain (after checking manually)
    % output{6, 3} = [12, 26, 28, 36, 38, 40, 44, 46];
    % 

end