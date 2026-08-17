function [bestModel, meanAccuracy, stdAccuracy, acc_full, acc_selected, featureTable] = ...
    runClassificationWithSMOTE(X, Y, subject, folder_path)
    
    % runClassificationWithSMOTE function
    % Implements:
    %  - SMOTE inside training folds (to prevent data leakage)
    %  - Repeated 5-fold cross-validation
    %  - Automated model selection with hyperparameter optimization (SVM, Naive Bayes, NN)
    %  - Stratified cross-validation for imbalanced EEG datasets
    
    % Inputs:
    %   X - [nSamples x nFeatures] EEG feature matrix
    %   Y - [nSamples x 1]         Class labels (categorical or numeric)
    %
    % Outputs:
    %   finalModel   - Best performing model across all repeated folds
    %   meanAccuracy - Mean accuracy across repetitions
    %   stdAccuracy  - Standard deviation of accuracy
    

    logFile = fullfile(folder_path, ['subject ', num2str(subject), ' - classification_log.txt']);
    fid = fopen(logFile, 'w'); % Open log file for writing
    

    %% Ensure Proper Data Format
    if istable(X), X = X{:,:}; end  % Convert table to numeric matrix
    if istable(Y), Y = Y{:,:}; end  % Convert table to vector
    Y = categorical(Y);  % Convert labels to categorical
    
    fprintf(fid, "Starting classification with feature selection...\n");
    

    %% Compute ANOVA F-Scores and Eta-Squared for Each Feature
    numFeatures = size(X, 2);
    fScores = zeros(numFeatures, 1);
    pValues = zeros(numFeatures, 1);
    etaSquared = zeros(numFeatures, 1);
    
    for i = 1:numFeatures
        featureData = X(:, i);
        [p, tbl] = anova1(featureData, Y, 'off'); % 'off' prevents figure output
        fScores(i) = tbl{2,5}; % Extract F-statistic
        pValues(i) = p; % Store p-value
        
        % Compute eta-squared
        SS_total = tbl{4,2}; % Sum of squares total
        SS_effect = tbl{2,2}; % Sum of squares for the effect
        etaSquared(i) = SS_effect / SS_total;
    end

    
    %% Rank Features Based on F-Score and Save
    % Create a table with feature indices, F-scores, and p-values
    featureTable = table((1:numFeatures)', fScores, pValues, -log(pValues), etaSquared, ...
        'VariableNames', {'FeatureIndex', 'FScore', 'PValue', '-Log(PValue)', 'Etta_squared'});
    
    % Sort features by F-score in descending order
    featureTable = sortrows(featureTable, 'FScore', 'descend');
    
    % Determine selected features (etta-squared > 0.05)
    featureTable.Selected = featureTable.Etta_squared > 0.05;
    
    % Save feature ranking for later use
    selectedIdx = featureTable.FeatureIndex(featureTable.Selected);
    save([folder_path, '\feature_selection.mat'], 'featureTable');
    
    fprintf(fid, 'Saved selected features to feature_selection.mat\n');
    fprintf(fid, 'Selected %d out of %d features (etta-squared > 0.05)\n', length(selectedIdx), numFeatures);
   
    
    %% Setup Repeated 5-Fold Cross-Validation
    numRepeats = 10;
    cv = cvpartition(Y, 'KFold', 5, 'Stratify', true);
    
    acc_full = zeros(numRepeats,1);
    acc_selected = zeros(numRepeats,1);
    bestModels = cell(numRepeats,2); % Store models (full & selected)
    

    %% Repeated Cross-Validation Loop
    for rep = 1:numRepeats
        if rep > 1
            cv = repartition(cv); % Shuffle folds for next repetition
        end
    
        foldAcc_full = zeros(cv.NumTestSets, 1);
        foldAcc_selected = zeros(cv.NumTestSets, 1);
    
        for fold = 1:cv.NumTestSets
            % Get train/test splits
            trainIdx = training(cv, fold);
            testIdx = test(cv, fold);
            X_train = X(trainIdx, :);
            Y_train = Y(trainIdx);
            X_test = X(testIdx, :);
            Y_test = Y(testIdx);
    
            % Apply SMOTE to training set only
            [X_train_bal, Y_train_bal] = smoteBalancing(X_train, Y_train, 1.0);
            X_train_bal_selected = X_train_bal(:, selectedIdx); % Apply feature selection
    
            % Inner CV for hyperparameter tuning
            opts = struct('CVPartition', cvpartition(Y_train_bal, 'KFold', 5), ...
                          'MaxObjectiveEvaluations', 30, ...
                          'UseParallel', true);
    
            learnerTypes = ["svm", "nb", "net"];
    
            % Train model on FULL features
            mdl_full = fitcauto(X_train_bal, Y_train_bal, ...
                                'Learners', learnerTypes, ...
                                'OptimizeHyperparameters', 'auto', ...
                                'HyperparameterOptimizationOptions', opts);
            pred_full = predict(mdl_full, X_test);
            foldAcc_full(fold) = sum(pred_full == Y_test) / length(Y_test);
    
            % Train model on SELECTED features (only if features were selected)
            if ~isempty(selectedIdx)
                mdl_selected = fitcauto(X_train_bal_selected, Y_train_bal, ...
                                        'Learners', learnerTypes, ...
                                        'OptimizeHyperparameters', 'auto', ...
                                        'HyperparameterOptimizationOptions', opts);
                pred_selected = predict(mdl_selected, X_test(:, selectedIdx));
                foldAcc_selected(fold) = sum(pred_selected == Y_test) / length(Y_test);
            else
                foldAcc_selected(fold) = 0; % No selected features, skip this model
            end

            
        end

        close all; % close optimization progress figures
    
        % Store mean accuracy for this repetition
        acc_full(rep) = mean(foldAcc_full);
        acc_selected(rep) = mean(foldAcc_selected);
        bestModels{rep,1} = mdl_full;
        if ~isempty(selectedIdx)
            bestModels{rep,2} = mdl_selected;
        end
    
        fprintf(fid, 'Repeat %d - Full: %.2f%% | Selected: %.2f%%\n', rep, acc_full(rep) * 100, acc_selected(rep) * 100);
    end
    

    %% Compare Final Models and Select the Best One
    meanAcc_full = mean(acc_full);
    meanAcc_selected = mean(acc_selected);
    stdAcc_full = std(acc_full);
    stdAcc_selected = std(acc_selected);
    
    fprintf(fid, '\nFinal Accuracy (Mean across 10x5 CV):\n');
    fprintf(fid, '  Full Features: %.2f%% (± %.2f)\n', meanAcc_full * 100, stdAcc_full * 100);
    fprintf(fid, '  Selected Features: %.2f%% (± %.2f)\n', meanAcc_selected * 100, stdAcc_selected * 100);
    
    if meanAcc_selected > meanAcc_full
        [~, bestRep] = max(acc_selected);
        bestModel = bestModels{bestRep,2};
        meanAccuracy = meanAcc_selected;
        stdAccuracy  = std(acc_selected);
        fprintf(fid, 'Selected Features Performed Better! Using the feature-reduced model.\n');
    else
        [~, bestRep] = max(acc_full);
        bestModel = bestModels{bestRep,1};
        meanAccuracy = meanAcc_full;
        stdAccuracy  = std(acc_full);
        
        fprintf(fid, 'Full Features Performed Better! Using the full feature model.\n');
    end

    save([folder_path, '\best_model.mat'], 'bestModel'); % Save best model
    
    fclose(fid); % Close log file
    

end


