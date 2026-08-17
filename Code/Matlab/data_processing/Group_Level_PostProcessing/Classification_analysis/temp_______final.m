function bestModel = temp_______final(X, Y)
% runEEGClassificationWithFeatureSelection
% - Performs ANOVA-based feature selection (F > 1.3)
% - Applies SMOTE for class balancing inside CV folds
% - Runs repeated 5x5 cross-validation
% - Trains models with full and selected features, selects the best one

%% 1️⃣ Ensure Proper Data Format
if istable(X), X = X{:,:}; end  % Convert table to numeric matrix
if istable(Y), Y = Y{:,:}; end  % Convert table to vector
Y = categorical(Y);  % Convert labels to categorical

%% 2️⃣ Compute ANOVA F-Scores for Each Feature
numFeatures = size(X, 2);
fScores = zeros(numFeatures, 1);
pValues = zeros(numFeatures, 1);

for i = 1:numFeatures
    featureData = X(:, i);
    [p, tbl] = anova1(featureData, Y, 'off'); % 'off' prevents figure output
    fScores(i) = tbl{2,5}; % Extract F-statistic
    pValues(i) = p; % Store p-value
end

%% 3️⃣ Select Features with F-Score > 1.3
selectedIdx = find(fScores > 1.3);
fprintf('Selected %d out of %d features (F > 1.3)\n', length(selectedIdx), numFeatures);
X_selected = X(:, selectedIdx); % Feature subset

%% 4️⃣ Setup Repeated 5-Fold Cross-Validation
numRepeats = 5;
cv = cvpartition(Y, 'KFold', 5, 'Stratify', true);

acc_full = zeros(numRepeats,1);
acc_selected = zeros(numRepeats,1);
bestModels = cell(numRepeats,2); % Store models (full & selected)

%% 5️⃣ Repeated Cross-Validation Loop
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
            pred_selected = predict(mdl_selected, X_test);
            foldAcc_selected(fold) = sum(pred_selected == Y_test) / length(Y_test);
        else
            foldAcc_selected(fold) = 0; % No selected features, skip this model
        end
    end

    % Store mean accuracy for this repetition
    acc_full(rep) = mean(foldAcc_full);
    acc_selected(rep) = mean(foldAcc_selected);
    bestModels{rep,1} = mdl_full;
    bestModels{rep,2} = mdl_selected;

    fprintf('Repeat %d - Full: %.2f%% | Selected: %.2f%%\n', rep, acc_full(rep) * 100, acc_selected(rep) * 100);
end

%% 6️⃣ Compare Final Models and Select the Best One
meanAcc_full = mean(acc_full);
meanAcc_selected = mean(acc_selected);

fprintf('\nFinal Accuracy (Mean across 5x5 CV):\n');
fprintf('  Full Features: %.2f%%\n', meanAcc_full * 100);
fprintf('  Selected Features: %.2f%%\n', meanAcc_selected * 100);

if meanAcc_selected > meanAcc_full
    [~, bestRep] = max(acc_selected);
    bestModel = bestModels{bestRep,2};
    fprintf('✅ Selected Features Performed Better! Using the feature-reduced model.\n');
else
    [~, bestRep] = max(acc_full);
    bestModel = bestModels{bestRep,1};
    fprintf('✅ Full Features Performed Better! Using the full feature model.\n');
end
end

%% Helper Function: SMOTE Balancing
function [X_smote, Y_smote] = smoteBalancing(X, Y, smoteRatio)
% SMOTE applied within training folds to prevent data leakage

Ycat = categorical(Y);
classLabels = categories(Ycat);
counts = countcats(Ycat);
maxCount = max(counts);

X_smote = X;
Y_smote = Ycat;

for c = 1:numel(classLabels)
    className = classLabels{c};
    idxClass = (Ycat == className);
    classCount = sum(idxClass);

    if classCount < maxCount
        numToGenerate = round((maxCount - classCount) * smoteRatio);
        if numToGenerate <= 0, continue; end

        X_minority = X(idxClass, :);
        nMinor = size(X_minority, 1);

        k = 5; % K-nearest neighbors for interpolation
        idx = knnsearch(X_minority, X_minority, 'K', k+1); 
        synthSamples = zeros(numToGenerate, size(X,2));

        for s = 1:numToGenerate
            baseIndex = randi(nMinor);
            neighborIndex = idx(baseIndex, randi(k) + 1);
            diff = X_minority(neighborIndex,:) - X_minority(baseIndex,:);
            gap = rand();
            synthSamples(s,:) = X_minority(baseIndex,:) + gap .* diff;
        end

        X_smote = [X_smote; synthSamples];
        Y_smote = [Y_smote; repmat(categorical({className}), numToGenerate, 1)];
    end
end
Y_smote = removecats(Y_smote);
end
