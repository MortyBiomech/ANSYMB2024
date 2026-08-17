function [X_smote, Y_smote] = smoteBalancing(X, Y, smoteRatio)

    % SMOTE Function - Generates Synthetic Minority Data Within Each Training Fold
    % SMOTEBALANCING Synthetic Minority Oversampling Technique applied inside training folds.
    
    % Convert labels to categorical
    Ycat = categorical(Y);
    
    % Identify class distributions
    classLabels = categories(Ycat);
    counts = countcats(Ycat);
    
    % Find max class count
    maxCount = max(counts);
    
    X_smote = X;
    Y_smote = Ycat;
    
    for c = 1:numel(classLabels)
        className = classLabels{c};
        idxClass = (Ycat == className);
        classCount = sum(idxClass);
    
        % Generate synthetic samples only if class is underrepresented
        if classCount < maxCount
            numToGenerate = round((maxCount - classCount) * smoteRatio);
            if numToGenerate <= 0
                continue;
            end
    
            % Extract minority class samples
            X_minority = X(idxClass, :);
            nMinor = size(X_minority, 1);
    
            % Use k-nearest neighbors to interpolate new samples
            k = 5; % Number of neighbors
            idx = knnsearch(X_minority, X_minority, 'K', k+1); 
            synthSamples = zeros(numToGenerate, size(X,2));
    
            for s = 1:numToGenerate
                baseIndex = randi(nMinor);
                neighborIndex = idx(baseIndex, randi(k) + 1);
                diff = X_minority(neighborIndex,:) - X_minority(baseIndex,:);
                gap = rand();
                synthSamples(s,:) = X_minority(baseIndex,:) + gap .* diff;
            end
    
            % Append synthetic samples
            X_smote = [X_smote; synthSamples];
            Y_smote = [Y_smote; repmat(categorical({className}), numToGenerate, 1)];
        end
    end
    
    % Remove empty categories in case no new data was added
    Y_smote = removecats(Y_smote);

end