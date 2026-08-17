function analyze_multiple_subjects_classification_accuracies(N_matrix, accuracy_matrix)
    % Number of subjects and classifiers
    num_subjects = size(accuracy_matrix, 1);  % Should be 14
    num_classifiers = size(accuracy_matrix, 2);  % Should be 2

    % Total number of tests (subjects * classifiers)
    m = num_subjects * num_classifiers;

    % Initialize storage for p-values and random chance accuracies
    p_values = zeros(num_subjects, num_classifiers);
    p_random = zeros(num_subjects, 1);

    % Loop through subjects
    for i = 1:num_subjects
        % Extract subject-specific class distributions
        N1 = N_matrix(i, 1);
        N2 = N_matrix(i, 2);
        N3 = N_matrix(i, 3);
        N_total = N1 + N2 + N3;

        % Compute random chance accuracy for this subject
        p_random(i) = max([N1, N2, N3]) / N_total;

        % Convert accuracy values to number of correct classifications
        observed_correct = round(accuracy_matrix(i, :) * N_total);

        % Compute p-values using the binomial test
        for j = 1:num_classifiers
            p_values(i, j) = 1 - binocdf(observed_correct(j) - 1, N_total, p_random(i));
        end
    end

    % Flatten p-values for correction
    p_values_flat = p_values(:);

    % Apply Bonferroni Correction
    alpha = 0.05;
    alpha_bonferroni = alpha / m;
    significant_bonferroni = p_values_flat < alpha_bonferroni;

    % Apply Holm-Bonferroni Correction
    [p_sorted, sort_idx] = sort(p_values_flat);
    thresholds_holm = alpha ./ (m - (1:m) + 1);
    significant_holm = false(size(p_values_flat));
    for i = 1:m
        if p_sorted(i) < thresholds_holm(i)
            significant_holm(sort_idx(i)) = true;
        else
            break;
        end
    end

    % Apply False Discovery Rate (FDR, Benjamini-Hochberg)
    thresholds_fdr = (1:m) / m * alpha;
    significant_fdr = false(size(p_values_flat));
    last_significant = find(p_sorted <= thresholds_fdr, 1, 'last');
    if ~isempty(last_significant)
        significant_fdr(sort_idx(1:last_significant)) = true;
    end

    % Reshape results back to subject x classifier format
    significant_bonferroni = reshape(significant_bonferroni, num_subjects, num_classifiers);
    significant_holm = reshape(significant_holm, num_subjects, num_classifiers);
    significant_fdr = reshape(significant_fdr, num_subjects, num_classifiers);

    % Display results
    for i = 1:num_subjects
        fprintf('Subject %d:\n', i);
        fprintf('  Random Chance Accuracy: %.2f%%\n', p_random(i) * 100);
        for j = 1:num_classifiers
            fprintf('  Classifier %d:\n', j);
            fprintf('    Observed Accuracy: %.2f%%\n', accuracy_matrix(i, j) * 100);
            fprintf('    p-value: %.5f\n', p_values(i, j));
            if significant_bonferroni(i, j)
                fprintf('    Bonferroni Significant: Yes\n');
            else
                fprintf('    Bonferroni Significant: No\n');
            end
            
            if significant_holm(i, j)
                fprintf('    Holm Significant: Yes\n');
            else
                fprintf('    Holm Significant: No\n');
            end
            
            if significant_fdr(i, j)
                fprintf('    FDR Significant: Yes\n\n');
            else
                fprintf('    FDR Significant: No\n\n');
            end
        end
    end
end
