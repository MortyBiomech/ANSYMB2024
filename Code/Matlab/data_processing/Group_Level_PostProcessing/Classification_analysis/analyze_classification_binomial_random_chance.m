function analyze_classification_binomial_random_chance(N1, N2, N3, observed_accuracy)
    % Total number of samples
    N_total = N1 + N2 + N3;
    
    % Calculate random chance accuracy (majority class probability)
    p_random = max([N1, N2, N3]) / N_total;
    
    % Convert observed accuracy to number of correct classifications
    observed_correct = round(observed_accuracy * N_total);
    
    % Compute p-value using the binomial test
    p_value = 1 - binocdf(observed_correct - 1, N_total, p_random);
    
    % Determine significance level
    if p_value < 0.05
        significance = 'Significant (p < 0.05)';
    else
        significance = 'Not Significant (p >= 0.05)';
    end
    
    % Display results
    fprintf('Random Chance Accuracy: %.2f%%\n', p_random * 100);
    fprintf('Observed Accuracy: %.2f%%\n', observed_accuracy * 100);
    fprintf('Binomial Test p-value: %.5f\n', p_value);
    fprintf('Result: %s\n', significance);
end