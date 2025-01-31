function chance = random_chance_binomial(n, p, alpha)

% % Parameters
% n         % Total number of trials
% p         % Probability of success under the null hypothesis
% alpha     % Significance level

% Initialize k to the expected number of successes
k = ceil(n * p);

% Calculate cumulative probability P(X >= k)
P = 1 - binocdf(k - 1, n, p);

% Increase k until P <= alpha
while P > alpha && k <= n
    k = k + 1;
    P = 1 - binocdf(k - 1, n, p);
end

% Check if k exceeds n
if k > n
    error('No value of k satisfies the condition for the given n and p.');
end

chance = (k/n)*100;

% Output the result
fprintf('Minimum number of successes k: %d\n', k);
fprintf('Cumulative probability P(X >= k): %.6f\n', P);


end