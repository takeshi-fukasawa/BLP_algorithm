%% Check Iaria and Wang (2025)'s function

%% Set parameters
J = 25;          % Number of products
K2 = 3;          % Number of coefficients with random coefficients
NbsSample = 1000; % Number of simulation draws

% --- DGP (Data Generating Process) ---
mean_utility = randn(J, 1); 

X2 = [ones(J,1), randn(J, K2-1)]; 

theta2 = tril(0.5 * randn(K2, K2)) + diag(rand(K2, 1)); 

random_coef = randn(K2, NbsSample); 


true_ms_sim = rand(J, 1);
market_share = true_ms_sim / (sum(true_ms_sim) + 0.1);

profile on
[diff_ms, dev_diff_ms] = GenerateMarketShare_parametric(...
    mean_utility, NbsSample, random_coef, X2, theta2, market_share);
profile off
profsave(profile('info'), 'others./profile_results')

