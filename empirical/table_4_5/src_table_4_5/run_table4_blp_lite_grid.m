function [beta_hat] = run_table4_blp_lite_grid(varargin)

    if nargin > 0
        M = varargin{1};
        sigma = varargin{2};
    else
        M = 2000;
        sigma = 1;
    end
    %% Create log file
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    logname = sprintf('table4_blp_grid_log_%s.log', timestamp);
    if exist(logname, 'file'), delete(logname); end
    diary(logname);
    diary on;
    fprintf('Log started at %s\n', timestamp);

    %% Load data
    fprintf('Loading data...\n');
    load SSS2018.mat;

    num_weeks = T; num_stores = N; num_prods = J;
    % Reshape data into [products × store-weeks]
    shares = reshape(Y, [num_prods, num_stores * num_weeks]);  % s_jt
    price = reshape(X(:,1), [num_prods, num_stores * num_weeks]);  % p_jt
    deal = reshape(X(:,2), [num_prods, num_stores * num_weeks]);  % deal_jt
    
    % Stack product characteristics (price, deal, interaction)
    X_jt = [price(:), deal(:), price(:) .* deal(:)]; 
    shares = shares(:);  
    
    %% Step 1: Draw random coefficients (η vectors)
    rng(12);    % For reproducibility 
    eta_draws = sigma * randn(3, M);  

    %% Step 2: Create candidate β vectors on the unit ball
    % Create grid 
    grid_values = -1:0.1:1;
    [beta1, beta2, beta3] = ndgrid(grid_values, grid_values, grid_values);
    beta_candidates = [beta1(:), beta2(:), beta3(:)];
    
    % Remove duplicates and normalize to unit ball
    beta_candidates = unique(beta_candidates, 'rows');
    beta_candidates = beta_candidates ./ vecnorm(beta_candidates, 2, 2);
    
    L = size(beta_candidates, 1);  % Number of candidate β vectors
    
    %% Step 3: Evaluate each candidate β
    Q_values = zeros(L, 1);  % Store criterion values
    
    for l = 1:L
        current_beta = beta_candidates(l, :)';
        
        % Compute utilities with random coefficients (δ + X'η)
        % Note: In this simple version, δ = X'β (no ξ or other terms)
        utility = X_jt * current_beta + X_jt * eta_draws;
        
        % Compute market shares using logit formula
        exp_utility = exp(utility);
        sum_exp_utility = sum(exp_utility, 1);  % Sum across products for each draw
        % i delete 1 here? (july 17)
        %pred_shares = mean(exp_utility ./ (1 + sum_exp_utility), 2);  % Average over draws
        pred_shares = mean(exp_utility ./ (sum_exp_utility), 2);  % Average over draws
        
        % Compute criterion function (Euclidean distance squared)
        Q_values(l) = sum((pred_shares - shares).^2);
    end
    
    %% Step 4: Find the optimal β
    [min_Q, min_idx] = min(Q_values);
    beta_hat = beta_candidates(min_idx, :);
    fprintf('Optimal beta found:\n');
    disp(beta_hat);

    filename = sprintf('beta_hat_results_%s.mat');
    
    % Save beta_hat and other relevant variables
    save(filename, 'beta_hat', 'min_Q', 'beta_candidates', 'Q_values');
     diary off;
end