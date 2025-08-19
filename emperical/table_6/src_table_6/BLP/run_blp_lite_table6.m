function sign_recovery_table = run_blp_lite_table6_0706()

    %% Setup Logging
    logname = sprintf('log_run_blp_lite_table6_%s.log', datestr(now,'yyyymmdd_HHMMSS'));
    if exist(logname, 'file'), delete(logname); end
    diary(logname); diary on;
    fprintf('Logging started: %s\n', datestr(now));
    fprintf('Log file: %s\n\n', logname);

    %% Parameters
    alphas = [0.15, 0.3, 0.5];
    M = 1000;          % Number of random coefficient draws
    sigma = 1;         % Scale of random coefficients
    B = 1000;           % Number of simulations

    %% Simulation parameters
    N = 205;           % # of stores/markets
    D = 3;             % # of product characteristics
    J = 4;             % # of products (including outside option)
    T = 10;            % # of time periods

    %% True parameters (known in simulation)
    true_beta = 2 * [-2; 1; 1];  % True beta: [-4; 2; 2]
    true_signs = sign(true_beta)';

    %% Initialize results storage
    results = struct();
    for a = 1:length(alphas)
        alpha = alphas(a);
        correct_sign_counts = zeros(B, 1);

        fprintf('\n========== Alpha = %.2f ==========\n', alpha);

        parfor b = 1:B
            %% Generate simulated data
            rng(b);  % Unique seed for reproducibility

            Z = rand(N,J);
            X = nan(N,D,J,T);

            % Product characteristics
            X(:, 1, :, :) = 4 * rand(N,J,T);  
            X(:, 2, :, :) = (1-alpha)*rand(N,J,T) + alpha*repmat(Z,1,1,T);
            X(:, 3, :, :) = X(:, 1, :, :) .* X(:, 2, :, :);

            A_scale = Z + 1;
            A_location = zeros(N, J);

            Xbeta0 = permute(sum(X .* repmat(true_beta', [N,1,J,T]), 2),[1,3,4,2]);
            A1 = repmat(A_scale, 1, 1, T);
            A2 = repmat(A_location, 1, 1, T);
            AXbeta0 = A1 .* (Xbeta0 + A2);

            expind = exp(AXbeta0);
            sum_expind = sum(expind, 2);
            shares = expind ./ sum_expind;

            %% Reshape for estimation
            shares = permute(shares, [2 1 3]);
            shares = reshape(shares, [J, N*T]);
            shares = shares(1:end-1, :);

            price = permute(X(:,1,:,:), [3 1 4 2]);
            price = reshape(price, [J, N*T]);
            price = price(1:end-1, :);

            deal = permute(X(:,2,:,:), [3 1 4 2]);
            deal = reshape(deal, [J, N*T]);
            deal = deal(1:end-1, :);

            X_jt = [price(:), deal(:), price(:).*deal(:)];
            shares = shares(:);

            %% Estimation
            grid_values = -1:0.1:1;
            [beta1, beta2, beta3] = ndgrid(grid_values, grid_values, grid_values);
            beta_candidates = [beta1(:), beta2(:), beta3(:)];
            beta_candidates = unique(beta_candidates, 'rows');
            beta_candidates = beta_candidates ./ vecnorm(beta_candidates, 2, 2);

            Q_values = zeros(size(beta_candidates, 1), 1);

            for l = 1:size(beta_candidates, 1)
                current_beta = beta_candidates(l, :)';
                utility = X_jt * current_beta + X_jt * (sigma * randn(3, M));
                exp_utility = exp(utility);
                sum_exp_utility = sum(exp_utility, 1);
                % change: i delete 1 here
                %pred_shares = mean(exp_utility ./ (1 + sum_exp_utility), 2);
                pred_shares = mean(exp_utility ./ (sum_exp_utility), 2);
                Q_values(l) = sum((pred_shares - shares).^2);
            end

            [~, min_idx] = min(Q_values);
            beta_hat = beta_candidates(min_idx, :);

            correct_sign_counts(b) = all(sign(beta_hat) == true_signs);

            if mod(b, 10) == 0
                fprintf('Alpha %.2f: Completed simulation %d/%d\n', alpha, b, B);
            end
        end

        results(a).alpha = alpha;
        results(a).sign_recovery_rate = mean(correct_sign_counts) * 100;
        results(a).correct_sign_counts = correct_sign_counts;
    end

    %% Display summary table
    fprintf('\n================ Sign Recovery Rates ==================\n');
    fprintf('Alpha\t\t%% Correct Signs\n');
    fprintf('------------------------------------------\n');
    for a = 1:length(alphas)
        fprintf('%.2f\t\t%.2f%%\n', results(a).alpha, results(a).sign_recovery_rate);
    end
    fprintf('=======================================================\n');

    %% Create output table
    sign_recovery_table = table(...
        alphas', ...
        [results.sign_recovery_rate]', ...
        'VariableNames', {'Alpha', 'Percentage_Correct_Signs'});

    %% Save results
    save('sign_recovery_results.mat', 'results', 'true_beta', 'sign_recovery_table');
    fprintf('Results saved to sign_recovery_results.mat\n');

    %% Close Logging
    fprintf('\nLogging finished: %s\n', datestr(now));
    diary off;
end
