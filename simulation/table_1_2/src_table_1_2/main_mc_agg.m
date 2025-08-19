function main_mc_agg(N, D, J, T, B_first, B, core)

    % Create filename
    filename = sprintf('output_N%d_D%d_J%d_T%d', N, D, J, T);
    
    %--- Step 1: Setup ---
    [N, D, J, T, B_first, B, B_last, sigma, beta0, core] = mc_agg_0_setup(N, D, J, T, B_first, B, core);

    %--- Step 2: Simulate data ---
    [X_B, y_B, EyA_B, theta0] = mc_agg_1_dgp(N, D, J, T, B, beta0, sigma);

    clearvars -except N D J T B B_first B_last core beta0 theta0 X_B y_B EyA_B filename;
    %save(sprintf('%s_SimData.mat', filename));

    %--- Step 3: Compute gamma moments ---
    [dX_B, dEy_B] = mc_agg_2_gamma(X_B, y_B, EyA_B, N, D, J, T, B, B_first, B_last);

    clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B filename;

    %--- Step 4: Estimate beta ---

    [beta_B, beta0] = mc_agg_3_beta_D3(N, D, J, T, B, B_first, B_last, core, beta0, theta0, dX_B, dEy_B);


    clearvars -except N D J T B beta0 theta0 beta_B filename;
    %save(sprintf('%s_Results.mat', filename));

    %--- Step 5: Evaluate results ---
    [Mid_bias, Ub_MD, Lb_MD, dUL_mean, Mid_MSE, b_mean, SD, RMSE, Mid_rsMSE, Mid_MND, Mid_SAB, dul_sum] = ...
        mc_agg_4_evaluation(beta_B, beta0);

    % Save as struct
    metrics_struct = struct( ...
        'Mid_bias', Mid_bias, ...
        'Ub_MD', Ub_MD, ...
        'Lb_MD', Lb_MD, ...
        'dUL_mean', dUL_mean, ...
        'Mid_MSE', Mid_MSE, ...
        'b_mean', b_mean, ...
        'SD', SD, ...
        'RMSE', RMSE, ...
        'Mid_rsMSE', Mid_rsMSE, ...
        'Mid_MND', Mid_MND, ...
        'Mid_SAB', Mid_SAB, ...
        'dul_sum', dul_sum ...
    );
    save(sprintf('%s_metrics_struct.mat', filename), 'metrics_struct');

    % Stop logging output
    diary off;

    % Clean up parallel pool (if used)
    %poolobj = gcp('nocreate');
    %evalc('delete(poolobj)');
    
end
