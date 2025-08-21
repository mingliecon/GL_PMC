function run_table_4_5(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_table_4_5(number of cores)"
    % Example, put "run_table_4_5(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 23/May/2025

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    logfilename = sprintf('log_run_table_4_5_%s.txt', timestamp);
    diary(logfilename);
    diary on;
    fprintf('Logging started at %s\n', datestr(now));
    fprintf('Log file: %s\n\n', logfilename);
    if nargin > 0
        core = varargin{1};
    else
        core = 4;
    end
    %clear; close all; clc;
    %rng(123)
    
    %% Step 1: Run both scripts
    fprintf('Starting generation of Table 4...\n');
    tic  % Start timer for Table 4
    run('Table_4(core).m');        
    
    run('Table_4_nlogn(core).m');  
    elapsedTime = toc;  % Stop timer for Table 4
    fprintf('Finished generation of Table 4. Time elapsed: %.2f seconds\n', elapsedTime);

    
    %% Step 2: Load both result files
    load('beta_B.mat');        % Results from Table_4_5.m
    load('beta_B_nlogn.mat');  % Results from Table_4_5_nlogn.m
    
    %% Step 3: Generate table 4
    
    price_beta_m = round(beta_B(2,1),4);
    price_beta_l = round(beta_B(1,1),4);
    price_beta_u = round(beta_B(3,1),4);
    promo_beta_m = round(beta_B(2,2),4);
    promo_beta_l = round(beta_B(1,2),4);
    promo_beta_u = round(beta_B(3,2),4);
    cross_beta_m = round(beta_B(2,3),4);
    cross_beta_l = round(beta_B(1,3),4);
    cross_beta_u = round(beta_B(3,3),4);
    
    price_beta_log_m = round(beta_B_nlogn(2,1),4);
    price_beta_log_l = round(beta_B_nlogn(1,1),4);
    price_beta_log_u = round(beta_B_nlogn(3,1),4);
    promo_beta_log_m = round(beta_B_nlogn(2,2),4);
    promo_beta_log_l = round(beta_B_nlogn(1,2),4);
    promo_beta_log_u = round(beta_B_nlogn(3,2),4);
    cross_beta_log_m = round(beta_B_nlogn(2,3),4);
    cross_beta_log_l = round(beta_B_nlogn(1,3),4);
    cross_beta_log_u = round(beta_B_nlogn(3,3),4);
    
    
    %%
    
    variables = {'Price_ijt', 'Promo_ijt', 'Price * Promo_ijt'};
    % Create a cell array to hold the data
    headers = {'Variable', 'β_e=0^m', 'β^i_e=0', 'β^u_e=0', 'β_e=0.14^m', 'β^i_e=0.14', 'β^u_e=0.14'};
    data = {
        variables{1}, price_beta_m, price_beta_l, price_beta_u, price_beta_log_m, price_beta_log_l, price_beta_log_u;
        variables{2}, promo_beta_m, promo_beta_l, promo_beta_u, promo_beta_log_m, promo_beta_log_l, promo_beta_log_u;
        variables{3}, cross_beta_m, cross_beta_l, cross_beta_u, cross_beta_log_m, cross_beta_log_l, cross_beta_log_u
    };
    
    % Convert to table
    T = cell2table(data, 'VariableNames', headers);
    
    %save table 4

    table4_save1 = '../result_table_4_5/table4.csv';
    writetable(T, table4_save1);
    table4_save2 = '../../result_empirical/table4.csv';
    writetable(T, table4_save2);
    
    % Display confirmation message
    %fprintf('Results saved to: %s\n', fullfile(pwd, csvFileName));
    
    % Optional: Display formatted table in command window (your original code)
    fprintf('\nTable 4: Empirical Application: Estimation Results\n\n');
    fprintf('%-25s %-15s %-25s %-15s %-25s\n', ...
        '', ...
        'β_{e=0}^m', '[β^i, β^u]_{e=0}', ...
        'β_{e=.14}^m', '[β^i, β^u]_{e=.14}');
    fprintf('%-25s %-15.4f [%-7.4f, %-7.4f] %-15.4f [%-7.4f, %-7.4f]\n', ...
        variables{1}, ...
        price_beta_m, price_beta_l, price_beta_u, ...
        price_beta_log_m, price_beta_log_l, price_beta_log_u);
    fprintf('%-25s %-15.4f [%-7.4f, %-7.4f] %-15.4f [%-7.4f, %-7.4f]\n', ...
        variables{2}, ...
        promo_beta_m, promo_beta_l, promo_beta_u, ...
        promo_beta_log_m, promo_beta_log_l, promo_beta_log_u);
    fprintf('%-25s %-15.4f [%-7.4f, %-7.4f] %-15.4f [%-7.4f, %-7.4f]\n', ...
        variables{3}, ...
        cross_beta_m, cross_beta_l, cross_beta_u, ...
        cross_beta_log_m, cross_beta_log_l, cross_beta_log_u);
    
    disp('Finish Generation of table 4.m...')

    %% Step 4: Run BLP method
    run('run_table4_blp_lite_grid.m');
    %%
    load("beta_hat_results_.mat")
    
    %% Step 5: Generate table 5
    run('cyclic_src/run_cyclic.m');
    %%
    run('ols_olsfe_mlogitfe_src/run_ols_olsfe_mlogitfe.m')
    %%
    load('beta_B.mat');        % Results from Table_4_5.m
    load('beta_B_nlogn.mat');  
    
    load('cyclic_result.mat')
    load('ols_olsfe_mlogitfe_result.mat')
    load("beta_hat_results_.mat")
    price_beta_m = round(beta_B(2,1),4);
    price_beta_l = round(beta_B(1,1),4);
    price_beta_u = round(beta_B(3,1),4);
    promo_beta_m = round(beta_B(2,2),4);
    promo_beta_l = round(beta_B(1,2),4);
    promo_beta_u = round(beta_B(3,2),4);
    cross_beta_m = round(beta_B(2,3),4);
    cross_beta_l = round(beta_B(1,3),4);
    cross_beta_u = round(beta_B(3,3),4);
    
    price_beta_log_m = round(beta_B_nlogn(2,1),4);
    price_beta_log_l = round(beta_B_nlogn(1,1),4);
    price_beta_log_u = round(beta_B_nlogn(3,1),4);
    promo_beta_log_m = round(beta_B_nlogn(2,2),4);
    promo_beta_log_l = round(beta_B_nlogn(1,2),4);
    promo_beta_log_u = round(beta_B_nlogn(3,2),4);
    cross_beta_log_m = round(beta_B_nlogn(2,3),4);
    cross_beta_log_l = round(beta_B_nlogn(1,3),4);
    cross_beta_log_u = round(beta_B_nlogn(3,3),4);
    
    
    %col1
    
    col1 = [round(price_beta_log_m, 4); 
            round(promo_beta_log_m, 4); 
            round(cross_beta_log_m, 4)];
    
    col2 = round(beta10, 4);
    col3 = round(ols_norm, 4);
    col4 = round(olsfe_norm, 4);
    col5 = round(mlogitfe_norm, 4);
    col6 = round(beta_hat(:), 4);
    
    
    % Create row names
    rowNames = {'Price_{ijt}'; 'Promo_{ijt}'; 'Price_{ijt} × Promo_{ijt}'};
    
   T = table(col1, col2, col3, col4, col5, col6, ...
    'RowNames', rowNames, ...
    'VariableNames', {'β_c=0.14^n', 'β_CyclicMono', 'β_OLS', 'β_OLS-FE', 'β_MLogit-FE', 'β_RCLM'});

    % Display the table
    disp('Table 5 with BLP estimates:')
    disp(T);

    writetable(T, '../result_table_4_5/table5.csv', 'WriteRowNames', false);
    writetable(T, '../../result_empirical/table5.csv', 'WriteRowNames', false);
    %clear; close all; clc;
    
    fprintf('Logging finished at %s\n', datestr(now));
    diary off;
end
