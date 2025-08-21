function run_table_6(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_table_6(number of cores)"
    % Example, put "run_table_6(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 17/May/2025


    if nargin > 0
        core = varargin{1};
    else
        core = 12;
    end
    
    fprintf('\nYou are running with %d cores\n\n', core)
    %clear; close all; clc;
    %rng(123)

    %% 1. Run CyclicMono and load its data files
    disp('Running CyclicMono analysis...');
    run('CyclicMono/run_cyclicmono(core).m');
    %%
    % Load CyclicMono data files
    
    cyclic_0_3 = load('CyclicMono/cyclic_alpha_0_3_beta_percent.mat');
    cyclic_0_5 = load('CyclicMono/cyclic_alpha_0_5_beta_percent.mat');
    cyclic_0_15 = load('CyclicMono/cyclic_alpha_0_15_beta_percent.mat');
    
    disp('CyclicMono data loaded successfully.');
    
    %% 2. Run GL and load its data files
    disp('Running GL analysis...');
    run('GL/run_GL(core).m');
    %%
    % Load GL data files
    
    gl_3 = load('GL/GL_0_3_beta_percent.mat');
    gl_0_5 = load('GL/GL_0_5_beta_percent.mat');
    gl_0_15 = load('GL/GL_0_15_beta_percent.mat');
    
    disp('GL data loaded successfully.');
    
    %% 3. Run OLS/Logit and load its data files
    disp('Running OLS/Logit analysis...');
    run('OLS_LOGIT/run_ols_logit(core).m');
    %%
    % Load OLS/Logit data files
    
    ols_logit_0_3 = load('OLS_LOGIT/alpha_0_3_results.mat');
    ols_logit_0_5 = load('OLS_LOGIT/alpha_0_5_results.mat');
    ols_logit_0_15 = load('OLS_LOGIT/alpha_0_15_results.mat');
    
    disp('OLS/Logit data loaded successfully.');

    %% 4. Run BLP and load its data files
    run('BLP/run_blp_lite_table6.m');
    %%
    load('BLP/sign_recovery_results.mat');
    %%
    beta_blp = sign_recovery_table.Percentage_Correct_Signs;
    %%

    %% 5. Create Table 6
    
    alpha_values = [0.15; 0.3; 0.5];
    
    beta_m = [gl_0_15.b_percent; gl_3.b_percent; gl_0_5.b_percent] * 100;
    
    beta_cyclic_mono = [cyclic_0_15.b_percent; cyclic_0_3.b_percent; cyclic_0_5.b_percent] * 100;
    
    beta_ols = [ols_logit_0_15.ols_pct*100; ols_logit_0_3.ols_pct*100; ols_logit_0_5.ols_pct*100];
    beta_ols_fe = [ols_logit_0_15.ofe_pct*100; ols_logit_0_3.ofe_pct*100; ols_logit_0_5.ofe_pct*100];
    beta_mlogit_fe = [ols_logit_0_15.mfe_pct*100; ols_logit_0_3.mfe_pct*100; ols_logit_0_5.mfe_pct*100];
    
    % Create Table 6 with the new betaBLP column
    Table6 = table(alpha_values, beta_m, beta_cyclic_mono, beta_ols, beta_ols_fe, beta_mlogit_fe, beta_blp, ...
        'VariableNames', {'alpha', 'betam', 'betaCyclicMono', 'betaOLS', 'betaOLS_FE', 'betaMLogit_FE', 'betaRCLM'});
    
    % Display the table
    disp('Table 6: Percentage of Correct Signs of Estimated Coefficients');
    disp(Table6)
        
    %% Save as table6.csv
    writetable(Table6, '../result_table_6/table6.csv');
    writetable(Table6, '../../result_empirical/table6.csv');
    %clear; close all; clc;
end