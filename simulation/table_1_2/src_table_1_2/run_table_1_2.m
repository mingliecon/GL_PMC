function run_table_1_2(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_table_1_2(number of cores)"
    % Example, put "run_table_1_2(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 05/June/2025

    if nargin > 0
        core = varargin{1};
    else
        core = 4;
    end
    fprintf('Running table_1_2 with cores=%d\n',core);
    %clear; close all; clc;


    %% Run Simulation
    main_mc_agg(1,3,3,2,1,1000,core);
    main_mc_agg(4,3,3,2,1,1000,core);
    main_mc_agg(10,3,3,2,1,1000,core);
    
    %% Load Result Data
    s1 = load('output_N1_D3_J3_T2_metrics_struct.mat');
    s2 = load('output_N4_D3_J3_T2_metrics_struct.mat');
    s3 = load('output_N10_D3_J3_T2_metrics_struct.mat');
    
      %% Create table 1
    
    m3 = s3.metrics_struct;
    table1_values = [
        round(m3.Mid_bias, 4);
        round(m3.Ub_MD, 4);
        round(m3.Lb_MD, 4);
        round(m3.dUL_mean, 4);
        round(m3.SD, 4);
        round(m3.RMSE, 4);
        round(repmat(m3.Mid_rsMSE, 1, 3), 4);  
        round(repmat(m3.Mid_MND, 1, 3), 4)     
    ];
    
    % Define row and column names
    table1_rowNames = {
        'mid bias';
        'upper bias';
        'lower bias';
        'mean(u-l)';
        'standard deviation';
        'root MSE (coord)';
        'root MSE (vector)';
        'mean norm deviation (MND)';
    };
    
    table1_colNames = {'beta_1', 'beta_2', 'beta_3'};  % for D = 3
    
    % Create and display the table
    T1 = array2table(table1_values, 'VariableNames', table1_colNames, 'RowNames', table1_rowNames);
    
    disp(' ')
    disp('Table 1')
    disp('------------------------------------------------')
    
    disp(T1)
    
    % Export to CSV
    writetable(T1, '../result_table_1_2/table1.csv', 'WriteRowNames', true);
    writetable(T1, '../../result_simulation/table1.csv', 'WriteRowNames', true);


    %% Create table 2
    m1 = s1.metrics_struct;
    m2 = s2.metrics_struct;
    % Calculate summary statistics
    bias_sum = [
        sum(abs(m3.Mid_bias)), 
        sum(abs(m2.Mid_bias)), 
        sum(abs(m1.Mid_bias))
    ];
    
    dUL_sum = [
        sum(m3.dUL_mean),  
        sum(m2.dUL_mean), 
        sum(m1.dUL_mean)
    ];
    
    rMSE = [
        m3.Mid_rsMSE,  
        m2.Mid_rsMSE, 
        m1.Mid_rsMSE
    ];
    
    MND = [
        m3.Mid_MND, 
        m2.Mid_MND, 
        m1.Mid_MND
    ];
    
    
    % 3. Assemble Top Table
    N_values = [10000; 4000; 1000];
    topData = round([bias_sum, dUL_sum, rMSE, MND], 4); % Combine column vectors
    topRowNames = {'N = 10,000'; 'N = 4,000'; 'N = 1,000'};
    topColNames = {'SumBias', 'SumMeanUL', 'rMSE', 'MND'};
    
    % Create table with row names
    topTable = array2table(topData, 'VariableNames', topColNames, 'RowNames', topRowNames);
    
    disp(' ')
    disp('Upper Part of Table 2')
    disp('------------------------------------------------')
    
    disp(topTable)
    writetable(topTable, '../result_table_1_2/table2_upper_part.csv', 'WriteRowNames', true);
    writetable(topTable, '../../result_simulation/table2_upper_part.csv', 'WriteRowNames', true);   

    N_values = [10000; 4000; 1000];
    
    % 4. Scaling & Ratios - Lower Part of Table 2
    % First two columns 
    scaling_sqrt = sqrt(N_values(1:2)/1000); 
    scaling_cbrt = (N_values(1:2)/1000).^(1/3);   
    
    % Last two columns 
    rMSE_ratio = rMSE(3) ./ rMSE(1:2);  
    MND_ratio = MND(3) ./ MND(1:2);    
    
    % Create the bottom table data with rounding
    bottomData = [
        round(scaling_sqrt, 2), round(scaling_cbrt, 2), ...
        round(rMSE_ratio, 2), round(MND_ratio, 2)
    ];
    
    % Column names matching your screenshot
    bottomColNames = {
        '(N/1,000)^1/2', ...
        '(N/1,000)^1/3', ...
        'rMSE_1000/rMSE_N', ...
        'MND_1000/MND_N'
    };
    
    % Row names
    bottomRowNames = {'N = 10,000', 'N = 4,000'};
    
    % Create the table
    bottomTable = array2table(bottomData, ...
        'VariableNames', bottomColNames, ...
        'RowNames', bottomRowNames);

    disp(' ')
    disp('Lower Part of Table 2')
    disp('------------------------------------------------')
    disp(bottomTable)
    
    % Save the bottom table to a CSV file
    writetable(bottomTable, '../result_table_1_2/table2_bottom_part.csv', 'WriteRowNames', true);
    writetable(bottomTable, '../../result_simulation/table2_bottom_part.csv', 'WriteRowNames', true);
    
end
