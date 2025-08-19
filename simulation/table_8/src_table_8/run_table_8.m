function run_table_8(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_table_8(number of cores)"
    % Example, put "run_table_8(12)" in command window and press enter
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
        core = 4;
    end
    fprintf('Running with cores=%d\n',core);

    %% 1. Run all the simulations with different combinations of D, J, T

    MC_AGG_20230909(10,3,3,2,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,3,4,2,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,3,3,4,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,3,4,4,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,4,3,2,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,4,4,2,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,4,3,4,1,1000,core);
    clearvars -except core
    MC_AGG_20230909(10,4,4,4,1,1000,core);
    clearvars -except core

  
 %% 2. Load data for each simulation, for creating table 8

    % Load first file and extract Mid_rsMSE
    data1 = load('results_D3_J3_T2.mat');
    D3_J3_T2_Mid_rsMSE = data1.Mid_rsMSE;
    D3_J3_T2_Mid_MND = data1.Mid_MND;

    % Load second file and extract Mid_rsMSE
    data2 = load('results_D3_J4_T2.mat');
    D3_J4_T2_Mid_rsMSE = data2.Mid_rsMSE;
    D3_J4_T2_Mid_MND = data2.Mid_MND;

    % Load third file and extract Mid_rsMSE
    data3 = load('results_D3_J3_T4.mat');
    D3_J3_T4_Mid_rsMSE = data3.Mid_rsMSE;
    D3_J3_T4_Mid_MND = data3.Mid_MND;

    % Load fourth file and extract Mid_rsMSE
    data4 = load('results_D3_J4_T4.mat');
    D3_J4_T4_Mid_rsMSE = data4.Mid_rsMSE;
    D3_J4_T4_Mid_MND = data4.Mid_MND;

    % Load fifth file and extract Mid_rsMSE
    data5 = load('results_D4_J3_T2.mat');
    D4_J3_T2_Mid_rsMSE = data5.Mid_rsMSE;
    D4_J3_T2_Mid_MND = data5.Mid_MND;

    % Load sixth file and extract Mid_rsMSE
    data6 = load('results_D4_J4_T2.mat');
    D4_J4_T2_Mid_rsMSE = data6.Mid_rsMSE;
    D4_J4_T2_Mid_MND = data6.Mid_MND;

    % Load seventh file and extract Mid_rsMSE
    data7 = load('results_D4_J3_T4.mat');
    D4_J3_T4_Mid_rsMSE = data7.Mid_rsMSE;
    D4_J3_T4_Mid_MND = data7.Mid_MND;

    % Load eighth file and extract Mid_rsMSE
    data8 = load('results_D4_J4_T4.mat');
    D4_J4_T4_Mid_rsMSE = data8.Mid_rsMSE;
    D4_J4_T4_Mid_MND = data8.Mid_MND;

    %% 3. Create table 8
    % Create rMSE table
    rMSE_table = round([
        D3_J3_T2_Mid_rsMSE, D3_J3_T4_Mid_rsMSE, D3_J4_T2_Mid_rsMSE, D3_J4_T4_Mid_rsMSE;
        D4_J3_T2_Mid_rsMSE, D4_J3_T4_Mid_rsMSE, D4_J4_T2_Mid_rsMSE, D4_J4_T4_Mid_rsMSE
    ],4);

    % Create MND table
    MND_table = round([
        D3_J3_T2_Mid_MND, D3_J3_T4_Mid_MND, D3_J4_T2_Mid_MND, D3_J4_T4_Mid_MND;
        D4_J3_T2_Mid_MND, D4_J3_T4_Mid_MND, D4_J4_T2_Mid_MND, D4_J4_T4_Mid_MND
    ],4);



    %% 4. Display formatted table 8 in command window
    fprintf('Table 8: Performance Varying D, J, T\n\n');
    
    % rMSE table
    fprintf('%8s | %6s | %6s | %6s | %6s\n', 'rMSE', 'J=3', 'J=3', 'J=4', 'J=4');
    fprintf('%8s | %6s | %6s | %6s | %6s\n', '', 'T=2', 'T=4', 'T=2', 'T=4');
    fprintf('--------------------------------------------\n');
    fprintf('%8s | %6.4f | %6.4f | %6.4f | %6.4f\n', 'D=3', rMSE_table(1,:));
    fprintf('%8s | %6.4f | %6.4f | %6.4f | %6.4f\n', 'D=4', rMSE_table(2,:));
    fprintf('\n');
    
    % MND table
    fprintf('%8s | %6s | %6s | %6s | %6s\n', 'MND', 'J=3', 'J=3', 'J=4', 'J=4');
    fprintf('%8s | %6s | %6s | %6s | %6s\n', '', 'T=2', 'T=4', 'T=2', 'T=4');
    fprintf('--------------------------------------------\n');
    fprintf('%8s | %6.4f | %6.4f | %6.4f | %6.4f\n', 'D=3', MND_table(1,:));
    fprintf('%8s | %6.4f | %6.4f | %6.4f | %6.4f\n', 'D=4', MND_table(2,:));
    fprintf('\n');

        %% 5. Save table 8
    vertical_combined = [
        rMSE_table(1,:);
        rMSE_table(2,:);
        MND_table(1,:);
        MND_table(2,:)
    ];
    
    % Add row and column labels
    row_names = {'D3_rMSE', 'D4_rMSE', 'D3_MND', 'D4_MND'};
    col_names = {'J3_T2', 'J3_T4', 'J4_T2', 'J4_T4'};
    
    % Convert to table
    combined_table = array2table(vertical_combined, ...
        'RowNames', row_names, ...
        'VariableNames', col_names);
    
    % Ensure directory exists
    if ~exist('../result_table_8', 'dir')
        mkdir('../result_table_8');
    end
    
    % Save to CSV
    writetable(combined_table, '../result_table_8/table8.csv', 'WriteRowNames', true);
    writetable(combined_table, '../../result_simulation/table8.csv', 'WriteRowNames', true);
    
    fprintf('Completed with last parameter = %d\n', core);
end