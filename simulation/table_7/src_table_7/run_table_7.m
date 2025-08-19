function run_table_7(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_table_7(number of cores)"
    % Example, put "run_table_7(12)" in command window and press enter
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

    %% 1. Run PointNID
    disp('Running PointNID...');
    run('PointNID/run_pointnid(core).m');
    %%
    % Load PointNID data
    
    NID_0_01 = load('PointNID/nlogn_0_01_result.mat');
    NID_0_1 = load('PointNID/nlogn_0_1_result.mat');
    NID_1 = load('PointNID/nlogn_1_result.mat');
    
    disp('NID data loaded successfully.');
    
    %% 2. Run PointID
    disp('Running PointID...');
    run('PointID/run_pointid(core).m');
    %% 
    % Load PointNID data
    
    ID_0_01 = load('PointID/idnlogn_0_01_result.mat');
    
    disp('ID data loaded successfully.');
    
    %% 3. Create Table 7
    id_0_01_Ub_rsMSE = ID_0_01.Ub_rsMSE;
    id_0_01_Ub_MND =  ID_0_01.Ub_MND;
    id_0_01_Lb_rsMSE = ID_0_01.Lb_rsMSE;
    id_0_01_Lb_MND = ID_0_01.Lb_MND;
    id_0_01_Mid_rsMSE = ID_0_01.Mid_rsMSE;
    id_0_01_Mid_MND = ID_0_01.Mid_MND;
    %%
    nid_0_01_Ub_rsMSE = NID_0_01.Ub_rsMSE;
    nid_0_01_Ub_MND =  NID_0_01.Ub_MND;
    nid_0_01_Lb_rsMSE = NID_0_01.Lb_rsMSE;
    nid_0_01_Lb_MND = NID_0_01.Lb_MND;
    nid_0_01_Mid_rsMSE = NID_0_01.Mid_rsMSE;
    nid_0_01_Mid_MND = NID_0_01.Mid_MND;
    %%
    nid_0_1_Ub_rsMSE = NID_0_1.Ub_rsMSE;
    nid_0_1_Ub_MND =  NID_0_1.Ub_MND;
    nid_0_1_Lb_rsMSE = NID_0_1.Lb_rsMSE;
    nid_0_1_Lb_MND = NID_0_1.Lb_MND;
    nid_0_1_Mid_rsMSE = NID_0_1.Mid_rsMSE;
    nid_0_1_Mid_MND = NID_0_1.Mid_MND;
    %%
    nid_1_Ub_rsMSE = NID_1.Ub_rsMSE;
    nid_1_Ub_MND =  NID_1.Ub_MND;
    nid_1_Lb_rsMSE = NID_1.Lb_rsMSE;
    nid_1_Lb_MND = NID_1.Lb_MND;
    nid_1_Mid_rsMSE = NID_1.Mid_rsMSE;
    nid_1_Mid_MND = NID_1.Mid_MND;
    
    %%
    % Create the table structure matching the screenshot
    point_ID = {'(i) yes'; '(ii) no'; ''; ''};
    c_hat = {'-'; '0.01'; '0.1'; '1'};
    
    % Extract data from your variables
    % Row 1: With Point ID (using ID_0_01 variables)
    row1_rMSE = round([id_0_01_Mid_rsMSE, id_0_01_Ub_rsMSE, id_0_01_Lb_rsMSE],4);
    row1_MND = round([id_0_01_Mid_MND, id_0_01_Ub_MND, id_0_01_Lb_MND],4);
    
    % Row 2: Without Point ID, c=0.01 (using NID_0_01 variables)
    row2_rMSE = round([nid_0_01_Mid_rsMSE, nid_0_01_Ub_rsMSE, nid_0_01_Lb_rsMSE],4);
    row2_MND = round([nid_0_01_Mid_MND, nid_0_01_Ub_MND, nid_0_01_Lb_MND],4);
    
    % Row 3: Without Point ID, c=0.1 (using NID_0_1 variables)
    row3_rMSE = round([nid_0_1_Mid_rsMSE, nid_0_1_Ub_rsMSE, nid_0_1_Lb_rsMSE],4);
    row3_MND = round([nid_0_1_Mid_MND, nid_0_1_Ub_MND, nid_0_1_Lb_MND],4);
    
    % Row 4: Without Point ID, c=1 (using NID_1 variables)
    row4_rMSE = round([nid_1_Mid_rsMSE, nid_1_Ub_rsMSE, nid_1_Lb_rsMSE],4);
    row4_MND = round([nid_1_Mid_MND, nid_1_Ub_MND, nid_1_Lb_MND],4);
    
    % Combine all data
    all_rMSE = round([row1_rMSE; row2_rMSE; row3_rMSE; row4_rMSE],4);
    all_MND = round([row1_MND; row2_MND; row3_MND; row4_MND],4);
    
    % Create the table
    Table7 = table(point_ID, c_hat, ...
        all_rMSE(:,1), all_rMSE(:,2), all_rMSE(:,3), ...
        all_MND(:,1), all_MND(:,2), all_MND(:,3), ...
        'VariableNames', {...
            'point_ID', 'c_hat', ...
            'beta_m_rMSE', 'beta_u_rMSE', 'beta_l_rMSE', ...
            'beta_m_MND', 'beta_u_MND', 'beta_l_MND' ...
        });
    
    % Format the display to match the screenshot (4 decimal places)
    format short
    disp('Table 7: Performance with and without Point ID: Further Examination');
    disp(Table7);
    
    
    %% Save the table as table6.csv
    writetable(Table7, '../result_table_7/table7.csv');
    writetable(Table7, '../../result_simulation/table7.csv');
    % Clear workspace and close figures
    %clear; close all; clc;
end