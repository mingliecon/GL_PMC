function run_simulation(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_simulation(number of cores)"
    % Example, put "run_simulation(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 25/May/2025

    %% set up number of cores and save current directory

    if nargin > 0
        core = varargin{1};
    else
        core = 4;
    end

    fprintf('\nYou are running with %d cores\n\n', core)
    orig_dir = pwd;
    save('_tempvars.mat', 'core', 'orig_dir');

    %% table 1_2
    disp('Generating table 1 and 2 ...');

    load('_tempvars.mat');

    % Move to subfunction (of table 1 and 2) directory
    cd(fullfile('table_1_2', 'src_table_1_2'));
    
    % Call function for table 3
    evalc('run_table_1_2()');
    
    % Return to original folder
    cd(orig_dir);

    disp('Table 1 and 2  generated');

    clear; close all; clc;

    %% table 7

    load('_tempvars.mat');
    % Move to subfunction (of table 7) directory
    disp('Generating table 7...');

    cd(fullfile('table_7', 'src_table_7'));
    
    % Call function for table 7
    evalc('run_table_7(core)');
    
    % Return to original folder
    cd(orig_dir);
    disp('Table 7 generated');
    
    %% table 8

    load('_tempvars.mat');
    % Move to subfunction (of table 8) directory
    disp('Generating table 8...');

    cd(fullfile('table_8', 'src_table_8'));
    
    % Call function for table 8
    evalc('run_table_8(core)');
    
    % Return to original folder
    cd(orig_dir);
    disp('Table 8 generated');

    % Call function for figure 2&3
    evalc('run_fig_2_3(core)');
    
    % Return to original folder
    cd(orig_dir);
    disp('Figure 2 & 3 generated');
    



end


    