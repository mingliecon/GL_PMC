function run_empirical(varargin)
    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_emperical(number of cores)"
    % Example, put "run_emperical(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 5/June/2025

    %% set up number of cores and save current directory

    if nargin > 0
        core = varargin{1};
    else
        core = 4;
    end
    fprintf('\nYou are running with %d cores\n\n', core)
    orig_dir = pwd;
    save('_tempvars.mat', 'core', 'orig_dir');

    %% table 3
    disp('Generating table 3 ...');

    load('_tempvars.mat');

    % Move to subfunction (of table 3) directory
    cd(fullfile('table_3', 'src_table_3'));
    
    % Call function for table 3
    evalc('run_table_3()');
    
    % Return to original folder
    cd(orig_dir);

    disp('Table 3 generated');

    clear; close all; clc;


    %% table 4 _5
    
    load('_tempvars.mat');
    % Move to subfunction (of table 4 and 5) directory
    cd(fullfile('table_4_5', 'src_table_4_5'));
    
    % Call function for table 4 and 5
    run_table_4_5(core);
    
    % Return to original folder
    cd(orig_dir);
    clear; close all; clc;

    %% table 6
    disp('Generating table 6...');
    load('_tempvars.mat');
    % Move to subfunction (of table 4 and 5) directory

    cd(fullfile('table_6', 'src_table_6'));
    
    % Call function for table 4 and 5
    evalc('run_table_6(core)');
    
    % Return to original folder
    cd(orig_dir);
    disp('Table 6 Generated...');
    clear; close all; clc;

end


    