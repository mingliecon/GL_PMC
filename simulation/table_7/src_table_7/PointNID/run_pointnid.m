   
function run_pointnid(varargin)
    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_pointnid(number of cores)"
    % Example, put "run_pointnid(12)" in command window and press enter
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
    


    simulation_nid(10,3,3,2,1,1000,1,core);
    %%
    simulation_nid(10,3,3,2,1,1000,.1,core);
    simulation_nid(10,3,3,2,1,1000,.01,core);
    %MC_AGG_20230909_nlogn(10,3,3,2,1,1000,0,core);
end