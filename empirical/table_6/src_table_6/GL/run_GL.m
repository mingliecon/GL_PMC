function run_GL(varargin)
    
    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "run_GL(number of cores)"
    % Example, put "run_GL(12)" in command window and press enter
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

    fprintf('Running GL with cores=%d\n', core);
    
    %rng(123)
    Gl_Table6(.15, core)
    Gl_Table6(.3, core)
    Gl_Table6(.5, core) 
end