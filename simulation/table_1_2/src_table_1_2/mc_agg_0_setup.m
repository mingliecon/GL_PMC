%% Setup
function [N, D, J, T, B_first, B, B_last, sigma, beta0, core] = mc_agg_0_setup(varargin)
    
    if ~isempty(varargin)
        temp = varargin;
        N = temp{1}*1000;
        D = temp{2};
        J = temp{3};
        T = temp{4};
        B_first = temp{5};
        B = temp{6};
        core = temp{7};
    else
        %default
        N = 10000; % # of individuals in the market
        D = 4;  % # of dimension of observable product characteristics (>= 3)
        J = 3;   % # of alternatives including outside option 0 (>= 3)
        T = 2;  % # of time periods (>= 2)
        B_first = 1;
        B = 1000;
        core = 4;
    end
    
    %evalc('parpool(core)');
    
    logname = sprintf("Table2_N=%d_D=%d_J=%d_T=%d_bigger_index.log",N,D,J,T);
    
    if exist(logname)
        diary off
        delete(logname);
    end
    
    diary(logname)
    
    % --------------------------------------------
    %B_first = 1;
    %B = 100; % # of monte carlo repititions
    B_last = B_first + B - 1;
    sigma = J;
    
    % B_first
    % B
    % B_last
    % N
    % D
    % J
    % T
    % sigma
    
    beta0 = [2; repmat(1,[D-1,1])];

end
%%


%[N, D, J, T, B_first, B, B_last, sigma, beta0, core] = mc_agg_0_setup()