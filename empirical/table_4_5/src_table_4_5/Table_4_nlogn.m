function Table_4_nlogn(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "Table_4_nlogn(number of cores)"
    % Example, put "Table_4_nlogn(12)" in command window and press enter
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting core number higher than physical cores will cause error
    %
    % Last Change Date: 23/May/2025

    if nargin > 0
        core = varargin{1};
    else
        core = 4;
    end
    %clear; close all; clc;
   rngSeed = 419;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream



    %% Initialize parallel pool
    if isempty(gcp('nocreate'))
        parpool(core); 
    end


    logname = 'table45.txt';
    if exist(logname)
        diary off
        delete(logname);
    end
    diary('table45.txt');
    
    load EMP_D3_V2_data;
    nlogn = 0.14;
    b = 1;
    X = X_B(:,:,:,:,b);
    int = X(:, 1, :, :) .* X(:, 2, :, :);
    X = cat(2, X, int);
    [N, D, J, T] = size(X);
    y = y_B(:,:,:,b); % choice by each DMA, all 0 
    EyA = EyA_B(:,:,:,b); % market share
    Ey = Ey_B(:,:,:,b); % the b-th draw from E
    B = 1;
    
    %% step get dX and gamma fcn estimator
    
    %% convert multidim objects to flat 2D matrices for Y EyA Ey X
    dY_long = [];
    Xflat = [];
    Xflat2 = [];
    ind_new = [];
    dX = nan(N,D,J,T);
    B_first = 1;
    B_last = B_first + B - 1;
    
    for j = 1:J
        count = 0;
        for t = 1:T-1
            for s = (t+1):T
                count = count + 1;
                ind_new = [ind_new; [1:N]', ones(N,1) * j, ones(N,1) * t, ones(N,1) * s];
                dY_long = [dY_long; y(:,j,t)-y(:,j,s), EyA(:,j,t)-EyA(:,j,s),Ey(:,j,t) - Ey(:,j,s)];
                X_sub = X(:,1:3,:,:); 
                X_t = reshape(X_sub(:,:,:,t),N,[]);
                X_s = reshape(X_sub(:,:,:,s),N,[]);
                Xflat = [Xflat; X_t, X_s]; % Xflat doesnt change with j! only change with i and t
                dX(:,:,j,count) = X(:,:,j,t) - X(:,:,j,s);
    
                X_sub2 = X(:,1:2,:,:); % for lasso later
                X_t2 = reshape(X_sub2(:,:,:,t),N,[]);
                X_s2 = reshape(X_sub2(:,:,:,s),N,[]);
                Xflat2 = [Xflat2; X_t2, X_s2];
            end
        end
    end
    
     %% generate csv for R to run; replace it with every b cycle
     %     % get E(dy|X) for each product j using 1st order sieve
     T0 = size(dX, 4);
     Edy = nan(N*T0*J,1);
    
     for j = 1 : J
         j1 = (j-1) * N * T0 + 1 ;
         j2 = j * N * T0 ;
         Y_j = dY_long(j1:j2,2);
         X_j = Xflat2(j1:j2,:);
         X_j_med = repmat(median(X_j, 1),[N*T0,1]); % 2*D*J columns
         X_j_med_bas = subplus(X_j - X_j_med).^2;
         %       X_j_25 = repmat(quantile(X_j, .25, 1),[N*T0,1]);
         %       X_j_75 = repmat(quantile(X_j, .75, 1),[N*T0,1]);
         X_j_cross = double(create_interaction_variables(X_j,'all',2));
         %        X_j_sieve = [X_j_cross, X_j.^2, X_j_med_bas];
         X_j_sieve = [X_j_cross, X_j.^2];
         [L,FitInfo] = lasso(X_j_sieve,Y_j,'CV',10);
         idxLambdaMinMSE = FitInfo.IndexMinMSE;
         coef = L(:,idxLambdaMinMSE);
         coef0 = FitInfo.Intercept(idxLambdaMinMSE);
         Edy(j1:j2) = X_j_sieve * coef + coef0;
     end
    
     %% collect results
    dX_B(:,:,:,:,b) = dX; %#ok<*PFOUS>
    dEy_B(:,:,:,b) = permute(reshape(Edy,[N,T0,J]),[1,3,2]);
    
    dX_B = cat(4,dX_B,-dX_B);
    dEy_B = cat(3,dEy_B,-dEy_B);
    
    dXEname = sprintf('dXdE_%d_%d',B_first,B_last);
    save(dXEname,'dX_B','dEy_B');
    clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B nlogn rngSeed
    
    %% estimate parameters
    M_Step = 50; % Number of points along each dimension
    Tol = 1e-2;
    %nlogn = 0;
    
    Out_Buffer_Polar = 1;  % Select an integer
    Extra_Buffer_Polar = 1;
    
    b_wrong = [];
    
    check_B = nan(B,1);
    
    theta_B = nan(3,D-1,B);
    
    theta_out_B = nan(2,D-1,B);
    
    beta_B = nan(3,D,B);
    
    beta_out_B = nan(2,D,B);
    
    b = 1;
    
    agg_theta_q = [];
    
    %Qmax_test = -1;
    %Qmin_test = +Inf;
    
    dX = dX_B(:,:,:,:,b);
    dE = dEy_B(:,:,:,b);
    
    range_control = 1 ;
    %determine the range: 1 means (-pi, pi], 2 means (0, 2*pi]
    
    range_l = -pi ;
    range_u = pi ;
    theta_l = [-0.5*pi, range_l] ;
    theta_u = [0.5*pi, range_u] ;
    
    Qmax = -1;
    Qmin = +Inf;
    
    
    Loop_count = 0;
    Loop_control = 1;
    Tol_polar = 2 * pi;
    theta_scale_ratio = 0;
    
    tic
    count_1 = 0
    while (Loop_control == 1)
        count_1 = count_1 + 1;
        Loop_count = Loop_count + 1
        % cut range of theta by M_step parts evenly
        Tol_polar = [(theta_u(1) - theta_l(1)) / (M_Step - 1), (theta_u(2) - theta_l(2)) / M_Step]
        
        theta = theta_l + [0:1:(M_Step-1)]' * Tol_polar;
        
        
        [Theta1, Theta2] = ndgrid(theta(:,1), theta(:,2));
        
        theta1_gr = reshape(Theta1,[],1);
        theta2_gr = reshape(Theta2,[],1);
        theta_gr = [theta1_gr, theta2_gr];
        
        Length_L1 = M_Step^(D-1); % total number of points in the target field
        
        b_gr = nan(Length_L1, D);
        b_gr = [cos(theta1_gr).*cos(theta2_gr), cos(theta1_gr).*sin(theta2_gr), sin(theta1_gr)];
        
        % clear theta theta1_gr theta2_gr Theta1 Theta2
        % evaluate these points using Qfunc
        Q_gr = nan(Length_L1,1);


        parfor k = 1 : Length_L1
            stream = RandStream('Threefry', 'Seed', rngSeed + k + 4000000); % Unique seed per iteration
            RandStream.setGlobalStream(stream);  


            Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
           % Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)', localStream);

        end
        
        agg_theta_q = [agg_theta_q; theta1_gr, theta2_gr, Q_gr];
        
        Qmax_new = max(Q_gr);
        Qmin_new = min(Q_gr);
        
        %if theta_scale_ratio < 0.1
        
        Qmax = max(Qmax, Qmax_new); % don't need max() by construction?
        Qmin = min(Qmin, Qmin_new);
        
        ind_qt = find(Q_gr <= quantile(Q_gr, 0.2));
        
        theta_gr_new = theta_gr(ind_qt, :);
        theta_l_new = min(theta_gr_new, [], 1);
        theta_u_new = max(theta_gr_new, [], 1); % these are boundary of points that fall in lower 25% quantile of Q
        
        ind_min= find(Q_gr == Qmin_new);
        theta_min = theta_gr(ind_min,:);
        
        theta_min_l = min(theta_min, [], 1);
        theta_min_u = max(theta_min, [], 1); % these are boundary of points that correspond to min of Q
        
        theta_scale = max(theta_min_u - theta_min_l);
        theta_scale_ratio = theta_scale / max(theta_u - theta_l); % compare size of Qmin region against whole region of theta Q < 0.25
        
        
        if Loop_count == 3 && (abs(theta_l_new(2)-range_l) <= pi/8 || abs(theta_u_new(2)-range_u) <= pi/8 )
            Loop_count = 0
            range_l = 0 ;
            range_u = 2*pi ;
            theta_l = [-0.5*pi, range_l] ;
            theta_u = [0.5*pi, range_u] ;
            range_control = 2;
            continue
        end
    %     
        
        
        if Loop_count >= 3 && Qmin > 0 % if theta size ratio is large, break out of loop
            fprintf('Loop 1 finished at Step %d \r', Loop_count)
            
            Qmin
            
            Loop_control = 2;
            b_min = b_gr(ind_min,:);
            
            break;
            
        elseif Loop_count >= 3 && Qmin == 0
            
            fprintf('Loop 1 finished at Step %d \r', Loop_count)
            Qmin
            Loop_control = 3;
            break;
            
        end
        
        if Loop_count >= 50 || count_1 >= 10
            Loop_control = 2;
            b_min = b_gr(ind_min,:);
            
            break;
        end
        
       
        theta_l =  theta_l_new; % shrink theta region to points that satisfy Q < .25
        theta_u =  theta_u_new;
        
        
    end
    
    ind_notmin = find(Q_gr ~= Qmin);
    theta_notmin = theta_gr(ind_notmin,:);
    theta_notmin_l = min(theta_notmin, [], 1);
    theta_notmin_u = max(theta_notmin, [], 1)
    
    theta_l_del = theta_min_l - theta_notmin_l;
    theta_u_del = theta_notmin_u - theta_min_u;
    theta_del = min(min(theta_l_del), min(theta_u_del));
    
    theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
    theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
    Tol_Step = min(Tol_polar) / 5; % this tol setup
    
    
    
    while Loop_control == 2
        
        Loop_count = Loop_count + 1 %initialize loop count for loopcontrol 2?
        
        theta_l_out = [max(theta_min_l(1) - (Out_Buffer_Polar + theta_l_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), -0.5 * pi), max(theta_min_l(2) - (Out_Buffer_Polar + theta_l_del_0(2) * Extra_Buffer_Polar) * Tol_polar(2), range_l)]
        theta_u_out = [min(theta_min_u(1) + (Out_Buffer_Polar + theta_u_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), 0.5 * pi), min(theta_min_u(2) + (Out_Buffer_Polar + theta_u_del_0(2) * Extra_Buffer_Polar)  * Tol_polar(2), range_u)]
        
        theta_l_in = theta_min_l %initialize inner circle
        theta_u_in = theta_min_u
        
        M_Step_L2 = floor(max((theta_u_out - theta_l_out) / Tol_Step)) + 1; %adaptive number of parts cut
        Tol_polar = [(theta_u_out(1) - theta_l_out(1)) / (M_Step_L2 - 1), (theta_u_out(2) - theta_l_out(2)) / M_Step_L2]
        
        
        theta = theta_l_out +  [0:1:(M_Step_L2-1)]' * Tol_polar ;
        
        [Theta1, Theta2] = ndgrid(theta(:,1), theta(:,2));
        
        theta1_gr = reshape(Theta1,[],1);
        theta2_gr = reshape(Theta2,[],1);
        theta_gr = [theta1_gr,theta2_gr];
        
        % clear theta theta1_gr theta2_gr Theta1 Theta2 % ind_temp1 ind_temp2 ind_temp
        
        b_gr = [cos(theta_gr(:,1)).*cos(theta_gr(:,2)), cos(theta_gr(:,1)).*sin(theta_gr(:,2)), sin(theta_gr(:,1))];
        
        Length_L2 = size(theta_gr,1);
        Q_gr = nan(Length_L2,1);
        

       
        parfor k = 1 : Length_L2
            stream = RandStream('Threefry', 'Seed', rngSeed + k + 5000000); % Unique seed per iteration
            RandStream.setGlobalStream(stream);  
      
      
            Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)');
        end
        
        agg_theta_q = [agg_theta_q; theta1_gr, theta2_gr, Q_gr];
        
        Qmax_new = max(Q_gr);
        Qmin_new = min(Q_gr);
        
        ind_min = find(Q_gr == Qmin_new);
        ind_notmin = find(Q_gr ~= Qmin_new);
        
        theta_min = theta_gr(ind_min,:);
        theta_min_l = min(theta_min, [], 1)
        theta_min_u = max(theta_min, [], 1)
        
        theta_notmin = theta_gr(ind_notmin,:);
        theta_notmin_l = min(theta_notmin, [], 1)
        theta_notmin_u = max(theta_notmin, [], 1)
        
        theta_l_del = theta_min_l - theta_notmin_l;
        theta_u_del = theta_notmin_u - theta_min_u;
        theta_del = min(min(theta_l_del), min(theta_u_del));
        
        theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
        theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
        
        if Qmin_new >= Qmin || Qmin_new == 0 % if next round Qmin of finer points is no smaller than current Qmin, enter final loop
            
            Qmin = min(Qmin, Qmin_new);
            fprintf('Loop 2 finished at Step %d \r', Loop_count)
            Tol_Step = Tol_Step / 2;
            Loop_control = 3;
            break;
            
        else
            
            Qmin = Qmin_new;
            
        end
        
        if Loop_count >= 50
            Qmin = min(Qmin, Qmin_new);
            fprintf('Loop 2 finished at Step %d \r', Loop_count)
            Tol_Step = Tol_Step / 2;
            Loop_control = 3;
            break;
        end
        
    end
    
    %     plot(theta_min(:,1),theta_min(:,2),'m.')
    
    
    theta_Qmin = theta_min;
    
    while Loop_control == 3
        
        Loop_count = Loop_count + 1
        
        theta_l_out = [max(theta_min_l(1) - (Out_Buffer_Polar + theta_l_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), -0.5 * pi), max(theta_min_l(2) - (Out_Buffer_Polar + theta_l_del_0(2) * Extra_Buffer_Polar) * Tol_polar(2), range_l)]
        theta_u_out = [min(theta_min_u(1) + (Out_Buffer_Polar + theta_u_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), 0.5 * pi), min(theta_min_u(2) + (Out_Buffer_Polar + theta_u_del_0(2) * Extra_Buffer_Polar)  * Tol_polar(2), range_u)]
        
        theta_l_in = theta_min_l; %initialize inner circle
        theta_u_in = theta_min_u;
        
        M_Step_L3 = floor((theta_u_out - theta_l_out) ./ Tol_polar) + 1; %adaptive number of parts cut
        
        for d = 1:(D-1)
            theta_cell{d} = theta_l_out(d) +  [0:1:(M_Step_L3(d)-1)]' * Tol_polar(d);
        end
        
        [Theta1, Theta2] = ndgrid(theta_cell{1}, theta_cell{2});
        
        theta1_gr = reshape(Theta1,[],1);
        theta2_gr = reshape(Theta2,[],1);
        theta_gr = [theta1_gr,theta2_gr];
        
        ind_temp1 = (min((theta_gr >= theta_l_in), [], 2) == 1);
        ind_temp2 = (min((theta_gr <= theta_u_in), [], 2) == 1);
        ind_temp = min(ind_temp1, ind_temp2);
        ind_keep = find(ind_temp ~= 1);
        
        Length_L3 = length(ind_keep);
        theta_gr = theta_gr(ind_keep,:);
        
        % clear theta theta1_gr theta2_gr Theta1 Theta2 % ind_temp1 ind_temp2 ind_temp 
        
        b_gr = [cos(theta_gr(:,1)).*cos(theta_gr(:,2)), cos(theta_gr(:,1)).*sin(theta_gr(:,2)), sin(theta_gr(:,1))];
        Q_gr = nan(Length_L3,1);
        

        parfor k = 1 : Length_L3
            stream = RandStream('Threefry', 'Seed', rngSeed + k + 6000000); % Unique seed per iteration
            RandStream.setGlobalStream(stream);  
          
            Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
        end
        
        agg_theta_q = [agg_theta_q; theta_gr, Q_gr];
        
        Qmax_new = max(Q_gr);
        Qmin_new = min(Q_gr);
        
        if Qmin_new < Qmin
            warning('smaller Q encountered')
            Qmin = Qmin_new;
        end
        
        ind_min = find(Q_gr == Qmin_new);
        ind_notmin = find(Q_gr ~= Qmin_new);
        
        theta_min_new = theta_gr(ind_min,:);
        theta_min_l_new = min(theta_min_new, [], 1)
        theta_min_u_new = max(theta_min_new, [], 1)
        
        theta_Qmin = [theta_Qmin; theta_min_new];
        
        theta_min_l = min(theta_Qmin,[],1);
        theta_min_u = max(theta_Qmin,[],1);
        
        theta_notmin = theta_gr(ind_notmin,:);
        theta_notmin_l = min(theta_notmin, [], 1);
        theta_notmin_u = max(theta_notmin, [], 1)
        
        theta_l_del = theta_min_l - theta_notmin_l;
        theta_u_del = theta_notmin_u - theta_min_u;
        theta_del = min(min(theta_l_del), min(theta_u_del));
        
        theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
        theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
        
        if theta_del > 0 && Tol_Step <= Tol
            Tol_Final = Tol_Step
            Loop_control = 4;
            break;
        end
        
        if theta_del > 0
            Extra_Buffer_Polar = 1;
            Tol_polar = max(Tol_polar /2, Tol);
            Tol_Step = max(Tol_polar)
        else
            Extra_Buffer_Polar = Extra_Buffer_Polar  + 1;
        end
        
        if Loop_count >= 100
            
            break;
        end
        
        
        
    end
    
    ind_nlogn = find(agg_theta_q(:,3) <= Qmin_new + nlogn);
    theta_nlogn = agg_theta_q(ind_nlogn,1:2);
    theta_nlogn_min = min(theta_nlogn, [], 1);
    theta_nlogn_max = max(theta_nlogn, [], 1);
    
    % Q0 = Qfunc(dX, dE, beta0);
    
    theta_B(:,:,b) = [theta_min_l; mean(theta_Qmin,1);theta_min_u];
    
    beta_Qmin = [cos(theta_nlogn(:,1)).*cos(theta_nlogn(:,2)), cos(theta_nlogn(:,1)).*sin(theta_nlogn(:,2)), sin(theta_nlogn(:,1))];
    b_notmin = [cos(theta_notmin(:,1)).*cos(theta_notmin(:,2)), cos(theta_notmin(:,1)).*sin(theta_notmin(:,2)), sin(theta_notmin(:,1))];
    
    
    beta_B(:,:,b) = [min(beta_Qmin,[],1); 1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
    
    %b_inbound = [min(b_Qmin,[],1); max(b_Qmin,[],1)];
    beta_out_B(:,:,b) = [min(b_notmin,[],1); max(b_notmin,[],1)];
    theta_out_B(:,:,b) = [theta_notmin_l; theta_notmin_u];
    
    %    num_b_wrong = length(b_wrong);
    
    clearvars -except N D J T B theta_B theta_out_B beta_B beta_out_B check_B B_first B_last nlogn agg_theta_q beta_Qmin rngSeed
    eval(sprintf('save Sim_%dD_results_%d_%d_%d.mat',D,B_first,B_last,nlogn))
    
    %% Performance Evaluation
    beta_B
    beta_B_nlogn = beta_B;
    save('beta_B_nlogn.mat', 'beta_B_nlogn');
    diary off
    
    poolobj = gcp('nocreate');
    evalc('delete(poolobj)');

end