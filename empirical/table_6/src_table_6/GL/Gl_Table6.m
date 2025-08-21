function Gl_Table6(varargin)

    %% Instruction of setting number of alpha and computational cores
    %
    % Default: alpha = 0.15, core = 4
    %
    %   Gl_Table6() runs with default parameters (alpha=0.15, cores=4)
    %   Gl_Table6(alpha) specifies alpha value (cores=4)
    %   Gl_Table6(alpha, core) specifies both alpha and core count
    %
    % Input Parameters:
    %   alpha - Learning rate parameter (default: 0.15)
    %   core - Number of parallel workers (default: 4)
    %
    % Example:
    %   Gl_Table6()               % Uses alpha=0.15, cores=4
    %   Gl_Table6(0.20)           % Uses alpha=0.20, cores=4
    %   Gl_Table6(0.20, 8)        % Uses alpha=0.20, cores=8
    % Change core value if your system has more available cores
    %
    % To check available cores: type "feature('numcores')" or
    % "maxNumCompThreads" in matlab command
    %
    %% WARNING: Setting number of core higher than physical cores will cause error
    %
    % Last Change Date: 17/May/2025


    %% Defaults
    alpha = 0.15;
    core = 4;
    
    %% Handle inputs
    if nargin > 0, alpha = varargin{1}; end
    if nargin > 1, core = varargin{2}; end
    
    %% Initialize parallel pool
    if isempty(gcp('nocreate'))
        parpool(core); 
    end

    %% Run analysis (replace with your actual code)
    fprintf('Running GL_table 6 with alpha=%.2f, cores=%d\n', alpha, core);

    %alpha = .15;
    N = 205; % # of individuals in the market
    D = 3;  % # of dimension of observable product characteristics (>= 3)
    J = 4;  % # of alternatives including outside option 0 (>= 3)
    T = 10;  % # of time periods (>= 2)
    B_first = 1;
    B = 1000;
    
    logname = sprintf("GL_Table6_N=%d_D=%d_J=%d_T=%d_alpha=%d.log",N,D,J,T,alpha);
    
    if exist(logname)
        diary off
        delete(logname);
    end
    
    diary(logname)
 
    
    B_last = B_first + B - 1;
    
    
    %% DGP
    rngSeed = 123;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream


    beta0 = 2 * [-2; repmat(1,[D-1,1])];


    %beta0norm = beta0 / norm(beta0);
    % 
    % if D == 3
    %     theta0 = [asin(beta0norm(2)), atan(beta0norm(2)/beta0norm(1))]
    % elseif D==4
    %     theta0 = [asin(beta0norm(3)), atan(beta0norm(3)/sqrt(beta0norm(1)^2 +beta0norm(2)^2) ), atan(beta0norm(2)/beta0norm(1))]
    % end
    
    
    X_B = nan(N,D,J,T,B);
    y_B = nan(N,J,T,B);
    EyA_B = nan(N,J,T,B);
    AXbeta0_B = nan(N,J,T,B);
    
    %rng(B_first);
    %spmd
        %rng(B_first + labindex - 1, 'twister');
    %end
%%%%%%%%%%%% change 29 may
    %% Initialize single stream for entire script
    % rngSeed = 42;  % Master seed
    % stream = parallel.pool.Constant(RandStream('Threefry', 'Seed', rngSeed));
    % 

    parfor b = B_first:B_last
         % localStream = stream.Value;
         % localStream.Substream = b;
        %rng(b); 
        stream = RandStream('Threefry', 'Seed', rngSeed + b); % Unique seed per iteration
        RandStream.setGlobalStream(stream);  
    
        % Z_i = (2 * sqrt(3) * rand(N,1)) - sqrt(3); % latent variable to induce correlation between X_i2jt and A_i2
        Z = rand(N,J);
        X = nan(N,D,J,T);
        %x1 = price x2 = promo x3 = price * promo
        %    X(:, 1, :, :) = 2 * randi(2, N, J, T) - 3 ;
        %    X(:, 2, :, :) = 2 * sqrt(2*sigma) * rand(N,J,T) - sqrt(6) + repmat(Z_i, 1, J, T);
        % X(:, [3:D], :, :) = 2 * rand(N, D-2, J, T) - 1;
        X(:, 1, :, :) = 4 * rand(N,J,T);
        X(:,2,:,:) = (1-alpha) * rand(N,J,T) + (alpha) * repmat(Z,1,1,T);
        X(:, 3, :, :) = X(:, 1, :, :) .* X(:, 2, :, :);
        % X(:,4:D,:,:) = zeros(N,D-3,J,T);
        % create positive correlation between X_promo and A_scale
        A_scale = Z + 1;
        %     A_location = 0.1 * rand(N, J);
        A_location = zeros(N, J); % shut down location A2, so the effect is purely from A1 correlated with X2
    
        %A_scale = 0.5 * rand(N,1) + 2;
        %     A_location(:, 1) = zeros(N,1);
        %     A_location(:, 2) = subplus(Z_i);
        %     A_location(:, [3:J]) = 0.5 * rand(N,J-2) - 0.25;
    
        Xbeta0 = permute(sum(X .* repmat(beta0', [N,1,J,T]), 2),[1,3,4,2]);
        A1 = repmat(A_scale, 1, 1, T);
        A2 = repmat(A_location, 1, 1, T);
    
        AXbeta0 = A1 .* (Xbeta0 + A2);
    
        % check 
        % 
        % 
        % 
        eps = evrnd(0,1,[N,J,T]); %
       % eps = -log(-log(rand(localStream, N, J, T)));

        iu = AXbeta0 + eps; % again N*J*T matrix
    
        % calculate y vector by extracting the max of each row and make it 1
        y = double(bsxfun(@eq,iu,max(iu,[],2)));
    
        % calculate true market shares
        expind = exp(AXbeta0);
        sum_expind = sum(expind, 2);
        EyA = expind ./sum_expind;
    
        %  clear eps iu expind sum_expind
    
        X_B(:,:,:,:,b) = X;
        y_B(:,:,:,b) = y;
        EyA_B(:,:,:,b) = EyA;
        AXbeta0_B(:,:,:,b) =AXbeta0;
    end
    
    
    y_B = EyA_B;% in emprical data we only have EyA_B, so in this exercise we
    % keep EyA_B only and all the results based on y_B should be interpreted as EyA_B.
    % if this is not desirable, comment this line.
    
    
    %clearvars -except N D J T B B_first B_last core beta0 theta0 X_B y_B EyA_B AXbeta0_B alpha nlogn beta0
    %eval(sprintf('save Sim_1data_%d_%d_%d.mat',B_first,B_last,alpha))
    
    %% get dXdE matrices
    
    Ey_B = zeros(N,J,T,B);
    B0 = B;
    T0 = T*(T-1)/2;
    
    dX_B = nan(N,D,J,T0,B0);
    dEy_B = nan(N,J,T0,B0);
    %dEyA_B = nan(N,J,T0,B0);
    %dETrues_B = nan(N,J,T0,B0);
    
    %Xflat_name = sprintfc('X%d', 1:2*D*J);
    %ind_name = {'N_var','J_var','T_var','S_var','dy','dEyA','dEy'};
    %headers = horzcat(ind_name,Xflat_name);
    
    %% step 1 get gamma fcn estimator
 
    parfor b = 1: B
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 10000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 
        %rng(b); %wenli 2
        X = X_B(:,1:2,:,:,b);
        X1 = X_B(:,:,:,:,b);
        y = y_B(:,:,:,b);
        EyA = EyA_B(:,:,:,b);
        Ey = Ey_B(:,:,:,b);
    %    dataname = sprintf('data_%d.csv', b + B_first - 1);
     %   importname = sprintf('NNet_%d.csv',b + B_first - 1);
    
        %% convert multidim objects to flat 2D matrices for Y EyA Ey X
        dY_long = [];
        Xflat = [];
        ind_new = [];
        dX = nan(N,D,J,T0);
    
        for j = 1:J
            count = 0;
            for t = 1:T-1
                for s = (t+1):T
                    count = count + 1;
                    ind_new = [ind_new; [1:N]', ones(N,1) * j, ones(N,1) * t, ones(N,1) * s];
                    dY_long = [dY_long; y(:,j,t)-y(:,j,s), EyA(:,j,t)-EyA(:,j,s),Ey(:,j,t) - Ey(:,j,s)];
                    X_t = reshape(X(:,:,:,t),N,[]);
                    X_s = reshape(X(:,:,:,s),N,[]);
                    Xflat = [Xflat; X_t, X_s]; % Xflat doesnt change with j! only change with i and t
                    dX(:,:,j,count) = X1(:,:,j,t) - X1(:,:,j,s);
                end
            end
        end
    
        Edy = nan(N*T0*J,1);
    
        for j = 1 : J
            j1 = (j-1) * N * T0 + 1 ;
            j2 = j * N * T0 ;
            Y_j = dY_long(j1:j2,1);
            X_j = Xflat(j1:j2,:);
            % X_j_med = repmat(median(X_j, 1),[N*T0,1]); % 2*D*J columns
            % X_j_med_bas = subplus(X_j - X_j_med).^2;
            %       X_j_25 = repmat(quantile(X_j, .25, 1),[N*T0,1]);
            %       X_j_75 = repmat(quantile(X_j, .75, 1),[N*T0,1]);
            X_j_cross = double(create_interaction_variables(X_j,'all',2));
            %        X_j_sieve = [X_j_cross, X_j.^2, X_j_med_bas];
            X_j_sieve = [X_j_cross, X_j.^2];
            %
            %opts = statset('UseParallel', false, 'Streams', RandStream('twister', 'Seed', b));
           
             opts = statset('UseParallel', false, 'Streams', stream);
            [L, FitInfo] = lasso(X_j_sieve, Y_j, 'CV', 10, 'Options', opts);
    

            
            %[L,FitInfo] = lasso(X_j_sieve,Y_j,'CV',10);
            idxLambdaMinMSE = FitInfo.IndexMinMSE;
            coef = L(:,idxLambdaMinMSE);
            coef0 = FitInfo.Intercept(idxLambdaMinMSE);
            Edy(j1:j2) = X_j_sieve * coef + coef0;
            %mdl_j = fitlm(X_j, Y_j,'quadratic');
            %Edy(j1:j2) = predict(mdl_j,X_j);
        end
    
    
        dX_B(:,:,:,:,b) = dX; %#ok<*PFOUS>
        dEy_B(:,:,:,b) = permute(reshape(Edy,[N,T0,J]),[1,3,2]);
    end
    
    dX_B = cat(4,dX_B,-dX_B);
    dEy_B = cat(3,dEy_B,-dEy_B);
    %dEyA_B = cat(3,dEyA_B,-dEyA_B);
    
    %dXEname = sprintf('dXdE_%d_%d',B_first,B_last);
    %save(dXEname,'dX_B','dEy_B');
    %clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B alpha nlogn
    %delete *.csv
    
    %% Sim_3D_parGrid
    
    beta0 = 2 * [-2; repmat(1,[D-1,1])];
    
    beta0 = beta0 /norm(beta0);
    
    % if beta0(1)>=0
    %     theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))];
    % elseif beta0(1)<0 && beta0(2)>=0
    %     theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))+pi];
    % else
    %     theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))-pi];
    % end
    
    M_Step = 50; % Number of points along each dimension
    Tol = 1e-2;
    
    Out_Buffer_Polar = 1;  % Select an integer
    
    
    b_wrong = [];
    
    check_B = nan(B,1);
    
    theta_B = nan(3,D-1,B);
    
    theta_out_B = nan(2,D-1,B);
    
    beta_B = nan(3,D,B);
    
    beta_out_B = nan(2,D,B);
  %%%%%% may 29

    parfor b = 1:B
        %rng(b); %wenli 4
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 20000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 
        Extra_Buffer_Polar = 1;
        fprintf('\nSimulation %d\n', b)
    
    
        %    Qmax_test = -1;
        %    Qmin_test = +Inf;
    
        dX = dX_B(:,:,:,:,b);
        dE = dEy_B(:,:,:,b);
    
        % b_u = ones(1,D);
        % b_l = -ones(1,D);
        %
        % b_min_l = [];
        % b_min_u = [];
        %
        % b_dif = 2;
        % b_dif_del = 2;
    
    
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
    
        %     tic
    
        while (Loop_control == 1)
    
            Loop_count = Loop_count + 1;
            % cut range of theta by M_step parts evenly
            Tol_polar = [(theta_u(1) - theta_l(1)) / (M_Step - 1), (theta_u(2) - theta_l(2)) / M_Step];
    
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
            for k = 1 : Length_L1
                Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
            end
    
            Qmax_new = max(Q_gr);
            Qmin_new = min(Q_gr);
    
            %if theta_scale_ratio < 0.1
    
            Qmax = max(Qmax, Qmax_new); % don't need max() by construction?
            Qmin = min(Qmin, Qmin_new);
    
            ind_qt = find(Q_gr <= quantile(Q_gr, 0.20));
    
            theta_gr_new = theta_gr(ind_qt, :);
            theta_l_new = min(theta_gr_new, [], 1);
            theta_u_new = max(theta_gr_new, [], 1); % these are boundary of points that fall in lower 25% quantile of Q
    
            ind_min= find(Q_gr == Qmin_new);
            theta_min = theta_gr(ind_min,:);
    
            theta_min_l = min(theta_min, [], 1);
            theta_min_u = max(theta_min, [], 1); % these are boundary of points that correspond to min of Q
    
            theta_scale = max(theta_min_u - theta_min_l);
            theta_scale_ratio = theta_scale / max(theta_u - theta_l); % compare size of Qmin region against whole region of theta Q < 0.25
    
    
            if all([Loop_count == 3, any([abs(theta_l_new(2)-range_l) <= pi/8, abs(theta_u_new(2)-range_u) <= pi/8])])
                Loop_count = 0;
                range_l = 0 ;
                range_u = 2*pi ;
                theta_l = [-0.5*pi, range_l] ;
                theta_u = [0.5*pi, range_u] ;
                range_control = 2;
                continue
            end
    
            %             Qmin
    
    
            if all([Loop_count >= 3, Qmin > 0]) % if theta size ratio is large, break out of loop
                %                 fprintf('Loop 1 finished at Step %d \r', Loop_count)
    
                %                 Qmin
    
                Loop_control = 2;
                b_min = b_gr(ind_min,:);
    
                break;
    
            elseif all([Loop_count >= 3, Qmin == 0])
    
                %                 fprintf('Loop 1 finished at Step %d \r', Loop_count)
                %                 Qmin
                Loop_control = 3;
                break;
    
            end
    
            theta_l =  theta_l_new; % shrink theta region to points that satisfy Q < .25
            theta_u =  theta_u_new;
    
        end
    
        ind_notmin = find(Q_gr ~= Qmin);
        theta_notmin = theta_gr(ind_notmin,:);
        theta_notmin_l = min(theta_notmin, [], 1);
        theta_notmin_u = max(theta_notmin, [], 1);
    
        theta_l_del = theta_min_l - theta_notmin_l;
        theta_u_del = theta_notmin_u - theta_min_u;
        theta_del = min(min(theta_l_del), min(theta_u_del));
    
        theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
        theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
        Tol_Step = min(Tol_polar) / 5; % this tol setup
    
    
        while Loop_control == 2
    
            Loop_count = Loop_count + 1; %initialize loop count for loopcontrol 2?
    
            theta_l_out = [max(theta_min_l(1) - (Out_Buffer_Polar + theta_l_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), -0.5 * pi), max(theta_min_l(2) - (Out_Buffer_Polar + theta_l_del_0(2) * Extra_Buffer_Polar) * Tol_polar(2), range_l)];
            theta_u_out = [min(theta_min_u(1) + (Out_Buffer_Polar + theta_u_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), 0.5 * pi), min(theta_min_u(2) + (Out_Buffer_Polar + theta_u_del_0(2) * Extra_Buffer_Polar)  * Tol_polar(2), range_u)];
    
            theta_l_in = theta_min_l; %initialize inner circle
            theta_u_in = theta_min_u;
    
            M_Step_L2 = floor(max((theta_u_out - theta_l_out) / Tol_Step)) + 1; %adaptive number of parts cut
            Tol_polar = [(theta_u_out(1) - theta_l_out(1)) / (M_Step_L2 - 1), (theta_u_out(2) - theta_l_out(2)) / M_Step_L2];
    
    
            theta = theta_l_out +  [0:1:(M_Step_L2-1)]' * Tol_polar ;
    
            [Theta1, Theta2] = ndgrid(theta(:,1), theta(:,2));
    
            theta1_gr = reshape(Theta1,[],1);
            theta2_gr = reshape(Theta2,[],1);
            theta_gr = [theta1_gr,theta2_gr];
    
            % clear theta theta1_gr theta2_gr Theta1 Theta2 % ind_temp1 ind_temp2 ind_temp
    
            b_gr = [cos(theta_gr(:,1)).*cos(theta_gr(:,2)), cos(theta_gr(:,1)).*sin(theta_gr(:,2)), sin(theta_gr(:,1))];
    
            Length_L2 = size(theta_gr,1);
            Q_gr = nan(Length_L2,1);
    
    
            for k = 1 : Length_L2
                Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
            end
    
            Qmax_new = max(Q_gr);
            Qmin_new = min(Q_gr);
    
            ind_min = find(Q_gr == Qmin_new);
            ind_notmin = find(Q_gr ~= Qmin_new);
    
            theta_min = theta_gr(ind_min,:);
            theta_min_l = min(theta_min, [], 1);
            theta_min_u = max(theta_min, [], 1);
    
            theta_notmin = theta_gr(ind_notmin,:);
            theta_notmin_l = min(theta_notmin, [], 1);
            theta_notmin_u = max(theta_notmin, [], 1) ;
    
            theta_l_del = theta_min_l - theta_notmin_l;
            theta_u_del = theta_notmin_u - theta_min_u;
            theta_del = min(min(theta_l_del), min(theta_u_del));
    
            theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
            theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
            if any([Qmin_new >= Qmin , Qmin_new == 0]) % if next round Qmin of finer points is no smaller than current Qmin, enter final loop
    
                Qmin = min(Qmin, Qmin_new);
                %             fprintf('Loop 2 finished at Step %d \r', Loop_count)
                Tol_Step = Tol_Step / 2;
                Loop_control = 3;
                break;
    
            else
    
                Qmin = Qmin_new;
    
            end
    
        end
    
        %     plot(theta_min(:,1),theta_min(:,2),'m.')
    
    
        theta_Qmin = theta_min;
    
        while Loop_control == 3
    
            Loop_count = Loop_count + 1;
    
            theta_l_out = [max(theta_min_l(1) - (Out_Buffer_Polar + theta_l_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), -0.5 * pi), max(theta_min_l(2) - (Out_Buffer_Polar + theta_l_del_0(2) * Extra_Buffer_Polar) * Tol_polar(2), range_l)];
            theta_u_out = [min(theta_min_u(1) + (Out_Buffer_Polar + theta_u_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), 0.5 * pi), min(theta_min_u(2) + (Out_Buffer_Polar + theta_u_del_0(2) * Extra_Buffer_Polar)  * Tol_polar(2), range_u)];
    
            theta_l_in = theta_min_l; %initialize inner circle
            theta_u_in = theta_min_u;
    
            M_Step_L3 = floor((theta_u_out - theta_l_out) ./ Tol_polar) + 1; %adaptive number of parts cut
    
            % for d = 1:(D-1)
            %     theta_cell(d) = theta_l_out(d) +  [0:1:(M_Step_L3(d)-1)]' * Tol_polar(d);
            % end
            %
    
            theta_cell1 = theta_l_out(1) +  [0:1:(M_Step_L3(1)-1)]' * Tol_polar(1);
            theta_cell2 = theta_l_out(2) +  [0:1:(M_Step_L3(2)-1)]' * Tol_polar(2);
    
            [Theta1, Theta2] = ndgrid(theta_cell1, theta_cell2);
    
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
    
            for k = 1 : Length_L3
                Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
            end
    
            Qmax_new = max(Q_gr);
            Qmin_new = min(Q_gr);
    
            if Qmin_new < Qmin
                %             warning('smaller Q encountered')
                Qmin = Qmin_new;
            end
    
            ind_min = find(Q_gr == Qmin_new);
            ind_notmin = find(Q_gr ~= Qmin_new);
    
            theta_min_new = theta_gr(ind_min,:);
            theta_min_l_new = min(theta_min_new, [], 1);
            theta_min_u_new = max(theta_min_new, [], 1);
    
            theta_Qmin = [theta_Qmin; theta_min_new];
    
            theta_min_l = min(theta_Qmin,[],1);
            theta_min_u = max(theta_Qmin,[],1);
    
            theta_notmin = theta_gr(ind_notmin,:);
            theta_notmin_l = min(theta_notmin, [], 1);
            theta_notmin_u = max(theta_notmin, [], 1);
    
            theta_l_del = theta_min_l - theta_notmin_l;
            theta_u_del = theta_notmin_u - theta_min_u;
            theta_del = min(min(theta_l_del), min(theta_u_del));
    
            theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
            theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
            if all([theta_del > 0, Tol_Step <= Tol])
                Tol_Final = Tol_Step;
                Loop_control = 4;
                break;
            end
    
            if theta_del > 0
                Extra_Buffer_Polar = 1;
                Tol_polar = max(Tol_polar /2, Tol);
                Tol_Step = max(Tol_polar);
            else
                Extra_Buffer_Polar = Extra_Buffer_Polar  + 1;
            end
    
        end
    
    
    
        Q0 = Qfunc(dX, dE, beta0);
    
        theta_B(:,:,b) = [theta_min_l; mean(theta_Qmin,1);theta_min_u];
    
        beta_Qmin = [cos(theta_Qmin(:,1)).*cos(theta_Qmin(:,2)), cos(theta_Qmin(:,1)).*sin(theta_Qmin(:,2)), sin(theta_Qmin(:,1))];
    
        beta_B(:,:,b) = [min(beta_Qmin,[],1); 1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
    
    end
    
    num_b_wrong = length(b_wrong);
    
    clearvars -except N D J T B beta0 theta0 theta_B theta_out_B beta_B beta_out_B check_B B_first B_last alpha nlogn
    eval(sprintf('save Sim_%dD_results_%d_%d_%d.mat',D,B_first,B_last,alpha))
    
    %% Performance Evaluation
    beta0
    
    beta_m = permute(beta_B(2,:,:),[2 3 1]);
    
    sign_beta_m = prod(subplus(beta_m .* repmat(beta0,[1 B])),1);
    
    b_percent = length(find(sign_beta_m>0))/B
    
    
    
    % Create a single filename with alpha value
    filename = sprintf('GL_%s_beta_percent.mat', strrep(num2str(alpha), '.', '_'));
    
    % Save all variables in one .mat file
    save(filename, 'b_percent');
    
    % Optional: Display confirmation
    fprintf('Results saved to %s\n', filename);
    
    
    
    
    
    diary off
    
    poolobj = gcp('nocreate');
    evalc('delete(poolobj)');

end
