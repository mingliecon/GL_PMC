function MC_AGG_20230909_nlogn(varargin)

    %----------------------------------------------------------------------
    % File: MC_AGG.m
    
    % Description: Generate simulation data
    
    % Coded by: Yuan Gao & Ming Li
    % Last modified: 12/01/2018
    
    if ~isempty(varargin)
        temp = varargin;
        N = temp{1}*1000;
        D = temp{2};
        J = temp{3};
        T = temp{4};
        B_first = temp{5};
        B = temp{6};
        nlogn = temp{7};
        core = temp{8};
    
    else
        %default
        N = 10000; % # of individuals in the market
        D = 3;  % # of dimension of observable product characteristics (>= 3)
        J = 3;   % # of alternatives including outside option 0 (>= 3)
        T = 2;  % # of time periods (>= 2)
        B_first = 1;
        B = 1000;
        %core = 64;
        nlogn = 0;
        core = 4;
    end

    rngSeed = 791;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream
    
    
    
    
    evalc('parpool(core)');
    
    logname = sprintf("Table4_N=%d_D=%d_J=%d_T=%d_B=%d_nlogn=%d_bigger_index_nlogn.log",N,D,J,T,B,nlogn);
    
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
    
    
    %% DGP
    
    beta0norm = beta0 / norm(beta0);
    
    if D == 3
        theta0 = [asin(beta0norm(2)), atan(beta0norm(2)/beta0norm(1))];
    elseif D==4
        theta0 = [asin(beta0norm(3)), atan(beta0norm(3)/sqrt(beta0norm(1)^2 +beta0norm(2)^2) ), atan(beta0norm(2)/beta0norm(1))];
    end
    
    
    X_B = nan(N,D,J,T,B);
    y_B = nan(N,J,T,B);
    EyA_B = nan(N,J,T,B);
    %AXbeta0_B = [];

    
    parfor b = 1:B
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 120000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 

    
        Z_i = 2*sqrt(3)*rand(N,1)-sqrt(3); % latent variable to induce correlation between X_i2jt and A_i2
    
        X = nan(N,D,J,T);
    
        X(:, 1, :, :) =  2*rand([N, J, T]) - 1;
        X(:, 2, :, :) = repmat(Z_i, 1, J, T) + sqrt(6) * randn([N,J,T]);
        X(:, 3:D, :, :) = randn([N,J,T]);
    
        A_scale = 0.5 * rand(N,1) + 2;
    
        A_location = nan(N, J);
        A_location(:, 1) = zeros(N,1);
        A_location(:, 2) = subplus(Z_i);
        A_location(:, 3:J) = 0.5 * rand(N,J-2) - 0.25;
    
        Xbeta0 = permute(sum(X .* repmat(beta0', [N,1,J,T]), 2),[1,3,4,2]);
        A1 = repmat(A_scale, 1, J, T);
        A2 = repmat(A_location, 1, 1, T);
    
        AXbeta0 = A1 .* (Xbeta0 + A2);
    
        eps = evrnd(0,1,[N,J,T]); %
        iu = AXbeta0 + eps; % again N*J*T matrix
    
        % calculate y vector by extracting the max of each row and make it 1
        y = double(bsxfun(@eq,iu,max(iu,[],2)));
    
        % calculate true market shares
        expind = exp(AXbeta0);
        sum_expind = sum(expind, 2);
        EyA = expind ./sum_expind;
    
        %clear eps iu expind sum_expind
    
        X_B(:,:,:,:,b) = X;
        y_B(:,:,:,b) = y;
        EyA_B(:,:,:,b) = EyA;
        %AXbeta0_B = cat(4, AXbeta0_B,AXbeta0);
    end
    
    %clearvars -except N D J T B B_first B_last core beta0 theta0 X_B y_B EyA_B nlogn
    %eval(sprintf('save Sim_1data_%d_%d_%d',B_first,B_last,nlogn))
    
    
    Ey_B = zeros(N,J,T,B);
    B0 = B;
    T0 = T*(T-1)/2;
    
    dX_B = nan(N,D,J,T0,B0);
    dEy_B = nan(N,J,T0,B0);
    %dEyA_B = nan(N,J,T0,B0);
    dETrues_B = nan(N,J,T0,B0);
    
    Xflat_name = sprintfc('X%d', 1:2*D*J);
    ind_name = {'N_var','J_var','T_var','S_var','dy','dEyA','dEy'};
    headers = horzcat(ind_name,Xflat_name);
    
    %% step 1 get gamma fcn estimator
    
    parfor b = 1: B
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 130000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 
    
        X = X_B(:,:,:,:,b);
        y = y_B(:,:,:,b);
        EyA = EyA_B(:,:,:,b);
        Ey = Ey_B(:,:,:,b);
        %    dataname = sprintf('data_%d.csv', b + B_first - 1);
        %    importname = sprintf('NNet_%d.csv',b + B_first - 1);
    
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
                    dX(:,:,j,count) = X(:,:,j,t) - X(:,:,j,s);
                end
            end
        end
        % generate N_var T_var S_var J_var dY dEYA dEY Xflat a (N*J*T*T-1) by (7+2*D*J) matrix
        %    data_flat = [ind_new, dY_long, Xflat];
    
        %% generate csv for R to run; replace it with every b cycle
        %     % get E(dy|X) for each product j using 1st order sieve
    
        Edy = nan(N*T0*J,1);
    
        for j = 1 : J
            j1 = (j-1) * N * T0 + 1 ;
            j2 = j * N * T0 ;
            Y_j = dY_long(j1:j2,1);
            X_j = Xflat(j1:j2,:);
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
            %mdl_j = fitlm(X_j, Y_j,'quadratic');
            %Edy(j1:j2) = predict(mdl_j,X_j);
        end
    
        %% collect results
        dX_B(:,:,:,:,b) = dX; %#ok<*PFOUS>
        dEy_B(:,:,:,b) = permute(reshape(Edy,[N,T0,J]),[1,3,2]);
    
    end
    
    dX_B = cat(4,dX_B,-dX_B);
    dEy_B = cat(3,dEy_B,-dEy_B);
    %dEyA_B = cat(3,dEyA_B,-dEyA_B);
    
    %dXEname = sprintf('dXdE_%d_%d_%d',B_first,B_last,nlogn);
    %save(dXEname,'dX_B','dEy_B');
    %clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B nlogn
    %delete *.csv
    
    %% Sim_3D_parGrid
    
    beta0 = [2; ones(D-1,1)];
    beta0 = beta0 /norm(beta0);
    
    
    if beta0(1)>=0
        theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))];
    elseif beta0(1)<0 && beta0(2)>=0
        theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))+pi];
    else
        theta0 = [asin(beta0(3)), atan(beta0(2)/beta0(1))-pi];
    end
    
    M_Step = 50; % Number of points along each dimension
    Tol = 1e-2;
    
    Out_Buffer_Polar = 1;  % Select an integer
    
    
    b_wrong = [];
    
    check_B = nan(B,1);
    
    theta_B = nan(3,D-1,B);
    
    theta_out_B = nan(2,D-1,B);
    
    beta_B = nan(3,D,B);
    
    beta_out_B = nan(2,D,B);
    
    
    parfor b = 1:B
    
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 141000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 
    
    
    
        agg_theta_q = [];
        Extra_Buffer_Polar = 1;
        %  fprintf('\nSimulation %d\n', b)
    
    
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
    
            agg_theta_q = [agg_theta_q; theta1_gr, theta2_gr, Q_gr];
    
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
    
            agg_theta_q = [agg_theta_q; theta1_gr, theta2_gr, Q_gr];
    
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
    
            agg_theta_q = [agg_theta_q; theta_gr, Q_gr];
    
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
    
        ind_nlogn = find(agg_theta_q(:,3) <= Qmin_new + nlogn);
        theta_nlogn = agg_theta_q(ind_nlogn,1:2);
        %theta_nlogn_min = min(theta_nlogn, [], 1);
        %theta_nlogn_max = max(theta_nlogn, [], 1);
    
    
        Q0 = Qfunc(dX, dE, beta0);
    
        beta_Qmin = [cos(theta_nlogn(:,1)).*cos(theta_nlogn(:,2)), cos(theta_nlogn(:,1)).*sin(theta_nlogn(:,2)), sin(theta_nlogn(:,1))];
        beta_B(:,:,b) = [min(beta_Qmin,[],1); 1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
    
    
    end
    
    %num_b_wrong = length(b_wrong);
    
    clearvars -except N D J T B beta0 theta0 theta_B theta_out_B beta_B beta_out_B check_B B_first B_last nlogn core
    eval(sprintf('save Sim_%dD_results_%d_%d_%d',D,B_first,B_last,nlogn))
    
    %% Performance Evaluation
    % for table 2
    % Mid_bias = mean(beta_B(2,:,:) - beta0', 3)
    %
    %
    % Ub_MD = mean(beta_B(3,:,:) - beta0', 3)
    Ub_MSE = mean((beta_B(3,:,:) - beta0').^2, 3);
    %
    % Lb_MD = mean(beta_B(1,:,:) - beta0', 3)
    Lb_MSE = mean((beta_B(1,:,:) - beta0').^2, 3);
    %
    Ub_rsMSE = sqrt(sum(Ub_MSE)) % root of sum of MSE
    Ub_MND = mean(vecnorm(permute(beta_B(3,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    %
    Lb_rsMSE = sqrt(sum(Lb_MSE)) % root of sum of MSE
    Lb_MND = mean(vecnorm(permute(beta_B(1,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    %
    % dUL_mean = mean(beta_B(3,:,:) - beta_B(1,:,:), 3)
    Mid_MSE = mean((beta_B(2,:,:) - beta0').^2, 3);
    % Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    % Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    %
    
    % for table 3
    %Mid_SAB = sum(abs(Mid_bias))
    %dul_sum = sum(dUL_mean)
    Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    % Create filename with nlogn value
    nlogn_str = strrep(num2str(nlogn), '.', '_'); % Replace '.' with '_' for valid filename
    filename = sprintf('idnlogn_%s_result.mat', nlogn_str);
    
    % Save all metrics in one .mat file
    save(filename, 'Ub_rsMSE', 'Ub_MND', 'Lb_rsMSE', 'Lb_MND', 'Mid_rsMSE', 'Mid_MND');
    
    % Optional: Display confirmation
    fprintf('Results saved to %s\n', filename);
    fprintf('Metrics calculated:\n');
    fprintf('Upper bound rsMSE: %.4f, MND: %.4f\n', Ub_rsMSE, Ub_MND);
    fprintf('Lower bound rsMSE: %.4f, MND: %.4f\n', Lb_rsMSE, Lb_MND);
    fprintf('Midpoint rsMSE: %.4f, MND: %.4f\n', Mid_rsMSE, Mid_MND);
    %Mid_MAB = max(abs(Mid_bias),[],2) % max absolute bias
    
    diary off
    
    poolobj = gcp('nocreate');
    evalc('delete(poolobj)');

end
