function MC_AGG_20230909(varargin)

    %----------------------------------------------------------------------
    % File: MC_AGG.m
    %for table 8
    
    % Description: Generate simulation data
    
    % Coded by: Yuan Gao & Ming Li
    % Last modified: 12/01/2018
    rngSeed = 797;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream
   
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
    
        N = 10000; % # of individuals in the market
        D = 4;  % # of dimension of observable product characteristics (>= 3)
        J = 3;   % # of alternatives including outside option 0 (>= 3)
        T = 2;  % # of time periods (>= 2)
        B_first = 1;
        B = 1000;
    
        core = 4;
    end
    
    evalc('parpool(core)');
    
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
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 300000); % Unique seed per iteration
        RandStream.setGlobalStream(stream); 

        Z_i = randn(N,1); % latent variable to induce correlation between X_i2jt and A_i2
    
        X = nan(N,D,J,T);
    
        X(:, 1, :, :) = 2 * rand(N, J, T) - 1 ;
        X(:, 2, :, :) = sqrt(2*sigma) * randn(N,J,T) + repmat(Z_i, 1, J, T);
        X(:, 3:D, :, :) = randn(N, D-2, J, T);
    
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
    
    clearvars -except N D J T B B_first B_last core beta0 theta0 X_B y_B EyA_B rngSeed
    eval(sprintf('save Sim_1data_%d_%d_%d',D,J,T))
    
    
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
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 310000); % Unique seed per iteration
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
            % opts = statset('UseParallel', false, 'Streams', RandStream('twister', 'Seed', b));
            % [L,FitInfo] = lasso(X_j_sieve, Y_j, 'CV', 10, 'Options', opts);
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
    
    dXEname = sprintf('dXdE_%d_%d_%d',D,J,T);
    %save(dXEname,'dX_B','dEy_B');
    clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B rngSeed
    %delete *.csv
    
    %% Sim_3D_parGrid
    
    beta0 = [2; ones(D-1,1)];
    beta0 = beta0 /norm(beta0);
    
    if D == 3
    
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
            stream = RandStream('Threefry', 'Seed', rngSeed + b + 320000); % Unique seed per iteration
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
            %
            beta_Qmin = [cos(theta_Qmin(:,1)).*cos(theta_Qmin(:,2)), cos(theta_Qmin(:,1)).*sin(theta_Qmin(:,2)), sin(theta_Qmin(:,1))];
            %b_notmin = [cos(theta_notmin(:,1)).*cos(theta_notmin(:,2)), cos(theta_notmin(:,1)).*sin(theta_notmin(:,2)), sin(theta_notmin(:,1))];
            %
            %
            beta_B(:,:,b) = [min(beta_Qmin,[],1); 1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
            %
    
    
        end
    
        clearvars -except N D J T B beta0 theta0 theta_B theta_out_B beta_B beta_out_B check_B B_first B_last core rngSeed
        eval(sprintf('save Sim_%d_%d_%d_results',D,J,T))
    
    else
    
        theta0 = zeros((D-1),1);
        tt = 1;
        for i = 1:(D-2)
            theta0(i) = asin(beta0(D-i+1)/tt);
            tt = tt * cos(theta0(i));
        end
        if beta0(1)>=0
            theta0(D-1) = atan(beta0(2)/beta0(1));
        elseif beta0(1)<0 && beta0(2)>=0
            theta0(D-1) = atan(beta0(2)/beta0(1))+pi;
        else
            theta0(D-1) = atan(beta0(2)/beta0(1))-pi;
        end
    
    
        M_Step = 20; % Number of points along each dimension
        Tol = 1e-2;
        Out_Buffer_Polar = 2; % Select a positive number
    
        b_wrong = [];
    
        check_B = nan(B,1);
    
        theta_B = nan(3,D-1,B);
    
        theta_out_B = nan(2,D-1,B);
    
        beta_B = nan(3,D,B);
    
        beta_out_B = nan(2,D,B);
    
    
        parfor b = 1:B
            stream = RandStream('Threefry', 'Seed', rngSeed + b + 330000); % Unique seed per iteration
            RandStream.setGlobalStream(stream); 
    
            Qmax = -1;
            Qmin = +Inf;
    
            dX = dX_B(:,:,:,:,b);
            dE = dEy_B(:,:,:,b);
    
            range_control = 1 ;
            %determine the range: 1 means (-pi, pi], 2 means (0, 2*pi]
    
            range_l = -pi ;
            range_u = pi ;
            theta_l = [-0.5 * pi * ones(1, D-2), range_l] ;
            theta_u = [0.5 * pi * ones(1, D-2), range_u] ;
    
    
            Loop_count = 0;
            Loop_control = 1;
            Tol_polar = 2* pi;
            theta_scale_ratio = 0;
    
            Extra_Buffer_Polar = 2;
    
            %     tic
    
            while (Loop_control == 1)
    
                Loop_count = Loop_count + 1;
    
                Tol_polar = (theta_u - theta_l) ./ [(M_Step - 1) * ones(1, D-2), M_Step];
    
                theta = theta_l + [0:1:(M_Step-1)]' * Tol_polar;
    
                %     theta = nan(M_Step, D-1);
                %     theta(:, [1:(D-2)]) = ones(M_Step,1) *  theta_l([1:(D-2)]) +  [0:1:(M_Step-1)]' / (M_Step-1)  * (theta_u([1:(D-2)]) - theta_l([1:(D-2)])) ;
                %     theta(:, D-1) = theta_l(D-1) +  [0:1:(M_Step-1)]' / M_Step  * (theta_u(D-1) - theta_l(D-1)) ;
                %
                theta_cell = cell(1,D-1);
                for d = 1:(D-1)
                    theta_cell{d} = theta(:, d);
                end
                Nd_grid = cell(1, D-1);  % create cell array to receive outputs of ndgrid
                [Nd_grid{:}] = ndgrid(theta_cell{:});  % fills all N Nd arrays
    
                Length_L1 = M_Step^(D-1);
                theta_gr = nan(Length_L1, D-1);
                for d = 1:(D-1)
                    theta_gr(:, d) = reshape(Nd_grid{d},[],1);
                end
    
                b_gr = nan(Length_L1, D);
    
                for d = 1:D
    
                    temp = ones(Length_L1,1);
                    for d_temp = 1:(d-1)
                        if d == 1
                            temp = temp;
                        elseif d > 1
                            temp = temp .*  cos(theta_gr(:, d_temp));
                        end
                    end
                    if d < D
                        temp = temp .* sin(theta_gr(:, d));
                    end
                    b_gr(:,D-d+1) = temp;
                    %b_gr = [cos(theta1_gr).*cos(theta2_gr).*cos(theta3_gr), cos(theta1_gr).*cos(theta2_gr) .* sin(theta3_gr), cos(theta1_gr).*sin(theta2_gr), sin(theta1_gr)];
                end
    
                %  clear theta theta1_gr theta2_gr Theta1 Theta2
    
                Q_gr = nan(Length_L1,1);
                for k = 1 : Length_L1
                    Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
                end
    
                Qmax_new = max(Q_gr);
                Qmin_new = min(Q_gr);
    
                Qmax = max(Qmax, Qmax_new); % don't need max() by construction?
                Qmin = min(Qmin, Qmin_new);
    
                ind_qt = find(Q_gr <= quantile(Q_gr, 1 / (2^(D-1))));
    
                theta_gr_new = theta_gr(ind_qt, :);
                theta_l_new = min(theta_gr_new, [], 1);
                theta_u_new = max(theta_gr_new, [], 1);
    
                theta_del = max(max(abs(theta_u_new - theta_u)), max(abs(theta_l_new - theta_l)));
    
                ind_min= find(Q_gr == Qmin_new);
                theta_min = theta_gr(ind_min,:);
    
                ind_notmin = find(Q_gr ~= Qmin_new);
                theta_notmin = theta_gr(ind_notmin,:);
    
                theta_min_l = min(theta_min, [], 1);
                theta_min_u = max(theta_min, [], 1);% these are boundary of points that correspond to min of Q
    
                theta_notmin_l = min(theta_notmin, [], 1);
                theta_notmin_u = max(theta_notmin, [], 1);
    
                theta_scale = max(theta_min_u - theta_min_l);
                theta_scale_ratio = theta_scale / max(theta_u - theta_l);
    
                theta_l_del = theta_min_l - theta_notmin_l;
                theta_u_del = theta_notmin_u - theta_min_u;
                theta_del = min(min(theta_l_del), min(theta_u_del));
    
                theta_l_del_0 = (theta_min_l - theta_notmin_l == 0);
                theta_u_del_0 = (theta_notmin_u - theta_min_u == 0);
    
    
                if range_control == 1 && Loop_count == 3 && (abs(theta_l_new(D-1)-range_l) < pi/8 || abs(theta_u_new(D-1)-range_u) < pi/8)
                    Loop_count = 0;
                    range_l = 0 ;
                    range_u = 2*pi ;
                    theta_l = [-0.5 * pi * ones(1, D-2), range_l] ;
                    theta_u = [0.5 * pi * ones(1, D-2), range_u] ;
                    range_control = 2;
                    continue
                end
    
                %&& theta_scale_ratio >= 0.1
    
                if Loop_count >= 3  && Qmin > 0 % if theta size ratio is large, break out of loop
                    %                 fprintf('Loop 1 finished at Step %d \r', Loop_count)
                    %                 Qmin
                    Loop_control = 2;
                    b_min = b_gr(ind_min,:);
    
                    break;
                elseif Loop_count >= 3 && Qmin == 0
    
                    %                 fprintf('Loop 1 finished at Step %d \r', Loop_count)
    
                    Loop_control = 3;
                    break;
                end
    
                theta_l =  theta_l_new;
                theta_u =  theta_u_new;
            end
    
            ind_notmin = find(Q_gr ~= Qmin);
            theta_notmin = theta_gr(ind_notmin,:);
            theta_notmin_l = min(theta_notmin, [], 1);
            theta_notmin_u = max(theta_notmin, [], 1) ;
    
            theta_l_del = theta_min_l - theta_notmin_l;
            theta_u_del = theta_notmin_u - theta_min_u;
            theta_del = min(min(theta_l_del), min(theta_u_del));
    
            theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
            theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
            Tol_Step = min(Tol_polar) / 5; % this tol setup
    
            theta_u_out = nan(1,D-1);
            theta_l_out = nan(1,D-1);
            while Loop_control == 2
    
                Loop_count = Loop_count + 1; %initialize loop count for loopcontrol 2?
    
                for d = 1:(D-2)
                    theta_l_out(d) = max(theta_min_l(d) - (Out_Buffer_Polar + theta_l_del_0(d) * Extra_Buffer_Polar) * Tol_polar(d), -0.5 * pi);
                    theta_u_out(d) = min(theta_min_u(d) + (Out_Buffer_Polar + theta_u_del_0(d) * Extra_Buffer_Polar) * Tol_polar(d), 0.5 * pi);
                end
    
                theta_l_out(D-1) = max(theta_min_l(D-1) - (Out_Buffer_Polar + theta_l_del_0(D-1) * Extra_Buffer_Polar)  * Tol_polar(D-1), range_l);
                theta_u_out(D-1) = min(theta_min_u(D-1) + (Out_Buffer_Polar + theta_u_del_0(D-1) * Extra_Buffer_Polar) * Tol_polar(D-1), range_u);
    
                theta_l_in = theta_min_l;
                theta_u_in = theta_min_u;
    
    
                %     M_Step_L2 = floor(max((theta_u_out - theta_l_out) / Tol_Step)) + 1; %adaptive number of parts cut
                %     Tol_polar = [(theta_u_out(1) - theta_l_out(1)) / (M_Step_L2 - 1), (theta_u_out(2) - theta_l_out(2)) / M_Step_L2]
                %
                %
                %     theta = theta_l_out +  [0:1:(M_Step_L2-1)]' * Tol_polar ;
                %
                %     [Theta1, Theta2] = ndgrid(theta(:,1), theta(:,2));
    
                M_Step_L2 = floor((theta_u_out - theta_l_out) ./ Tol_polar) + 1;
                %Tol_polar = (theta_u_out - theta_l_out) ./ [(M_Step_L3 - 1) * ones(1, D-2), M_Step_L3];
    
                %theta = theta_l_out +  [0:1:(M_Step_L3-1)]' * Tol_polar ;
    
                theta_cell = cell(1,D-1);
                
                for d = 1:(D-1)
                    theta_cell{d} = theta_l_out(d) +  [0:1:(M_Step_L2(d)-1)]' * Tol_polar(d);
                end
    
    
                Nd_grid = cell(1, D-1);  % create cell array to receive outputs of ndgrid
                [Nd_grid{:}] = ndgrid(theta_cell{:});  % fills all N Nd arrays
    
                Length_L2 = prod(M_Step_L2);
                theta_gr = nan(Length_L2, D-1);
                for d = 1:(D-1)
                    theta_gr(:, d) = reshape(Nd_grid{d},[],1);
                end
    
    
                b_gr = nan(Length_L2, D);
    
    
                for d = 1:D
    
                    temp = ones(Length_L2,1);
                    for d_temp = 1:(d-1)
                        if d == 1
                            temp = temp;
                        elseif d > 1
                            temp = temp .*  cos(theta_gr(:, d_temp));
                        end
                    end
                    if d < D
                        temp = temp .* sin(theta_gr(:, d));
                    end
                    b_gr(:,D-d+1) = temp;
                    %b_gr = [cos(theta1_gr).*cos(theta2_gr).*cos(theta3_gr), cos(theta1_gr).*cos(theta2_gr) .* sin(theta3_gr), cos(theta1_gr).*sin(theta2_gr), sin(theta1_gr)];
                end
    
    
                %     clear theta_gr
    
                Q_gr = nan(Length_L2,1);
                for k = 1 : Length_L2
                    Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
                end
    
    
                % clear theta theta1_gr theta2_gr Theta1 Theta2 % ind_temp1 ind_temp2 ind_temp
    
                %     b_gr = [cos(theta_gr(:,1)).*cos(theta_gr(:,2)), cos(theta_gr(:,1)).*sin(theta_gr(:,2)), sin(theta_gr(:,1))];
                %
                %     Length_L2 = size(theta_gr,1);
                %     Q_gr = nan(Length_L2,1);
                %
                %
                %     for k = 1 : Length_L2
                %         Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
                %     end
    
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
    
                if Qmin_new >= Qmin || Qmin_new == 0 % if next round Qmin of finer points is no smaller than current Qmin, enter final loop
    
                    Qmin = min(Qmin, Qmin_new);
                    %             fprintf('Loop 2 finished at Step %d \r', Loop_count)
                    Tol_Step = Tol_Step / 2;
                    Loop_control = 3;
                    break;
    
                else
    
                    Qmin = Qmin_new;
    
                end
    
            end
    
    
    
    
            theta_Qmin = theta_min;
    
            %Tol_Step = min(Tol_polar) / 2;
    
    
            while (Loop_control == 3)
    
                Loop_count = Loop_count + 1;
    
                for d = 1:(D-2)
                    theta_l_out(d) = max(theta_min_l(d) - (Out_Buffer_Polar + theta_l_del_0(d) * Extra_Buffer_Polar) * Tol_polar(d), -0.5 * pi);
                    theta_u_out(d) = min(theta_min_u(d) + (Out_Buffer_Polar + theta_u_del_0(d) * Extra_Buffer_Polar) * Tol_polar(d), 0.5 * pi);
                end
    
                theta_l_out(D-1) = max(theta_min_l(D-1) - (Out_Buffer_Polar + theta_l_del_0(D-1) * Extra_Buffer_Polar)  * Tol_polar(D-1), range_l);
                theta_u_out(D-1) = min(theta_min_u(D-1) + (Out_Buffer_Polar + theta_u_del_0(D-1) * Extra_Buffer_Polar) * Tol_polar(D-1), range_u);
    
                theta_l_in = theta_min_l;
                theta_u_in = theta_min_u;
    
    
                %%%%%%%%%%
                M_Step_L3 = floor((theta_u_out - theta_l_out) ./ Tol_polar) + 1;
                %Tol_polar = (theta_u_out - theta_l_out) ./ [(M_Step_L3 - 1) * ones(1, D-2), M_Step_L3];
    
                %theta = theta_l_out +  [0:1:(M_Step_L3-1)]' * Tol_polar ;
    
                theta_cell = cell(1,D-1);
                for d = 1:(D-1)
                    theta_cell{d} = theta_l_out(d) +  [0:1:(M_Step_L3(d)-1)]' * Tol_polar(d);
                end
                Nd_grid = cell(1, D-1);  % create cell array to receive outputs of ndgrid
                [Nd_grid{:}] = ndgrid(theta_cell{:});  % fills all N Nd arrays
    
                Length_L3 = prod(M_Step_L3);
                theta_gr = nan(Length_L3, D-1);
                for d = 1:(D-1)
                    theta_gr(:, d) = reshape(Nd_grid{d},[],1);
                end
    
    
                %%%%%%%%%%%%%
    
    
    
                ind_temp1 = (min((theta_gr >= theta_l_in), [], 2) == 1);
                ind_temp2 = (min((theta_gr <= theta_u_in), [], 2) == 1);
                ind_temp = min(ind_temp1, ind_temp2);
                ind_keep = find(ind_temp ~= 1);
    
                Length_L3 = length(ind_keep);
                theta_gr = theta_gr(ind_keep,:);
    
                b_gr = nan(Length_L3, D);
    
    
                for d = 1:D
    
                    temp = ones(Length_L3,1);
                    for d_temp = 1:(d-1)
                        if d == 1
                            temp = temp;
                        elseif d > 1
                            temp = temp .*  cos(theta_gr(:, d_temp));
                        end
                    end
                    if d < D
                        temp = temp .* sin(theta_gr(:, d));
                    end
                    b_gr(:,D-d+1) = temp;
                    %b_gr = [cos(theta1_gr).*cos(theta2_gr).*cos(theta3_gr), cos(theta1_gr).*cos(theta2_gr) .* sin(theta3_gr), cos(theta1_gr).*sin(theta2_gr), sin(theta1_gr)];
                end
    
    
    
    
                Q_gr = nan(Length_L3,1);
                for k = 1 : Length_L3
                    Q_gr(k) = Qfunc(dX, dE, b_gr(k,:)'); % use Qfunc in this same folder!!
                end
    
                Qmax_new = max(Q_gr);
                Qmin_new = min(Q_gr);
    
                %         if Qmin_new < Qmin
                %             warning('smaller Q encountered')
                %         end
    
                %     Qmax = max(Qmax, Qmax_new);
                %     Qmin = min(Qmin, Qmin_new);
    
                ind_min = find(Q_gr == Qmin);
                ind_notmin = find(Q_gr ~= Qmin);
    
                theta_min_new = theta_gr(ind_min,:);
                theta_min_l_new= min(theta_min_new, [], 1);
                theta_min_u_new = max(theta_min_new, [], 1);
    
                theta_Qmin = [theta_Qmin; theta_min_new];
                theta_min_l = min(theta_Qmin, [], 1);
                theta_min_u = max(theta_Qmin, [], 1);
    
                theta_notmin = theta_gr(ind_notmin,:);
                theta_notmin_l = min(theta_notmin, [], 1);
                theta_notmin_u = max(theta_notmin, [], 1);
    
                theta_l_del = theta_min_l - theta_notmin_l;
                theta_u_del = theta_notmin_u - theta_min_u;
                theta_del = min(min(theta_l_del), min(theta_u_del));
    
                theta_l_del_0 = (theta_min_l - theta_notmin_l <= 0);
                theta_u_del_0 = (theta_notmin_u - theta_min_u <= 0);
    
                %   clear theta_gr
    
                if theta_del > 0 && Tol_Step <= Tol
                    Loop_control = 4;
                    Tol_Final = Tol_Step;
                    break;
                end
    
                if theta_del > 0
                    Extra_Buffer_Polar = 2;
                    Tol_polar  = max(Tol_polar /2, Tol);
                    Tol_Step = max(Tol_polar);
                else
                    Extra_Buffer_Polar = Extra_Buffer_Polar  + 1;
                end
    
            end
    
    
            Q0 = Qfunc(dX, dE, beta0);
            theta_B(:,:,b) = [theta_min_l; mean(theta_Qmin,1);theta_min_u];
            beta_Qmin = [cos(theta_Qmin(:,1)).*cos(theta_Qmin(:,2)).*cos(theta_Qmin(:,3)), cos(theta_Qmin(:,1)).*cos(theta_Qmin(:,2)).*sin(theta_Qmin(:,3)), cos(theta_Qmin(:,1)).*sin(theta_Qmin(:,2)),sin(theta_Qmin(:,1))];
            beta_B(:,:,b) = [min(beta_Qmin,[],1);1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
    
        end
    end
    
    clearvars -except N D J T B beta0 theta0 theta_B theta_out_B beta_B beta_out_B check_B B_first B_last core rngSeed
    eval(sprintf('save Sim_%d_%d_%d_results',D,J,T))
    
    %num_b_wrong = length(b_wrong);
    
    
    %% Performance Evaluation
    % for table 2
    %Mid_bias = mean(beta_B(2,:,:) - beta0', 3)
    
    
    %Ub_MD = mean(beta_B(3,:,:) - beta0', 3)
    %Ub_MSE = mean((beta_B(3,:,:) - beta0').^2, 3)
    
    %Lb_MD = mean(beta_B(1,:,:) - beta0', 3)
    %Lb_MSE = mean((beta_B(1,:,:) - beta0').^2, 3)
    
    %Ub_rsMSE = sqrt(sum(Ub_MSE)) % root of sum of MSE
    %Ub_MND = mean(vecnorm(permute(beta_B(3,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    %Lb_rsMSE = sqrt(sum(Lb_MSE)) % root of sum of MSE
    %Lb_MND = mean(vecnorm(permute(beta_B(1,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    %dUL_mean = mean(beta_B(3,:,:) - beta_B(1,:,:), 3)
    Mid_MSE = mean((beta_B(2,:,:) - beta0').^2, 3); %#ok<SYNER>
    %Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    %Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    
    % for table 3
    %Mid_SAB = sum(abs(Mid_bias))
    %dul_sum = sum(dUL_mean)
    Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    % Create filename using only D, J, T
    filename = sprintf('results_D%d_J%d_T%d.mat', D, J, T);
    
    % Save all metrics and parameters
    save(filename, 'Mid_rsMSE', 'Mid_MND');
    %Mid_MAB = max(abs(Mid_bias),[],2) % max absolute bias
    
    diary off
    
    poolobj = gcp('nocreate');
    evalc('delete(poolobj)');
    clearvars -except core
end