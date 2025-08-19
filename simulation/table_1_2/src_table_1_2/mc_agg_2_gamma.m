%% Section 2: gamma
function [dX_B, dEy_B] = mc_agg_2_gamma(X_B, y_B, EyA_B, N, D, J, T, B, B_first, B_last)
    
    rngSeed = 794;  % Master seed
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
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 210000); % Unique seed per iteration
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
    
    dXEname = sprintf('dXdE_%d_%d_N=%d',B_first,B_last,N);
    save(dXEname,'dX_B','dEy_B');

    %clearvars -except N D J T B B_first B_last core beta0 theta0 dX_B dEy_B
    %delete *.csv
end