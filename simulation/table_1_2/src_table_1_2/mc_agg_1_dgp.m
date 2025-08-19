%% Section 0: DGP
function [X_B, y_B, EyA_B, theta0] = mc_agg_1_dgp(N, D, J, T, B, beta0, sigma, core)

    rngSeed = 794;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream
       
    beta0norm = beta0 / norm(beta0);
    
    % if D == 3
    theta0 = [asin(beta0norm(2)), atan2(beta0norm(2), beta0norm(1))];
    % elseif D == 4
    % 
    %     theta0 = [asin(beta0norm(3)), atan(beta0norm(3)/sqrt(beta0norm(1)^2 +beta0norm(2)^2) ), atan(beta0norm(2)/beta0norm(1))];
    % 
    % elseif D == 5
    %     theta0 = [ ...
    %         asin(beta0norm(4)), ...
    %         atan2(beta0norm(4), sqrt(beta0norm(1)^2 + beta0norm(2)^2 + beta0norm(3)^2)), ...
    %         atan2(beta0norm(3), sqrt(beta0norm(1)^2 + beta0norm(2)^2)), ...
    %         atan2(beta0norm(2), beta0norm(1))];
    % end

    
    X_B = nan(N,D,J,T,B);
    y_B = nan(N,J,T,B);
    EyA_B = nan(N,J,T,B);
    %AXbeta0_B = [];
    %rng(1);


    parfor b = 1:B
        stream = RandStream('Threefry', 'Seed', rngSeed + b + 200000); % Unique seed per iteration
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


end

%[X_B, y_B, EyA_B, theta0] = mc_agg_1_dgp(N, D, J, T, B, beta0, sigma)