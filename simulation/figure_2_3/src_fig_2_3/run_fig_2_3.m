function run_fig_2_3(varargin)

    %% Instruction of setting number of computational cores
    %
    % Default: 4 cores
    %
    % Run this file with default number of cores: 
    % press "Run Section" in "EDITOR"
    %
    % Change number of cores if your system has more available cores:
    %
    % In matlab command window, run: "un_fig_2_3(number of cores)"
    % Example, put "un_fig_2_3(12)" in command window and press enter
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

    fprintf('Running fig_2_3 with cores=%d\n', core);
    
    %% Initialize parallel pool
    if isempty(gcp('nocreate'))
        parpool(core); 
    end
    %% DGP
    rngSeed = 789;  % Master seed
    globalStream = RandStream('Threefry', 'Seed', rngSeed);
    RandStream.setGlobalStream(globalStream);  % All workers use this stream

    N = 10^7; % # of individuals in the market
    D = 3;  % # of dimension of observable product characteristics (>= 3)
    J = 3;   % # of alternatives including outside option 0 (>= 3)
    T = 2;  % # of time periods (>= 2)
    B_first = 1;
    B = 1;
    
    logname = sprintf("Fig23_N=%d_D=%d_J=%d_T=%d_bigger_index.log",N,D,J,T);
    
    if exist(logname)
        diary off
        delete(logname);
    end
    
    diary(logname)
    
    B_last = B_first + B - 1;
    sigma = J;
    
    
    beta0 = [2; repmat(1,[D-1,1])];
    
    
    %% DGP
    
    beta0norm = beta0 / norm(beta0);
    X = randn(N,D,J,T);
    
    Xbeta0 = permute(sum(X .* repmat(beta0', [N,1,J,T]), 2),[1,3,4,2]);
    
    eps = evrnd(0,1,[N,J,T]); %
    iu = Xbeta0 + eps; % again N*J*T matrix
    
    % calculate y vector by extracting the max of each row and make it 1
    y = double(bsxfun(@eq,iu,max(iu,[],2)));
    
    % calculate true market shares
    expind = exp(Xbeta0);
    sum_expind = sum(expind, 2);
    Ey = expind ./sum_expind;
    
    T0 = T*(T-1)/2;
    
    %% Main part: loop over b
    
    for b = 1 : 1
    
        dEy = [];
        ind_new = [];
        dX = nan(N,D,J,T0);
    
        for j = 1:J
            count = 0;
            for t = 1:T-1
                for s = (t+1):T
                    count = count + 1;
                    ind_new = [ind_new; [1:N]', ones(N,1) * j, ones(N,1) * t, ones(N,1) * s];
                    dEy = [dEy; Ey(:,j,t) - Ey(:,j,s)];
                    dX(:,:,j,count) = X(:,:,j,t) - X(:,:,j,s);
                end
            end
        end
    
    
        %% diagnosis and generate dE objects (N by J by T0)
    
        dEy = permute(reshape(dEy,[N,T0,J]),[1,3,2]);
    
    end
    
    dX_B = cat(4,dX,-dX);
    dEy_B = cat(3,dEy,-dEy);
    
    
    
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
    Extra_Buffer_Polar = 1;
    
    b_wrong = [];
    
    check_B = nan(B,1);
    
    theta_B = nan(3,D-1,B);
    
    theta_out_B = nan(2,D-1,B);
    
    beta_B = nan(3,D,B);
    
    beta_out_B = nan(2,D,B);
    
    
    for b = 1:B
    
        %     fprintf('\nSimulation %d\n', b + B_first - 1)
    
    
        Qmax_test = -1;
        Qmin_test = +Inf;
    
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
            parfor k = 1 : Length_L1
                stream = RandStream('Threefry', 'Seed', rngSeed + k + 7800000); % Unique seed per iteration
                RandStream.setGlobalStream(stream);  
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
    
    
            if Loop_count == 3 && (abs(theta_l_new(2)-range_l) <= pi/8 || abs(theta_u_new(2)-range_u) <= pi/8 )
                Loop_count = 0;
                range_l = 0 ;
                range_u = 2*pi ;
                theta_l = [-0.5*pi, range_l] ;
                theta_u = [0.5*pi, range_u] ;
                range_control = 2;
                continue
            end
    
            %             Qmin
    
    
            if Loop_count >= 3 && Qmin > 0 % if theta size ratio is large, break out of loop
                %                 fprintf('Loop 1 finished at Step %d \r', Loop_count)
    
                %                 Qmin
    
                Loop_control = 2;
                b_min = b_gr(ind_min,:);
    
                break;
    
            elseif Loop_count >= 3 && Qmin == 0
    
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
    
    
            parfor k = 1 : Length_L2
                stream = RandStream('Threefry', 'Seed', rngSeed + k + 8807000); % Unique seed per iteration
                RandStream.setGlobalStream(stream); 
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
    
        %     plot(theta_min(:,1),theta_min(:,2),'m.')
    
    
        theta_Qmin = theta_min;
    
        while Loop_control == 3
    
            Loop_count = Loop_count + 1;
    
            theta_l_out = [max(theta_min_l(1) - (Out_Buffer_Polar + theta_l_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), -0.5 * pi), max(theta_min_l(2) - (Out_Buffer_Polar + theta_l_del_0(2) * Extra_Buffer_Polar) * Tol_polar(2), range_l)];
            theta_u_out = [min(theta_min_u(1) + (Out_Buffer_Polar + theta_u_del_0(1) * Extra_Buffer_Polar) * Tol_polar(1), 0.5 * pi), min(theta_min_u(2) + (Out_Buffer_Polar + theta_u_del_0(2) * Extra_Buffer_Polar)  * Tol_polar(2), range_u)];
    
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
                stream = RandStream('Threefry', 'Seed', rngSeed + k + 9800300); % Unique seed per iteration
                RandStream.setGlobalStream(stream); 
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
    
            if theta_del > 0 && Tol_Step <= Tol
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
    
        %
        plot(theta_Qmin(:,1),theta_Qmin(:,2),'r.')
        hold on
        plot(theta_notmin(:,1),theta_notmin(:,2),'b.')
    
        %bdry = boundary(theta_Qmin(:,1),theta_Qmin(:,2),0);
        posx = min(theta_Qmin(:,1));
        posy = min(theta_Qmin(:,2));
        lenx = max(theta_Qmin(:,1)) - posx;
        leny = max(theta_Qmin(:,2)) - posy;
        rectangle('Position', [posx posy lenx leny],'LineWidth',2,'EdgeColor','k')
        %plot(theta_Qmin(bdry,1),theta_Qmin(bdry,2),'k','LineWidth',2)
    
        plot(theta0(1),theta0(2),'k.','MarkerSize',30)
    
        %  eval(sprintf('title(''theta: b = %d, Q_{min} = %d, Q_0 = %d'')', b+B_first-1, Qmin, Q0))
        %eval(sprintf('saveas(gcf,''Sim_3D_Grid_sieve_theta_%d.png'')', b+B_first-1))
        saveas(gcf, '../result_fig_2_3/figure2.png')
        saveas(gcf, '../result_simulation/figure2.png')
        hold off
        clf
    
        theta_B(:,:,b) = [theta_min_l; mean(theta_Qmin,1);theta_min_u];
        theta_out_B(:,:,b) = [theta_notmin_l; theta_notmin_u];
        %
    
        beta_Qmin = [cos(theta_Qmin(:,1)).*cos(theta_Qmin(:,2)), cos(theta_Qmin(:,1)).*sin(theta_Qmin(:,2)), sin(theta_Qmin(:,1))];
        b_notmin = [cos(theta_notmin(:,1)).*cos(theta_notmin(:,2)), cos(theta_notmin(:,1)).*sin(theta_notmin(:,2)), sin(theta_notmin(:,1))];
        %
        beta_B(:,:,b) = [min(beta_Qmin,[],1); 1/2 * (min(beta_Qmin,[],1) + max(beta_Qmin,[],1)); max(beta_Qmin,[],1)];
        %
        % b_inbound = [min(b_Qmin,[],1); max(b_Qmin,[],1)];
        beta_out_B(:,:,b) = [min(b_notmin,[],1); max(b_notmin,[],1)];
    
        plot3(beta_Qmin(:,1),beta_Qmin(:,2),beta_Qmin(:,3),'r.')
        %axis([0.8, 0.82, 0.4, 0.42, 0.4, 0.42]);
        hold on
        plot3(b_notmin(:,1),b_notmin(:,2),b_notmin(:,3),'b.')
        bdry_b = boundary(b_notmin(:,1),b_notmin(:,2),b_notmin(:,3),0);
        hold on
        %plot3(b_notmin(bdry_b,1),b_notmin(bdry_b,2),b_notmin(bdry_b,3))
        plot3(beta0(1),beta0(2),beta0(3),'k.','MarkerSize',30)
        hold on
        points = beta_Qmin
        % Convert back to spherical to find bounding box
        [az, el, ~] = cart2sph(points(:,1), points(:,2), points(:,3));
    
        % Bounding box in spherical coordinates
        az_min = min(az); az_max = max(az);
        el_min = min(el); el_max = max(el);
    
        % Define the 4 corners in spherical coordinates
        corners_sph = [az_min, el_min;
            az_max, el_min;
            az_max, el_max;
            az_min, el_max];
    
        % Convert rectangle corners back to Cartesian
        [xr, yr, zr] = sph2cart(corners_sph(:,1), corners_sph(:,2), 1);
        plot3([xr; xr(1)], [yr; yr(1)], [zr; zr(1)], 'k-', 'LineWidth', 2);
    
        % eval(sprintf('title(''beta: b = %d, Q_{min} = %d, Q_0 = %d'')', b+B_first-1, Qmin, Q0))
        %eval(sprintf('saveas(gcf,''Sim_3D_trues_beta_%d.png'')', b+B_first-1))
        saveas(gcf, '../result_fig_2_3/figure3.png')
        saveas(gcf, '../../result_simulation/figure3.png')
        hold off
        clf
    
    
    end
    
    num_b_wrong = length(b_wrong);
    
    clearvars -except N D J T B beta0 theta0 theta_B theta_out_B beta_B beta_out_B B_first B_last theta_Qmin theta_notmin beta_Qmin b_notmin
    eval(sprintf('save Sim_%dD_results_%d_%d',D,B_first,B_last))
    
    diary off
end