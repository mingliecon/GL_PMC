%% part 4: Evaluation

function [Mid_bias, Ub_MD, Lb_MD, dUL_mean, Mid_MSE, b_mean, SD, RMSE, Mid_rsMSE, Mid_MND, Mid_SAB, dul_sum] = mc_agg_4_evaluation(beta_B, beta0)
       % for table 2
   
    Mid_bias = mean(beta_B(2,:,:) - beta0', 3)
    
    
    Ub_MD = mean(beta_B(3,:,:) - beta0', 3)
    %Ub_MSE = mean((beta_B(3,:,:) - beta0').^2, 3)
    
    Lb_MD = mean(beta_B(1,:,:) - beta0', 3)
    %Lb_MSE = mean((beta_B(1,:,:) - beta0').^2, 3)
    
    %Ub_rsMSE = sqrt(sum(Ub_MSE)) % root of sum of MSE
    %Ub_MND = mean(vecnorm(permute(beta_B(3,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    %Lb_rsMSE = sqrt(sum(Lb_MSE)) % root of sum of MSE
    %Lb_MND = mean(vecnorm(permute(beta_B(1,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    dUL_mean = mean(beta_B(3,:,:) - beta_B(1,:,:), 3)
    Mid_MSE = mean((beta_B(2,:,:) - beta0').^2, 3);
    b_mean = mean(beta_B(1,:,:), 3);
    SD = mean((beta_B(2,:,:) - b_mean).^2, 3).^ 0.5
    RMSE = mean((beta_B(2,:,:) - beta0').^2, 3).^ 0.5
    
    Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    
    % for table 3
    Mid_SAB = sum(abs(Mid_bias))
    dul_sum = sum(dUL_mean)
    Mid_rsMSE = sqrt(sum(Mid_MSE)) % root of sum of MSE
    Mid_MND = mean(vecnorm(permute(beta_B(2,:,:) - beta0', [3 2 1]), 2,2),1) % mean normed deviation
    
    
    %Mid_MAB = max(abs(Mid_bias),[],2) % max absolute bias
    
    %diary off
    
    %poolobj = gcp('nocreate');
    %evalc('delete(poolobj)');

    
end

