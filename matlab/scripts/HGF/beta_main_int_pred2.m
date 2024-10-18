%get the beta values function
function [beta] = beta_main_int_pred2(results)

for i=1:size(results,1) %loop sbjects
    for j=1 %loop models
        if isfield(results{i,j},'p_prc')
            beta(i,j) = results{i,j}.p_obs.be0;
            beta(i,j+1) = results{i,j}.p_obs.be1;
            beta(i,j+2) = results{i,j}.p_obs.be2;
            beta(i,j+3) = results{i,j}.p_obs.be3;            
            beta(i,j+4) = results{i,j}.p_obs.be4;
            beta(i,j+5) = results{i,j}.p_obs.be5;
            beta(i,j+6) = results{i,j}.p_obs.be6;            
            beta(i,j+7) = results{i,j}.p_obs.be7;
            beta(i,j+8) = results{i,j}.p_obs.be8;     
%             beta(i,j+9) = results{i,j}.p_obs.be9;
%             beta(i,j+10) = results{i, j}.p_obs.sigma_rt;
        else
            display('heelo');
        end
    end
    
end

boxplot(beta)
title('betas')
end
