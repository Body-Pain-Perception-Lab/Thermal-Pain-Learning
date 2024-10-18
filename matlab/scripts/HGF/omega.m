%get the beta values function
function [omega] = omega(results)

for i=1:size(results,1) %loop sbjects
    for j=1 %loop models
        if isfield(results{i,j},'p_prc')
            omega(i,j) = results{i,j}.p_obs.tau;
            omega(i,j+1) = results{i,j}.p_prc.om(2);
            omega(i,j+2) = results{i,j}.p_prc.om(3);
        else
            display('heelo');
        end
    end
    
end

boxplot(omega)
title('omegas')
end