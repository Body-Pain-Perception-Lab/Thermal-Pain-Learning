%model comparison
function [posterior,out] = model_compar(results)

for i=1:size(results,1) %loop sbjects
    for j=1:size(results,2) %loop models
   
        if isfield(results{i,j},'optim')
            lmes(i,j) = results{i,j}.optim.LME;
            bic (i,j) = results{i,j}.optim.BIC;
        else
            lmes(i,j) = NaN;
            bic(i,j) = NaN;
        end
    end
end

% lmes(15,:) = [];
lmes1 = real(real(lmes)')
[posterior,out] = VBA_groupBMC(lmes1)
end