function [est] = fit_hgf_subject_model_final_with_tgi1(data, perceptual_config, response_config, optim_config )
% fit_hgf
% Loop through all subjects

%% setup

est = cell (size(data.u,2),length(perceptual_config));

for model_idx =1:length(perceptual_config)
    
    
    for subject_idx=1:size(data.u,2)
        
        u(:,1) = data.u(:,subject_idx);
        u(:,2) = data.u_cues(:,subject_idx); %TGI
        des_prob = data.u_densProb(:,subject_idx);
        
        y = data.y(:,subject_idx);
        
%         if (data.v(subject_idx) == 0)
%             continue
%         end
               
        
        %% fit parameter values for subject
        try
            est{subject_idx,model_idx} = tapas_fitModel(y, u, perceptual_config{model_idx}, response_config, optim_config);
            
        catch
            est{subject_idx,model_idx} = NaN;
        end
        
    display(subject_idx)    
    end
end


end