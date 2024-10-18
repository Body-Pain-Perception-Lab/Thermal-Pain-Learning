%% Run HGF multiple subjects and compare models
clear
%%

addpath('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\HGF_wksh\modified_perceptual_model')
addpath('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\HGF_wksh\Perceptual_model_scripts')
addpath('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\HGF_wksh\Niia\HGFscripts_works\Perceptual model')
addpath('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\HGF_wksh\Perceptual_model_scripts')
addpath('C:\Users\Jespe\OneDrive\Skrivebord\tapas\tapas-master\HGF')
addpath('C:\Users\Jespe\OneDrive\Skrivebord\response models\firstmodel')
%%
load('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\HGF_wksh\Perceptual_model_scripts\data.mat');
%%
% Initialize output - One model for each subject
clear est
% qq = 1:size(data.u,2);
qq = 1:10;

% qq(45)= [];
% qq(44)= [];
% qq(15)= [];

%%
for subject_idx=qq
        
        u(:,1) = data.u(:,subject_idx);
        u(:,2) = data.u_cues(:,subject_idx); %TGI
        des_prob = data.u_densProb(:,subject_idx);
        y = sqrt(data.y_conf(:,3,subject_idx));
        y(:,2) = data.stim1(:,subject_idx);
        y(:,1) = (y(:,1)-nanmean(y(:,1)))/nanstd(y(:,1));
%         y(:,3) = data.predacc1(:,subject_idx);
        hap = data.om2(subject_idx);
        hap(:,2) = data.om3(subject_idx);
        
worked =0;
    
    while ~worked
        try
             est{subject_idx} = tapas_fitModel(y, u, tapas_hgf_binary_pu_tgi_config_test_seq(subject_idx, hap), rm_RT_config_stim_PE_int, tapas_quasinewton_optim_config);
             worked=1;
        end
    end  
    end


%%
p = genpath('C:\Users\Jespe\OneDrive\Skrivebord\HGF_wksh\vba\VBA-toolbox-master')
addpath(p)
