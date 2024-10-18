


p_prc = [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2.5 -6, log(.5), 0 , 1];
% b0, b1
p_rsp = [0, 0, 5];


config_prc = tapas_ehgf_binary_pu_config;

% 1. simulate data
sim = tapas_simModel(u,...
                     'tapas_ehgf_binary_pu',...
                     p_prc,...
                     'rm_RT',...
                     p_rsp);
                 
                 
              

%%
dataFile = load('HGFwksh_CWT_multsub_noise.mat');
u = dataFile.u(:,1);
y = dataFile.y_rts(:,1);
y = log(y);
%%
% 2. fit a single subject
%first model

addpath('C:\Users\Jespe\OneDrive\Skrivebord\response models\firstmodel')


fit = tapas_fitModel(y,u,...
                     'tapas_ehgf_binary_pu_config',...   %now we're giving in priors
                     'rm_RT_config');
                 
                 
sim_post = tapas_simModel(u,...
                     'tapas_ehgf_binary_pu',...
                     fit.p_prc.p,...
                     'rm_RT',...
                     fit.p_obs.p);                 
                 
histogram(y,10) ; hold on;
histogram(sim_post.y,10);


fit.p_obs.be1



%%
%surprise model:
addpath('C:\Users\Jespe\OneDrive\Skrivebord\response models\surprisemodel')

%fitting with another parameter
fit = tapas_fitModel(y,u,...
                     'tapas_ehgf_binary_pu_config',...   %now we're giving in priors
                     'rm_RT_2_config');
                 
                 
sim_post = tapas_simModel(u,...
                     'tapas_ehgf_binary_pu',...
                     fit.p_prc.p,...
                     'rm_RT_2',...
                     fit.p_obs.p);                 
                 
histogram(y,10) ; hold on;
histogram(sim_post.y,10);


fit.p_obs.be1



%%
%try is with new priors sd on them are (0.1) and we ahve 3 beta estiamtes

addpath('C:\Users\Jespe\OneDrive\Skrivebord\response models\trymodel')
%fitting with another parameter
fit = tapas_fitModel(y,u,...
                     'tapas_ehgf_binary_pu_config',...   %now we're giving in priors
                     'rm_RT_2_config');
                 
                 
sim_post = tapas_simModel(u,...
                     'tapas_ehgf_binary_pu',...
                     fit.p_prc.p,...
                     'rm_RT_2',...
                     fit.p_obs.p);                 
                 
histogram(y,10) ; hold on;
histogram(sim_post.y,10);


fit.p_obs.be1