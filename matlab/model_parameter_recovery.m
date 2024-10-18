
rng(123)

% Get the directory of the main script
mainScriptPath = fileparts(which('model_parameter_recovery'));

disp(mainScriptPath);

% Construct the path to the 'scripts' folder
scriptsFolderPath = fullfile(mainScriptPath, 'scripts');

tapas = fullfile(mainScriptPath, 'tapas');

% Add the folder and all its subfolders to the MATLAB path
tapas_path = genpath(tapas);

% Generate a path string with all subfolders under 'scripts'
allSubfoldersPath = genpath(scriptsFolderPath);

% Add all subfolders to the MATLAB path
addpath(allSubfoldersPath);
addpath(tapas_path);

%%
matlabDataDir = fullfile(mainScriptPath, 'matlab data');

% Check if the directory exists
if isfolder(matlabDataDir)
    % Get the list of files in the directory
    fileList = dir(matlabDataDir);
    
    % Remove the '.' and '..' entries
    fileList = fileList(~ismember({fileList.name}, {'.', '..'}));
    
    % Check if the directory contains any files
    if isempty(fileList)
        % Run the get_matlab_data function
        get_matlab_data(mainScriptPath);
    else
        disp('The "matlab data" directory is not empty. Either this function as already been run or you should remove the matlab data directory');
    end
else
    % Create the directory and run the function
    mkdir(matlabDataDir);
    disp('The "matlab data" directory has been created.');
    get_matlab_data(mainScriptPath);
end
%%
% data preperation and analysis
%loading the data:
%two different contingencies wether one has a even or uneven stormdb
%number, further two way to incode TGI trials either as 1's and 0's and
%oppisite of the binary choice responses or as 0.5. We'll use the 1's and
%0's

%code to make data from raw files this is where tgi trials are u != y:
%prep_hgf_data_pain_learning_test_tgi10(74)

%if one wants to use the coding 0.5 for the tgi trials use:
% prep_hgf_data_pain_learning_real_y(74)
%this is where the ratings determine the u for tgi-trial
%then load the files that was created

removers = readmatrix('removers2.csv');
removers = removers(1:9,2);
removers = unique(removers,"row")

prep_hgf_data_pain_learning_feel(273,removers,mainScriptPath);

datafeel = load('HGFwksh_TPL_feelmultsub_TGI.mat');
data = datafeel;
%% now that the data is in we can see the two different contingency spaces:
plot(data.u(:,1), '.', 'Color', [0 0.6 0], 'MarkerSize', 11)
xlabel('Trial number')
ylabel('u')
axis([1, 320, -0.1, 1.1])
%% 
plot(data.u(:,2), '.', 'Color', [0 0.6 0], 'MarkerSize', 11)
xlabel('Trial number')
ylabel('u')
axis([1, 320, -0.1, 1.1])
%% Parameter recovery simulations
u(:,1) = data.u(:,1);
u(:,2) = data.u_cues(:,1);

seed = 111;
n_sim = 200;
%hgf

sim_omega2 =  -12+(4+12)*rand(1,n_sim);
sim_omega3 = -12+(4+12)*rand(1,n_sim);
sim_kappa = 0.05+(2-0.05)*rand(1,n_sim);
%rw
sim_alpha = 0.05+(0.9-0.05)*(rand(1,n_sim));
%sutton
sim_mu =  0.1+(10-0.05)*(rand(1,n_sim));
%pearce hall
sim_S = 0.05+(0.9-0.05)*(rand(1,n_sim));
%modified pearce hall
sim_eta = 0.05+(0.9-0.05)*(rand(1,n_sim));
sim_kappa_modph = 0.05+(0.9-0.05)*(rand(1,n_sim));

%response
sim_zeta = 0.5+(15-0.5)*(rand(1,n_sim));

rec2_omega2 = size(sim_omega2);
rec2_zeta_hgf2 = size(sim_omega2);


rec_omega2 = size(sim_omega2);
rec_omega3 = size(sim_omega3);
rec_kappa = size(sim_kappa);
rec_alpha = size(sim_alpha);
rec_mu = size(sim_mu);
rec_S = size(sim_S);
rec_eta = size(sim_eta);
rec_kappa_modph = size(sim_kappa_modph);

rec_zeta_hgf = size(sim_zeta);
rec_zeta_rw = size(sim_zeta);
rec_zeta_su1 = size(sim_zeta);

rec_zeta_ph = size(sim_zeta);
rec_zeta_modph = size(sim_zeta);


for i = 1:n_sim
    p_obs = sim_zeta(i);
    p_prc_hgf = [NaN 0 1 NaN 1 1 NaN 0 0 1 sim_kappa(i) NaN sim_omega2(i) sim_omega3(i)];
    p_prc_rw = [0.5,sim_alpha(i)];
    p_prc_su1 = [sim_mu(i),1,0.5,0.005];
    p_prc_ph = [0.5,0.5,sim_S(i)];
    p_prc_modph = [0.5,0.5,0.5,sim_eta(i),sim_kappa_modph(i)];

    %setting kappa too close to -inf makes an easy way to make a 2-level
    %HGF as all information flow upwards gets disconnected
    p_prc_hgf2 = [NaN 0 1 NaN 1 1 NaN 0 0 1 0 NaN sim_omega2(i) -6];

%2 level hgf
try

    exist('tapas_simModel', 'file')
    exist('tapas_hgf_binary_pu_tgi_test', 'file')
    exist('tapas_unitsq_sgm', 'file')
    sim = tapas_simModel(u,...
                         'tapas_hgf_binary_pu_tgi_test',...
                         p_prc_hgf2,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);

    est = tapas_fitModel(sim.y, sim.u, tapas_hgf_binary_tgi_config_wide_priors_2_level,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    rec2_omega2(i) = est.p_prc.om(2);
    rec2_zeta_hgf2(i) = est.p_obs.ze(1);
catch
    ("Bad HGF simulation");
    rec2_omega2(i) = NaN;
    rec2_zeta_hgf2(i) = NaN;
end
% rescorla wagner
try
    sim = tapas_simModel(u,...
                         'tapas_rw_binary',...
                         p_prc_rw,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);
    
    est = tapas_fitModel(sim.y, sim.u, tapas_rw_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    
    rec_alpha(i) = est.p_prc.al(1);
    rec_zeta_rw(i) = est.p_obs.ze(1);  
catch
    ("Bad RW simulation");
    rec_alpha(i) = NaN;
    rec_zeta_rw(i) = NaN;
end
% pearcehall
try
    sim = tapas_simModel(u,...
                         'tapas_ph_binary',...
                         p_prc_ph,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);
    
    est = tapas_fitModel(sim.y, sim.u, tapas_ph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    
    rec_S(i) = est.p_prc.S(1);
    rec_zeta_ph(i) = est.p_obs.ze(1);  
catch
    ("Bad PH simulation");
    rec_S(i) = NaN;
    rec_zeta_ph(i) = NaN;
end
% mod-pearcehall
try
    sim = tapas_simModel(u,...
                         'tapas_modph_binary',...
                         p_prc_modph,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);
    
    est = tapas_fitModel(sim.y, sim.u, tapas_modph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    
    rec_eta(i) = est.p_prc.eta(1);
    rec_kappa_modph(i) = est.p_prc.kap(1);
    rec_zeta_modph(i) = est.p_obs.ze(1);  
catch
    ("Bad MODPH simulation");
    rec_eta(i) = NaN;
    rec_kappa_modph(i) = NaN;
    rec_zeta_modph(i) = NaN;  
end
% sutton k1
try
    sim = tapas_simModel(u,...
                         'tapas_sutton_k1_binary',...
                         p_prc_su1,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);
    
    est = tapas_fitModel(sim.y, sim.u, tapas_sutton_k1_binary_test_optimal_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    
    rec_mu(i) = est.p_prc.mu(1);
    rec_zeta_su1(i) = est.p_obs.ze(1);

catch
    
    ("Bad SUTTON K1 simulation");
    rec_mu(i) = NaN;
    rec_zeta_su1(i) = NaN;

end
%3 level hgf
try
    sim = tapas_simModel(u,...
                         'tapas_hgf_binary_pu_tgi_test',...
                         p_prc_hgf,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);
    
    est = tapas_fitModel(sim.y, sim.u, tapas_hgf_binary_tgi_config_wide_priors,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    
    rec_omega2(i) = est.p_prc.om(2);
    rec_omega3(i) = est.p_prc.om(3);
    rec_kappa(i) = est.p_prc.ka(2);
    rec_zeta_hgf(i) = est.p_obs.ze(1);
catch
    ("Bad 3-level HGF simulation");
    rec_omega2(i) = NaN;
    rec_omega3(i) = NaN;
    rec_kappa(i) = NaN;
    rec_zeta_hgf(i) = NaN;
end
end
% saving the results
hgf_pr2 = [sim_omega2; rec2_omega2; sim_zeta; rec2_zeta_hgf2];
hgf_pr = [sim_omega2; rec_omega2; sim_omega3; rec_omega3; sim_kappa; rec_kappa; sim_zeta; rec_zeta_hgf];
rw_pr = [sim_alpha; rec_alpha; sim_zeta; rec_zeta_rw];
rw_su1 = [sim_mu; rec_mu; sim_zeta; rec_zeta_su1];

ph = [sim_S; rec_S; sim_zeta; rec_zeta_ph];
modph = [sim_kappa_modph; rec_kappa_modph; sim_eta; rec_eta; sim_zeta; rec_zeta_modph];

csvwrite("created files/hgf2_parameter_recovery.csv",hgf_pr2)
csvwrite("created files/hgf_parameter_recovery.csv",hgf_pr)
csvwrite("created files/rw_parameter_recovery.csv",rw_pr)
csvwrite("created files/su1_parameter_recovery.csv",rw_su1)
csvwrite("created files/ph_parameter_recovery.csv",ph)
csvwrite("created files/modph_parameter_recovery.csv",modph)
%%
("Done with parameter recovery. Moving to model recovery")

%%
%model recovery
u(:,1) = data.u(:,1);
u(:,2) = data.u_cues(:,1);
seed = 111;
n_sim = 200;
%hgf_parameters
sim_omega2 = -8+(0.1+8)*rand(1,n_sim);
%sim_omega3 = unifrnd(-12,4,1,n_sim);
%sim_kappa = unifrnd(0,3,1,n_sim);
%rw
sim_alpha = 0.1 + (0.9-0.1) * rand(1,n_sim);
%sutton
sim_mu = 0.1 + (10-0.1)*rand(1,n_sim);
%ph
sim_S = 0.1 + (0.9-0.1)*rand(1,n_sim);
%response
sim_zeta = 0.5+(15-0.5)*rand(1,n_sim);


vec_data_hgf = [0,0,0,0];
vec_data_rw = [0,0,0,0];
vec_data_su1 = [0,0,0,0];
vec_data_ph = [0,0,0,0];


for i = 1:n_sim
    p_obs = sim_zeta(i);
    p_prc_hgf = [NaN 0 1 NaN 1 1 NaN 0 0 1 0 NaN sim_omega2(i) -6];

    p_prc_rw = [0.5,sim_alpha(i)];
    p_prc_su1 = [sim_mu(i),1,0.5,0.005];
    p_prc_ph = [0.5,0.5,sim_S(i)];

try
   %simulated models
    sim_hgf = tapas_simModel(u,...
                         'tapas_hgf_binary_pu_tgi_test',...
                         p_prc_hgf,...
                         'tapas_unitsq_sgm',...
                         p_obs,...
                         seed);

    sim_rw = tapas_simModel(u,...
                     'tapas_rw_binary',...
                     p_prc_rw,...
                     'tapas_unitsq_sgm',...
                     p_obs,...
                     seed);

    sim_su1 = tapas_simModel(u,...
                     'tapas_sutton_k1_binary',...
                     p_prc_su1,...
                     'tapas_unitsq_sgm',...
                     p_obs,...
                     seed);
                 
    sim_ph = tapas_simModel(u,...
                 'tapas_ph_binary',...
                 p_prc_ph,...
                 'tapas_unitsq_sgm',...
                 p_obs,...
                 seed);             
    
catch
    continue
end

try   
    %fit models to hgf data
    est_rw_data_hgf = tapas_fitModel(sim_hgf.y, sim_hgf.u, tapas_rw_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_hgf_data_hgf = tapas_fitModel(sim_hgf.y, sim_hgf.u, tapas_hgf_binary_tgi_config_wide_priors_2_level,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);    
    est_su1_data_hgf = tapas_fitModel(sim_hgf.y, sim_hgf.u, tapas_sutton_k1_binary_model_recovery_optimal_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_ph_data_hgf = tapas_fitModel(sim_hgf.y, sim_hgf.u, tapas_ph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    %best model:
    vec_data_hgf1 = [est_rw_data_hgf.optim.LME,est_hgf_data_hgf.optim.LME,est_su1_data_hgf.optim.LME,est_ph_data_hgf.optim.LME];
    
    vec_data_hgf = [vec_data_hgf;vec_data_hgf1];


    %fit models to rw data
    est_rw_data_rw = tapas_fitModel(sim_rw.y, sim_rw.u, tapas_rw_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_hgf_data_rw = tapas_fitModel(sim_rw.y, sim_rw.u, tapas_hgf_binary_tgi_config_wide_priors_2_level,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);    
    est_su1_data_rw = tapas_fitModel(sim_rw.y, sim_rw.u, tapas_sutton_k1_binary_model_recovery_optimal_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_ph_data_rw = tapas_fitModel(sim_rw.y, sim_rw.u, tapas_ph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    vec_data_rw1 = [est_rw_data_rw.optim.LME,est_hgf_data_rw.optim.LME,est_su1_data_rw.optim.LME,est_ph_data_rw.optim.LME];

    vec_data_rw = [vec_data_rw;vec_data_rw1];


    %fit models to su1 data
    est_rw_data_su1 = tapas_fitModel(sim_su1.y, sim_su1.u, tapas_rw_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_hgf_data_su1 = tapas_fitModel(sim_su1.y, sim_su1.u, tapas_hgf_binary_tgi_config_wide_priors_2_level,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);    
    est_su1_data_su1 = tapas_fitModel(sim_su1.y, sim_su1.u, tapas_sutton_k1_binary_model_recovery_optimal_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_ph_data_su1 = tapas_fitModel(sim_su1.y, sim_su1.u, tapas_ph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    vec_data_su11 = [est_rw_data_su1.optim.LME,est_hgf_data_su1.optim.LME,est_su1_data_su1.optim.LME,est_ph_data_su1.optim.LME];

    vec_data_su1 = [vec_data_su1;vec_data_su11];
    
    
   %fit models to pearce hall data
    est_rw_data_ph = tapas_fitModel(sim_ph.y, sim_ph.u, tapas_rw_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_hgf_data_ph = tapas_fitModel(sim_ph.y, sim_ph.u, tapas_hgf_binary_tgi_config_wide_priors_2_level,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);    
    est_su1_data_ph = tapas_fitModel(sim_ph.y, sim_ph.u, tapas_sutton_k1_binary_model_recovery_optimal_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    est_ph_data_ph = tapas_fitModel(sim_ph.y, sim_ph.u, tapas_ph_binary_config,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);
    
    vec_data_phh = [est_rw_data_ph.optim.LME,est_hgf_data_ph.optim.LME,est_ph_data_su1.optim.LME,est_ph_data_ph.optim.LME];

    vec_data_ph = [vec_data_ph;vec_data_phh];


catch
    continue
end

end
% saving the values
csvwrite("created files/simulated_zeta.csv",sim_zeta)
csvwrite("created files/simulated_omega2.csv",sim_omega2)
csvwrite("created files/simulated_alpha.csv",sim_alpha)
csvwrite("created files/simulated_mu.csv",sim_mu)
csvwrite("created files/simulated_S.csv",sim_S)


csvwrite("created files/hgf_data_model_recovery.csv",vec_data_hgf)
csvwrite("created files/rw_data_model_recovery.csv",vec_data_rw)
csvwrite("created files/su1_data_model_recovery.csv",vec_data_su1)
csvwrite("created files/ph_data_model_recovery.csv",vec_data_ph)

%%
("Done with model recovery.")
