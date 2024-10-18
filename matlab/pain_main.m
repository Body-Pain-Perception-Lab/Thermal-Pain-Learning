
rng(123)
% Get the directory of the main script
mainScriptPath = fileparts(which('pain_main'));

% Construct the path to the 'scripts' folder
scriptsFolderPath = fullfile(mainScriptPath, 'scripts');

% Generate a path string with all subfolders under 'scripts'
allSubfoldersPath = genpath(scriptsFolderPath);

% Add all subfolders to the MATLAB path
addpath(allSubfoldersPath);
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

%%
%now we might want to fit our model the participants' actual responses:
%(i.e fit the model to actual data)
data = datafeel;
id = 2; 
u = data.u(:,id);
u1 = data.u(:,id);
u1(:,2) = data.u_cues(:,id);
y = data.y(:,id);

est = tapas_fitModel(y, u1, tapas_hgf_binary_pu_tgi_config_test,tapas_unitsq_sgm_config, tapas_quasinewton_optim_config);

tapas_hgf_binary_pu_tgi_plotTraj(est)
data.stormdb1(id)
%%
%fitting all models to all paraticipants
hgf_feel_2 = fit_hgf_subject_model_final_with_tgi_all1(datafeel, {'tapas_hgf_binary_pu_tgi_test_2_level_config'}, tapas_unitsq_sgm_config(), tapas_quasinewton_optim_config());

sk_feel = fit_hgf_subject_model_final_with_tgi_all1(datafeel, "tapas_sutton_k1_binary_config_real", tapas_unitsq_sgm_config(), tapas_quasinewton_optim_config());

rw_feel = fit_hgf_subject_model_final_with_tgi_all1(datafeel, "tapas_rw_binary_config_real", tapas_unitsq_sgm_config(), tapas_quasinewton_optim_config());

ph_feel = fit_hgf_subject_model_final_with_tgi_all1(datafeel, "tapas_ph_binary_config", tapas_unitsq_sgm_config(), tapas_quasinewton_optim_config());

%%
% model comparison
bor = hgf_feel_2
bor(:,2) = rw_feel
bor(:,3) = sk_feel
bor(:,4) = ph_feel
clear lmes
[posterior,out] = model_compar(bor)
%%
% saving
qq(:,1) = out.Ef
qq(:,2) = out.ep
csvwrite('created files/model_comparison2.csv',qq)
%Winning model is the 2 level hgf so we extract those trajectories

%%
% make trajectory data
 d1 = [];
 response = hgf_feel_2;
 data = datafeel;
for i = 1:length(data.stormdb1)
   q = response(i); 
   a = q{1, 1}.traj.mu;
   b = q{1, 1}.traj.sa;
   c = q{1, 1}.traj.muhat;
   c1 = q{1, 1}.traj.sahat;
   c2 = q{1, 1}.traj.v;
   c3 = q{1, 1}.traj.w;
   c4 = q{1, 1}.traj.da;
   c5 = q{1, 1}.traj.ud;
   c6 = q{1, 1}.traj.psi;
   c7 = q{1, 1}.traj.epsi;
   c8 = q{1, 1}.traj.wt;
   d = (repelem(data.stormdb1(i),306))';
   e = [a,b,c,c1,c2,c3,c4,c5,c6,c7,c8,d];
   f = (repelem(1,306))';
   g = [e,f];
   d2 = g;
   d1 = [d1;d2];
end

%save it
csvwrite("created files/binary_response_feel_22.csv",d1)

%%
%extracting subject specific parameters of the models to Brain analysis
d1 = [];
 response = hgf_feel_2;
 data = datafeel;
for i = 1:length(data.stormdb1)
   q = response(i); 
   a = q{1, 1}.p_prc.om(2);
   b = q{1, 1}.p_obs.ze(1);
   c = data.stormdb1(i);
   d2 = [a,b,c];
   d1 = [d1;d2];
end
%saving
csvwrite("created files/decision_making_parameter_compar2.csv",d1)