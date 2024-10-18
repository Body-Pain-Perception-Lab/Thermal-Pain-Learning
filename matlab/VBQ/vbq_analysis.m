%% setup paths
clear all; close all
addpath('/home/jespere/Thermal-Pain-Learning/matlab/VBQ')
addpath('/mnt/fast_scratch/toolboxes/MATLAB/spm12')
scriptdir= '/home/jespere/VMP_pipelines/code/mpm/analysis';

addpath(genpath(scriptdir))

basepath = '/home/jespere/Thermal-Pain-Learning_final/matlab/VBQQ';
csvfile =  fullfile(basepath, 'fmri_data.csv');


%% read in regressors
thermaldata = import_data(csvfile);
thermaldata = removevars(thermaldata, 'Var1');


%% format Regressor file

names = thermaldata.Properties.VariableNames;
R = table2array(thermaldata);

rowsWithNaNs = any(isnan(R),2);
R(rowsWithNaNs,:) = [];

%% save regressor file

save(fullfile(basepath, 'regressors.mat'), 'R', 'names')

%% clean up, load, and run analysis
analysis_name = 'VBQ_thermosensory';
regressors_file = 'regressors.mat';
%%
regressors = fullfile(basepath, 'regressors.mat');
mpm_vbq_runstats(analysis_name, basepath, regressors)

%% make contrasts0
analysis_dirs = dir(fullfile(basepath,analysis_name));

r_file_path = fullfile(basepath,[regressors_file(1:end-4),'_tiv.mat'] );

spm_file_path = fullfile(basepath,analysis_name, analysis_dirs(5).name,'SPM.mat');

make_regression_contrasts(r_file_path, spm_file_path)

spm_file_path = fullfile(basepath,analysis_name, analysis_dirs(3).name,'SPM.mat');

make_regression_contrasts(r_file_path, spm_file_path)

spm_file_path = fullfile(basepath,analysis_name, analysis_dirs(4).name,'SPM.mat');

make_regression_contrasts(r_file_path, spm_file_path)