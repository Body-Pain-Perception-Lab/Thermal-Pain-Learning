function demoStimuli_Launcher(scr, vars)
%
% Project: Demonstration of the stimuli uses in the pain learning task
%
% Francesca Fardo 
% Last edit: 03/07/2021

%% Initial settings
% Close existing workspace
close all; clc;

%% Define task specific vars 
vars.exptName = 'demoStim_v1-1';
vars.date_time = datestr(now,'ddmmyyyy_HHMMSS');
vars.DataFileName = strcat(vars.exptName, '_',vars.subIDstring, '_', vars.date_time);    % name of data file to write to
vars.UniqueFileName =  strcat(vars.exptName, '_',vars.subIDstring);
%% setup path
addpath(genpath('code'));
 
%% Do checks 
vars = TcsCheckStimulator(vars);

%% Start experiment
vars.startTask = tic;           % task start time
demoStimuli_main(vars, scr);       % task script
endTask = toc(vars.startTask);  % task end time
disp(['Demo Stimuli duration: ', num2str(round(endTask/60,1)), ' minutes']) % task duration

%% Restore path
rmpath(genpath('code'));
cd(vars.projdir)
