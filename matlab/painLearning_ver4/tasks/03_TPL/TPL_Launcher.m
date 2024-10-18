function TPL_Launcher(scr, vars)
%function TPL(scr, vars)
%
% Project: TPL task
%
% Input: subNo      4-digit subject ID number
%
% Sets paths, and calls main.m
%
% Previous versions:
%
% ======================================================
%
% -------------- PRESS ESC TO EXIT ---------------------
%
% ======================================================
%
% Niia Nikolova
% Last edit: 21/07/2020


%% Initial settings
% Close existing workspace
close all; clc;

%% Define task specific vars
vars.exptName = 'TPL_v1-1';
vars.date_time = datestr(now,'ddmmyyyy_HHMMSS');
vars.DataFileName = strcat(vars.exptName, '_',vars.subIDstring, '_', vars.date_time);    % name of data file to write to
vars.UniqueFileName =  strcat(vars.exptName, '_',vars.subIDstring);
%% setup path
addpath(genpath('code'));
 
%% do checks
vars = TcsCheckStimulator(vars);

%% Confirm temperatures per subject
vars.cueModality = 2; % 1 = visual, 2 = auditory 
vars.stim_cold = 18;
vars.stim_warm = 42;

%% Run experiment
vars.startTask = tic;           % task start time
TPL_main(vars, scr);            % task script
endTask = toc(vars.startTask);  % task end time
disp(['TPL duration: ', num2str(round(endTask/60,1)), ' minutes']) % task duration

%% Restore path
rmpath(genpath('code'));
cd(vars.projdir)
