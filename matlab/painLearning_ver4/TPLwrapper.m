function TPLwrapper(whichPart)
% function TPLwrapper(whichPart)
%
% Project: Thermosensation Pain Learning (TPL) task behavioural, for VMP 2.0 (summer/Fall 2021)
%
% Input: whichPart  optional argument to only run one of the TPL component tasks
%       1   QST Thresholding
%       2   PSI Thresholding
%       3   TPL Task
%
% Sets paths, and calls functions
%
% Niia Nikolova and Francesca Fardo
% Last edit: 30/06/2021

%% TPL tasks wrapper

% Close existing workspace
close all; clc;

%% Define general vars across tasks
vars.projdir = pwd;
vars.devFlag  = 0;              % Development flag 1. Set to 1 when developing the task, will optimize stim size for laptop, not hide cursor
vars.stimFlag = 1;              % Development flag 2. Set to 0 when developing the task without a stimulator
vars.subNo = input('What is the subject number (e.g. 0001)?   ');
vars.language = input('Which language: English (1) or Danish (2)?   ');
vars.subIDstring = sprintf('%04d', vars.subNo);
vars.startTrialN = 1;

% check for data dir
if ~exist('data', 'dir')
    mkdir('data')
end

% Define subject No if the value is missing 
if isempty(vars.subNo)
    vars.subNo = 9999; % debugging                                            
end

% Check if the subject folder already exists in data dir
vars.OutputFolder = fullfile(vars.projdir, 'data', ['sub_',vars.subIDstring], filesep);
if ~exist(vars.OutputFolder, 'dir') 
    mkdir(vars.OutputFolder)
end

%% Set up paths
addpath(vars.OutputFolder);
addpath(genpath('stimuli'));
addpath(genpath('..\LibTcsMatlab2021a'));

participant.MetaDataFileName = strcat(vars.subIDstring, '_metaData'); 
participant.partsCompleted = zeros(1,3);

%% Check that PTB is installed
[oldLevelScreen, oldLevelAudio] = checkPTBinstallation;
%% Open a PTB window
scr.ViewDist = 56; 
[scr] = displayConfig(scr);
AssertOpenGL;
if vars.devFlag
    [scr.win, scr.winRect] = PsychImaging('OpenWindow', scr.screenID, scr.BackgroundGray, [0 0 1000 1000]); %,[0 0 1920 1080] mr screen dim
else
    [scr.win, scr.winRect] = PsychImaging('OpenWindow', scr.screenID, scr.BackgroundGray); %,[0 0 1920 1080] mr screen dim
end
% PsychColorCorrection('SetEncodingGamma', scr.win, 1/scr.GammaGuess);

% Set text size, dependent on screen resolution
if any(logical(scr.winRect(:)>3000))       % 4K resolution
    scr.TextSize = 65;
else
    scr.TextSize = 28;
end
Screen('TextSize', scr.win, scr.TextSize);

% Set priority for script execution to realtime priority:
scr.priorityLevel = MaxPriority(scr.win);
Priority(scr.priorityLevel);

% Determine stim size in pixels
scr.dist        = scr.ViewDist;
scr.width       = scr.MonitorWidth;
scr.resolution  = scr.winRect(3:4);                    % number of pixels of display in horizontal direction

%% 01 Run thresholding task
if ((nargin < 1) || (whichPart==1)) %&& (participant.partsCompleted(taskN) == 0)
    vars.taskN = 1;
    cd(fullfile('.', 'tasks', '01_qstThreshold'))
    addpath(genpath('code'))
    qstThreshold_Launcher(scr, vars); % Launcher
    participant.partsCompleted(1) = 1;
    % Save metadata
    save(fullfile(vars.OutputFolder, ['sub_', participant.MetaDataFileName]), 'participant');
    % Continue to next task 
    if (nargin < 1)
        goOn1 = input('QST Threshold task completed. Continue to Demo? 1-yes, 0-no ');
        if ~goOn1
            return
        end
    end
end

%% 02 Demo Stimuli and ratings

if ((nargin < 1) || (whichPart==2)) %&& (participant.partsCompleted(2) == 0)
    vars.taskN = 2;
    cd(fullfile('.', 'tasks', '02_demoStimuli'))    
    demoStimuli_Launcher(scr,vars); % Launcher
    % if vars.RunSuccessfull
    participant.partsCompleted(2) = 1;
    % Save metadata
    save(fullfile(vars.OutputFolder, ['sub_', participant.MetaDataFileName]), 'participant');
    % Continue to next task 
    if (nargin < 1)
        goOn1 = input('Demo stimuli completed. Continue to TPL task? 1-yes, 0-no ');
        if ~goOn1
            return
        end
    end
end

%% 02 PSI Threshold
% 
% if ((nargin < 1) || (whichPart==2)) %&& (participant.partsCompleted(2) == 0)
%     vars.taskN = 2;
%     cd(fullfile('.', 'tasks', '02_psiThreshold'))    
%     psiThreshold_Launcher(scr,vars); % Launcher
%     % if vars.RunSuccessfull
%     participant.partsCompleted(2) = 1;
%     % Save metadata
%     save(fullfile(vars.OutputFolder, ['sub_', participant.MetaDataFileName]), 'participant');
%     % Continue to next task 
%     if (nargin < 1)
%         goOn1 = input('PSI Threshold task completed. Continue to TPL? 1-yes, 0-no ');
%         if ~goOn1
%             return
%         end
%     end
% end
% 
%% 03 Run TPL
if ((nargin < 1) || (whichPart==3)) %&& (participant.partsCompleted(taskN) == 0)
    vars.taskN = 3;
    % Run the task
    cd(fullfile('.', 'tasks', '03_TPL'))    
    TPL_Launcher(scr, vars); % Launcher
    % if vars.RunSuccessfull
    participant.partsCompleted(3) = 1;
    % Save metadata
    save(fullfile(vars.OutputFolder, ['sub_', participant.MetaDataFileName]), 'participant');
end

%% Finish up
% Copy data files to 1_VMP_aux
% copy2VMPaux(participant.subNo);

% Close screen etc
rmpath(genpath('code'));
rmpath(vars.OutputFolder);
sca;
ShowCursor;
fclose('all');
Priority(0);
ListenChar(0);          % turn on keypresses -> command window

%% Restore PTB verbosity
Screen('Preference', 'Verbosity', oldLevelScreen);
PsychPortAudio('Verbosity', oldLevelAudio);

%% Copy data to aux
try
    copy2VMPaux_TPL(vars.subNo)
    disp('Data successfully copied to AUX')
catch
    disp('Data were not copied to AUX')
end
