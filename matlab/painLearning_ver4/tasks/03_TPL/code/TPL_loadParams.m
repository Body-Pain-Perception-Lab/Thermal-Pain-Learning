%% Define parameters
%
% Project: TPL task
%
% Sets key parameters, called by main.m
%
% Francesca Fardo, Niia Nikolova
% Last edit: 24/06/2021

%% Key flags
vars.expN           = 2;
vars.emulate        = 0;                % 0 scanning, 1 testing
vars.Ratings        = 1;                % Ratings? (1 yes, 0 no)  
vars.nRatings       = 3;                % number of ratung scales
vars.InputDevice    = 2;                % Response method for conf rating. 1 - keyboard 2 - mouse
useEyeLink          = 0;                % Use EyeLink to record gaze & pupil    
vars.fixCrossFlag   = 1;                % Show fixation point
vars.pluxSynch      = 1;                % Get triggers for respiration
vars.breaksSelfPace = 1;                % Self paced breaks, yes = 1, no = 0
%% Fixed stimulation parameters
vars.timepoint = 0;         % timer 
vars.tcsData = [];          % variable to store stimulation data   
vars.skinT = 32;            % baseline temperature
vars.speed_ramp = 20;       % rate of temperature change from baseline to target
vars.return_ramp = vars.speed_ramp/2;      % rate of temperature change from target to baseline
vars.bp = 0;                % 1 = button press required
vars.save_file = 1;         % 1 = save tcs data

%% Stimulation temperatures (variable across participants)
% Start calculating morphs to use (this continues after block setup)
%[noThreshFlag, thresh, PMFstims]= getParticipantThreshold(vars);                 % get PMFstims = stimLevels at [.3 .5 .7 p(correct)]
%NJitterLevels = length(PMFstims);
%vars.NMorphJitters = NJitterLevels;
%vars.stim_cold = 20;%PMFstims(1:NJitterLevels/2);
%vars.stim_warm = 40;%PMFstims((NJitterLevels/2)+1:end);

%% Add TGI catch trials 
%catchTrialStim = '100';
%vars.catchTrialArray = [ones(1,vars.NCatchTrials/2), zeros(1,vars.NCatchTrials/2)];
%vars.catchTrialArray = mixArray(vars.catchTrialArray); 

%% Cueing
% Set temporal evolution of cue probablilities using
% createCueProbabilities.m

% Fetch trialSequence & cueSequence from saved file         
% [cueProbabilityOutput, blockParams, breaks] = createCueProbabilities(vars);

% Randomly select #1, 4 (later also 6, maybe 2[3 reversals, more
% difficult])
%randSequence = round(rand);
%if randSequence
if mod(vars.subNo,2)
    chosenSequence = 'sequence1_tgi.mat';
    vars.seqN = 1;
else
    chosenSequence = 'sequence4_tgi.mat';
    vars.seqN = 4;
end

[Output] = load(chosenSequence);
cueProbabilityOutput=Output.cueProbabilityOutput;
blockParams = Output.blockParams;
breaks = Output.breaks;
vars.cueProbabilityOutput = cueProbabilityOutput;
vars.NTrialsTotal        = length(cueProbabilityOutput);
vars.conditionSequence   = cueProbabilityOutput(:,2);            % conditions (1:5), % 1 cue_0 valid   2 cue_1 valid    3 cue_0 invalid    4 cue_1 invalid      5 non-predictive 
vars.faceSequence        = cueProbabilityOutput(:,4);            % sequence of face genders [0|1]
vars.cueSequence         = cueProbabilityOutput(:,5);            % sequence of cues [0|1]
vars.trialSequence       = cueProbabilityOutput(:,6);            % 1 valid, 2 invalid
vars.TempQual            = cueProbabilityOutput(:,10);           % 1 /warm, 0 cold
vars.cue0Prediction      = cueProbabilityOutput(:,12);           % %(is cue 0 predictive of warm (1) or cold (2) stimuli, on non-predictive(0)?)
vars.predictionTrialNext = cueProbabilityOutput(:,13);          % 1 if there is a prediction trial after the ITI of this trial
vars.ratings_yn          = cueProbabilityOutput(:,17); 
ratingsIndex             = find(vars.ratings_yn==1);
vars.changelocation      = ratingsIndex(5:5:length(ratingsIndex)-5);
vars.changelocation      = sort([vars.changelocation; breaks]);
vars.locationNumber = [];
a=[[1;vars.changelocation] [vars.changelocation(1); vars.changelocation(2:end); 307]];
replocation = a(:,2)-a(:,1);
locationN = repmat(1:3,1,10);
for i = 1:length(replocation)
    b = repmat(locationN(i),1,replocation(i));
    vars.locationNumber = [vars.locationNumber b];
end

vars.breaks                         = breaks;                  	% break AFTER this trial
vars.blockParams                    = blockParams;
vars.desiredBlockProbabilities      = blockParams(1,:);
vars.effectiveBlockProbabilities    = blockParams(2,:);
vars.blockLengths                   = blockParams(3,:);
stimJitterRepsByBlock               = blockParams(4,:);

vars.trialByTrialBlockVector    = [];
vars.TempStimuli                = [];                           % (Ntrials/2) x 2 array, [cold, warm]
blockCount                      = 1;

% Make a vector with block lengths & FaceMorphVals by trial
for thisBlock = 1:length(vars.blockLengths)
    tempBlockN = (blockCount .* ones(1,vars.blockLengths(thisBlock)));
    vars.trialByTrialBlockVector = [vars.trialByTrialBlockVector, tempBlockN];
    blockCount = blockCount + 1;
    clear tempBlockN

    thisBlockColdStims = mixArray(repmat(vars.stim_cold, 1, stimJitterRepsByBlock(thisBlock)));
    thisBlockWarmStims = mixArray(repmat(vars.stim_warm, 1, stimJitterRepsByBlock(thisBlock)));
    thisBlockAllStims = vertcat(thisBlockColdStims, thisBlockWarmStims);
    
    vars.TempStimuli = horzcat(vars.TempStimuli, thisBlockAllStims);
    
    % Add some extra stim morph levels b/c it doesn't always work out that
    % we have exactly = numbers of happy and angry stims - find a better
    % way to do this!
    if thisBlock == length(vars.blockLengths)
        thisBlockColdStims = mixArray(repmat(vars.stim_cold, 1, stimJitterRepsByBlock(thisBlock)));
        thisBlockWarmStims = mixArray(repmat(vars.stim_warm, 1, stimJitterRepsByBlock(thisBlock)));
        thisBlockAllStims = vertcat(thisBlockColdStims, thisBlockWarmStims);
        
        vars.TempStimuli = horzcat(vars.TempStimuli, thisBlockAllStims);
    end
    
end


%% Paths
% Faces
vars.TaskPath = fullfile('.', 'code', 'task');          
vars.StimFolder = fullfile('..', '..', 'stimuli', filesep);   %fullfile('.', 'stimuli', filesep);
vars.StimSize = 9;%9;in scanner                                      % DVA    (old 7, changed 07/07)                                  
%vars.StimsInDir = dir([vars.StimFolder, '*.tif']);      % list contents of 'stimuli' folder    

% Cues                                 
vars.CuesInDir = dir([vars.StimFolder, 'cue*']);      % list contents in 'stimuli' folder    

%% Task timing
vars.startExp =         tic;
vars.fixedTiming        = 0;    % Flag to force fixed timing for affect response  1 fixed, 0 self-paced (Conf rating always fixed, otherwise infinite!)
vars.RepeatMissedTrials = 0;
vars.CueT               = 1;
vars.jitter1_min        = 0.1;    % long variable ISI, 2-3 or 2-4 sec
vars.jitter1_max        = 0.2; 
vars.jitter1            = randInRange(vars.jitter1_min, vars.jitter1_max, [vars.NTrialsTotal,1]);
vars.RespT              = 3;  % max 3 sec to provide expectation
vars.PredD              = 1.5;  % Prediction display   
vars.RecT               = 5;    % record tcs data for 4 sec, then auditory cue
vars.StimT              = 10;   % stimulate for max 30 sec (but abort stimulation earlier
vars.ITI_min            = 1;    % long variable ISI, 2-3 or 2-4 sec
vars.ITI_max            = 2; 
vars.ITI                = randInRange(vars.ITI_min, vars.ITI_max, [vars.NTrialsTotal,1]);
vars.ConfT              = 5;    % present each rating scale for max X seconds
vars.ITI                = [];
vars.ITI_min_rateN      = 5;    % short variable ITI 
vars.ITI_max_rateN      = 6; 
vars.ITI_min_rateY      = 1;    % short variable ITI 
vars.ITI_max_rateY      = 2; 
for i = 1:length(vars.ratings_yn)
    if vars.ratings_yn == 0
        vars.ITI(i) = randInRange(vars.ITI_min_rateN, vars.ITI_max_rateN, 1);
    else
        vars.ITI(i) = randInRange(vars.ITI_min_rateY, vars.ITI_max_rateY, 1);
    end
end

vars.breakT             = 60;   % sec
%trialDur_min            = vars.CueT + vars.StimT + vars.RespT + vars.ConfT + vars.ISI_min + vars.ITI_min - 2;   % -2 if we expect ppts to respond quickly
%trialDur_max            = vars.CueT + vars.StimT + vars.RespT + vars.ConfT + vars.ISI_max + vars.ITI_max;
%trialDuration           = [trialDur_min, trialDur_max];

%% Plux synch variables
% Colours: White, Black
scr.pluxWhite     = WhiteIndex(scr.screenID);
scr.pluxBlack     = BlackIndex(scr.screenID);
%scr.pluxWhite     = WhiteIndex(scr.screenID);

% Duration
scr.pluxDur         = [2; 25];% [2;4]       % 2 frames - stim, 4 frames - response

% Size
% rows: x y width height
% provide size in cm and convert to pix
pluxRectTemp        = [0; 0; 1; 1];
multFactorW         = scr.resolution(1) ./ scr.MonitorWidth;
multFactorH         = scr.resolution(2) ./ scr.MonitorHeight;
scr.pluxRect(3)     = pluxRectTemp(3) .* multFactorW;
scr.pluxRect(4)     = pluxRectTemp(4) .* multFactorH;
scr.pluxRect = CenterRectOnPoint(scr.pluxRect,scr.pluxRect(3)/2,scr.resolution(2) - scr.pluxRect(4)/2);

 %% MR params
% vars.TR                 = 1.34;          % Seconds per volume
% vars.Dummies            = 4;             % Dummy volumes at start
% vars.Overrun            = 4;             % Dummy volumes at end
% vars.expectedDurationMax = vars.NTrialsTotal * trialDuration(2) / 60;
% vars.VolsPerExpmt       = round(vars.expectedDurationMax*60 /vars.TR) + vars.Dummies + vars.Overrun;
% 
% if vars.fixedTiming
%     disp(['Desired number of volumes: ', num2str(vars.VolsPerExpmt)]);
% else
%     disp(['Desired number of volumes (upper limit of session duration): ', num2str(vars.VolsPerExpmt)]);
% end
% disp('Press any key to continue.');
% pause;

%% Instructions
textSize = 35;
% switch vars.Ratings
%     
%     case 1
%         
%         switch vars.InputDevice
%             
%             case 1 % Keyboard
%                 vars.InstructionTask = 'Do you expect a cold or warm sensation? \n \n COLD - Left arrow key                         WARM - Right arrow key \n \n \n \n Then, you will be asked to rate the quality and intensity of your sensations. \n \n How cold? (1), How warm? (2), How burning? (3), and How painful? (4). \n \n Press SPACE to start...';
%                 vars.InstructionRating{1} = 'How cold does it feel? \n \n Unsure (1)     Sure (2)     Very sure (3)';
%                 vars.InstructionRating{2} = 'How warm does it feel? \n \n Unsure (1)     Sure (2)     Very sure (3)';
%                 vars.InstructionRating{3} = 'How burning does it feel? \n \n Unsure (1)     Sure (2)     Very sure (3)';
%                 vars.InstructionRating{4} = 'How painful does it feel? \n \n Unsure (1)     Sure (2)     Very sure (3)';
% 
%             case 2 % Mouse
%                 vars.InstructionTask = 'Do you expect a cold or warm sensation? \n \n COLD - Left click                          WARM - Right click \n \n \n \n Press SPACE to start...';
%                 vars.InstructionRatings = 'Please complete the ratings using the mouse. Left click to confirm.';
%         end
%     case 0
%         switch vars.InputDevice
%             
%             case 1 % Keyboard
%                 vars.InstructionTask = 'Decide if the face presented on each trial is angry or happy. \n \n ANGRY - Left arrow key                         HAPPY - Right arrow key \n \n \n \n Press SPACE to start...';
%             case 2 % Mouse
%                 vars.InstructionTask = 'Decide if the face presented on each trial is angry or happy. \n \n ANGRY - Left button                         HAPPY - Right button \n \n \n \n Press SPACE to start...';
%         end
% end

switch vars.language
    case 1 % English
        vars.InstructionTask = 'Please position the thermode to location 1. \n \n Do you expect a cold or warm sensation? \n \n COLD - Left click                          WARM - Right click \n \n \n \n Press SPACE to start...';
        vars.InstructionQ = 'Cold (L)     or     Warm (R)?';
        vars.InstructionRatings = 'Please complete the ratings using the mouse. Left click to confirm.';
        vars.InstructionPause{1} = 'Take a short break... \n \n The experiment will continue in ...';
        vars.InstructionPause{2} = 'Take a short break... \n \n Press SPACE to continue ...';
        vars.InstructionEnd = 'You have completed the session. Thank you!';
        value = {'cold', 'warm', 'burning', 'painful'};
        valueBold = {'COLD', 'WARM', 'BURNING', 'PAINFUL'};
        for ii = 1:4
            %vars.InstructionRating{ii} = ['Do you feel a ' value{ii} ' sensation now? \n \n  If so, please rate how intensily you perceive this sensation. \n \n' valueBold{ii} '?'];
             vars.InstructionRating{ii} = ['Please rate the most intense ' value{ii} ' sensation  you perceived. \n \n' valueBold{ii} '?'];
        end
        vars.ConfEndPoins = {'Not at all', 'Extremely'};
        vars.InstructionMove = 'Please move the thermode to location';
        
    case 2 % Danish
        vars.InstructionTask = 'Placer venligst sensoren på lokation 1. \n \n Forventer du en kold eller varm fornemmelse? \n \n KOLD - Venstreklik                          VARM - Højreklik \n \n \n \n Tryk på mellemrumstatsten for at begynde...';
        vars.InstructionQ = 'Kold (V)     eller     Varm (H)?';
        vars.InstructionRatings = 'Brug venligst musen til at bedømme. Venstreklik for at bekræfte.';
        vars.InstructionPause{1} = 'Tag en kort pause... \n \n Experimentet fortsætter om ...';
        vars.InstructionPause{2} = 'Tag en kort pause... \n \n Tryk på MELLEMRUM for at fortsætte ...';
        vars.InstructionEnd = 'Du har nu gennemført denne session. Mange tak!';
        value = {'kold', 'varm', 'brændende', 'smertefuld'};
        valueBold = {'KOLD', 'VARM', 'BRÆNDENDE', 'SMERTEFULD'};
        for ii = 1:4
            vars.InstructionRating{ii} = ['Bedøm venligst den mest intense ' value{ii} ' fornemmelse du bemærkede. \n \n' valueBold{ii} '?'];
        end
        vars.ConfEndPoins = {'Overhovedet ikke', 'Ekstremt'};
        vars.InstructionMove = 'Placer venligst sensoren på lokation';
end

