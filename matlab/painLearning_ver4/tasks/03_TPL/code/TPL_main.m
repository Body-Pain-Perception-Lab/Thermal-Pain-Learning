function TPL_main(vars, scr)
%function main(vars, scr)
%
% Project: TPL task
%
% Main experimental script. 
%
% Presents a cue, followed by a prediction display whoch queries whether the participant
% expects a warm or cold stimulation. This is followed by an actual stimulation and possibly 
% ratings (on cold, warm and burning sensations). The association between cues and stimuli is 
% predetermined accross probabilistic blocks.
%
% Cold and warm temperatures are determined by the participants PMF (?) 
%
% Input:
%   vars        struct with key parameters (most are deifne in loadParams.m)
%   scr         struct with screen / display settings
%
%
% Francesca Fardo, Niia Nikolova 
% Last edit: 04/08/2021

% Load the parameters
TPL_loadParams;

%% Fill in Results structure
% Results.sequence                       = vars.cueProbabilityOutput;
% Results.breaks                         = vars.breaks;                  	% break AFTER this trial
% Results.blockParams                    = vars.blockParams;
% Results.desiredBlockProbabilities      = vars.desiredBlockProbabilities;
% Results.effectiveBlockProbabilities    = vars.effectiveBlockProbabilities;
% Results.blockLengths                   = vars.blockLengths;
% Results.trialByTrialBlockVector        = vars.trialByTrialBlockVector; 
% Results.subID                          = repmat(vars.subNo,vars.NTrialsTotal,1); % subject number
% Results.trialN                         = 1:vars.NTrialsTotal;   
% Results.baseline                       = repmat(vars.skinT,vars.NTrialsTotal,1); 
% Results.cueModality                    = repmat(vars.cueModality,vars.NTrialsTotal,1);
% Results.cueSequence                    = vars.cueSequence;
% Results.conditionSequence              = vars.conditionSequence;

%% Define Results struct
uniqueFilename = strcat(vars.OutputFolder, vars.UniqueFileName, '.mat');
if ~exist(uniqueFilename)
    DummyDouble = ones(vars.NTrialsTotal,1).*NaN;
    %DummyString = strings(vars.NTrialsTotal,1);
    Results = struct('SOT_trial',       {DummyDouble}, ... % trial starts
                     'SOT_jitter',      {DummyDouble}, ... % jitter between cue off and prediction on
                     'SOT_cueOn',       {DummyDouble}, ... % cue on
                     'SOT_cueOff',      {DummyDouble}, ... % cue off
                     'SOT_predOn',      {DummyDouble}, ... % prediction starts
                     'SOT_predResp',    {DummyDouble}, ... % prediction timing
                     'SOT_stimOn',      {DummyDouble}, ... % stimulation on timing
                     'SOT_stimOff',     {DummyDouble}, ... % stimulation off timing
                     'SOT_vasOn',       {DummyDouble}, ... % vas on timing
                     'SOT_vasOff',      {DummyDouble}, ... % vas off timing
                     'SOT_iti',         {DummyDouble}, ... % iti timing
                     'predResp',        {DummyDouble}, ... % prediction timing
                     'predRT',          {DummyDouble}, ... % prediction reaction time
                     'predAcc',         {DummyDouble}, ... % prediction accuracy
                     'stim',            {DummyDouble}, ... % stimulation type 0 = c, 1 = w, 2 = tgi
                     'targetT',         {DummyDouble}, ... % stimulation temperatures
                     'vasYN',           {DummyDouble}, ... % whether the trial includes ratings
                     'vasResp',         {DummyDouble}, ... % vas response
                     'vasRT',           {DummyDouble}, ... % vas reaction time
                     'jitterDuration',  {DummyDouble}, ... % jitter duration
                     'stimDuration',    {DummyDouble}, ... % stimulation duration
                     'itiDuration',     {DummyDouble}, ... % iti duration
                     'trialDuration',   {DummyDouble}, ... % duration of each trial
                     'trialValidity',   {DummyDouble}, ... % whether RT and VAS are completed
                     'trialSuccess',    {DummyDouble});    % 1 = no missing data, 0 = some missing data
                     
                     
else
    vars.confirmedSubjN = input('Subject already exists. Do you want to continue anyway (yes = 1, no = 0)?    ');
    if vars.confirmedSubjN
        load(uniqueFilename,'Results')
        vars.startTrialN = input('Define the trial number to restart from?   ');
        vars.date_time = datestr(now,'ddmmyyyy_HHMMSS');
        vars.DataFileName = strcat(vars.exptName, '_',vars.subIDstring, '_', vars.date_time);    % name of data file to write to
    else
        return
    end
end

%% Ask to confirm temperatures
if vars.stimFlag 
    changeCold = input(['Do you want to change the cold temperature of ', num2str(vars.stim_cold), ' C. Type 1 (yes) or 0 (no)?  ']);
    if changeCold; vars.stim_cold = input('Define new cold temperature:    '); end
    changeWarm = input(['Do you want to change the warm temperature of ', num2str(vars.stim_warm), ' C. Type 1 (yes) or 0 (no)?  ']);
    if changeWarm; vars.stim_warm = input('Define new warm temperature   '); end
    disp(['Running task with : ' num2str(vars.stim_cold) ' C (cold)'])
    disp(['Running task with : ' num2str(vars.stim_warm) ' C (warm)'])    
end
%% Keyboard & keys configuration
[keys] = keyConfig();

% Reseed the random-number generator
SetupRand;

%% If this participant does not have a TPL (FAD) threshold, and we want to abort
% if noThreshFlag
%     disp('No threshold found. Terminating TPL task.');
%     vars.RunSuccessfull = 0;
%     vars.Aborted = 1;
%     experimentEnd(keys, Results, scr, vars);
%     return  
% end

%% Run tutorial?
TPL_Tutorial_YN = input('Do you want to run the tutorial (yes=1, no=0)?   ');
if TPL_Tutorial_YN
    TPL_Tutorial(scr, keys, vars);
end

%% Prepare to start
try
    %% Check if window is already open (if not, open screen window) 
    if ~isfield(scr, 'win')
        AssertOpenGL;
        [scr.win, scr.winRect] = PsychImaging('OpenWindow', scr.screenID, scr.BackgroundGray); %,[0 0 1920 1080] mr screen dim
%         PsychColorCorrection('SetEncodingGamma', scr.win, 1/scr.GammaGuess);
        % Set text size, dependent on screen resolution
        if any(logical(scr.winRect(:)>3000))       % 4K resolution
            scr.TextSize = 65;
        else
            scr.TextSize = textSize;
        end
        Screen('TextSize', scr.win, scr.TextSize);
        
        % Set priority for script execution to realtime priority:
        scr.priorityLevel = MaxPriority(scr.win);
        Priority(scr.priorityLevel);
        
        % Determine stim size in pixels
        scr.dist        = scr.ViewDist;
        scr.width       = scr.MonitorWidth;
        scr.resolution  = scr.winRect(3:4);   
    end
                   
    %% Plux parameters
    scr.hz          = Screen('NominalFrameRate', scr.win); 
    scr.pluxDurSec  =  scr.pluxDur / scr.hz;
    
    %% Dummy calls to prevent delays
    vars.ValidTrial = zeros(1,4);
    vars.RunSuccessfull = 0;
    vars.Aborted = 0;
    vars.Error = 0;
    thisTrialCorrect = 0;
    WaitSecs(0.1);
    GetSecs;
    vars.Resp = NaN;
    vars.ConfResp = NaN;%888;
    vars.abortFlag = 0;
    WaitSecs(0.500);
    [~, ~, keys.KeyCode] = KbCheck;
    
    %% Initialise EyeLink
    if useEyeLink
        vars.EyeLink = 1;
        
        % check for eyelink data dir
        if ~exist('./data/eyelink', 'dir')
            mkdir('./data/eyelink')
        end
        
        [vars] = ELsetup(scr, vars);
    end

    
    %% Show task instructions
    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
    DrawFormattedText(scr.win, uint8([vars.InstructionTask]), 'center', 'center', scr.TextColour);
    [~, ~] = Screen('Flip', scr.win);
    new_line;

    % Wait for Space press or ESC
    while keys.KeyCode(keys.Space) == 0
         [~, ~, keys.KeyCode] = KbCheck;
         WaitSecs(0.001);
         
         if keys.KeyCode(keys.Escape)==1
            % Save, mark the run
            vars.RunSuccessfull = 0;
            vars.Aborted = 1;
            experimentEnd(keys, Results, scr, vars);
            return
         end         
    end
    
    Results.SessionStartT = GetSecs;            % session start
    disp('Starting experiment')

    if useEyeLink
        Eyelink('message','STARTEXP');
    end
    
    %% Run through trials
    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
    [~, ~] = Screen('Flip', scr.win);               
    WaitSecs(3);            % pause before experiment start
    thisTrial = vars.startTrialN; % trial counter (user defined)
    %warmCounter = 1;
    %coldCounter = 1;
    %tgiCounter = 1;
    endOfExpt = 0;
    
    %% Draw plux trigger -- START
    if vars.pluxSynch
        Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
        [~, ~] = Screen('Flip', scr.win);
        WaitSecs(2)
    end
    [~, ~] = Screen('Flip', scr.win);
    
    %% Experiment
    while endOfExpt ~= 1       % General stop flag for the loop     
        
        Results.SOT_trial(thisTrial) = GetSecs - Results.SessionStartT; % trial starts
        if useEyeLink
            % EyeLink:  this trial
            startStimText = ['Trial ' num2str(thisTrial) ' start'];
            Eyelink('message', startStimText);
        end
        [~, TrialOn] = Screen('Flip',scr.win);
        
        % Will the outcome be C or W?
        thisTrialStim = vars.TempQual(thisTrial);         % 1 Warm, 0 Cold, 2 TGI
        Results.stim(thisTrial) = thisTrialStim;          % save this info in Results
        if thisTrialStim == 0
            tempPresented = 'cold';
        elseif thisTrialStim == 1
            tempPresented = 'warm';
        elseif thisTrialStim == 2
            tempPresented = 'tgi';
        end
        
        if vars.pluxSynch
            % if were in the first pluxDurationSec seconds, draw the rectangle
            if thisTrialStim == 0 &&((GetSecs - TrialOn) <= scr.pluxDurSec(2)) % Cold
                Screen('FillRect', scr.win, scr.pluxBlack, scr.pluxRect);
            else % Warm or TGI
                Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
            end
        end
        Screen('Flip',scr.win);
        
        %% Present cue
        thisCue = vars.cueSequence(thisTrial);
        disp(['Trial # ', num2str(thisTrial), '. Cue: ', num2str(thisCue)]);
        switch vars.cueModality 
            case 1 % visual
                thisTrialCue = ['cue_', num2str(thisCue), '.tif'];

                % Read stim image for this trial into matrix 'imdata'
                CueFilePath = strcat(vars.StimFolder, thisTrialCue);
                ImDataOrig = imread(char(CueFilePath));
                ImData = imresize(ImDataOrig, [StimSizePix NaN]);           % Adjust image size to StimSize dva in Y dir

                % Make texture image out of image matrix 'imdata'
                ImTex = Screen('MakeTexture', scr.win, ImData);

                % Draw texture image to backbuffer
                Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                Screen('DrawTexture', scr.win, ImTex);
                [~, CueOn] = Screen('Flip', scr.win);

                Results.SOT_cueOn(thisTrial) = CueOn - Results.SessionStartT;
                if useEyeLink
                    % EyeLink:  cue on
                    startStimText = ['Trial ' num2str(thisTrial) ' cue on'];
                    Eyelink('message', startStimText);
                end

                % While loop to show stimulus until CueT seconds elapsed.
                while (GetSecs - CueOn) <= vars.CueT

                    % Draw the cue screen
                    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                    Screen('DrawTexture', scr.win, ImTex);

                    % Draw plux trigger -- CUE
                    if vars.pluxSynch
                        % if were in the first pluxDurationSec seconds, draw the rectangle
                        if thisCue == 0     &&((GetSecs - CueOn) <= scr.pluxDurSec(1)) 
                            Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
                        elseif thisCue == 1 &&((GetSecs - CueOn) <= scr.pluxDurSec(1))
                            Screen('FillRect', scr.win, scr.pluxBlack, scr.pluxRect);
                        end
                    end

                    % Flip screen
                    Screen('Flip', scr.win);

                    % KbCheck for Esc key
                    if keys.KeyCode(keys.Escape)==1
                        % Save, mark the run
                        vars.RunSuccessfull = 0;
                        vars.Aborted = 1;
                        experimentEnd(keys, Results, scr, vars);
                        return
                    end

                    [~, ~, keys.KeyCode] = KbCheck;
                    WaitSecs(0.001);
                end

                [~, ~] = Screen('Flip', scr.win);            % clear screen
                Screen('Close', ImTex);                      % Close the image texture

              case 2 %auditoryCue
                if ~thisCue
                    cueFreq = 400;
                else
                    cueFreq = 1600;
                end
                Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                [~, CueOn] = Screen('Flip', scr.win);
                Results.SOT_cueOn(thisTrial) = CueOn - Results.SessionStartT;
                playBeep(cueFreq) % auditory cue
                
        end
        
        % cue off
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        [~, CueOff] = Screen('Flip', scr.win);
        Results.SOT_cueOff(thisTrial) = CueOff - Results.SessionStartT;
        
        if useEyeLink
            % EyeLink:  cue off
            startStimText = ['Trial ' num2str(thisTrial) ' cue off'];
            Eyelink('message', startStimText);
        end
        
        %% ISI (between cue off and display of prediction question)
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        if vars.fixCrossFlag
            scr = drawFixation(scr);
        end
        [~, StartJitter] = Screen('Flip', scr.win);
        
        Results.SOT_jitter(thisTrial) = StartJitter - Results.SessionStartT;
        Results.jitterDuration = vars.jitter1(thisTrial);
        if useEyeLink
            % EyeLink:  ITI
            startStimText = ['Trial ' num2str(thisTrial) ' jitter start'];
            Eyelink('message', startStimText);
        end
        
        % Present the gray screen for ITI duration
        while (GetSecs - StartJitter1) <= vars.jitter1(thisTrial)
            
            if keys.KeyCode(keys.Escape)==1
                % Save, mark the run
                vars.RunSuccessfull = 0;
                vars.Aborted = 1;
                experimentEnd(keys, Results, scr, vars);
                return
            end
        end
        
        [~, ~, keys.KeyCode] = KbCheck;
        WaitSecs(0.001);
       
        %% Prediction
       
        % Cold (L arrow) or Warm (R arrow)?
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        DrawFormattedText(scr.win, uint8([vars.InstructionQ]), 'center', 'center', scr.TextColour);
        
        [~, vars.StartRT] = Screen('Flip', scr.win);
        
        if useEyeLink
            % EyeLink: response
            startStimText = ['Trial ' num2str(thisTrial) ' response screen on'];
            Eyelink('message', startStimText);
        end
        
        % Fetch the participant's response, via keyboard or mouse
        [vars] = getResponse(keys, scr, vars);
        
        Results.SOT_predOn(thisTrial) = vars.StartRT - Results.SessionStartT;
        Results.SOT_predResp(thisTrial) = vars.EndRT - Results.SessionStartT;
        
        if vars.abortFlag               % Esc was pressed
            Results.predResp(thisTrial) = 9;
            % Save, mark the run
            vars.RunSuccessfull = 0;
            vars.Aborted = 1;
            experimentEnd(keys, Results, scr, vars);
            return
        end
        
        % Time to stop? (max # trials reached)
        if (thisTrial == vars.NTrialsTotal)
            endOfExpt = 1;
        end
        
        % Compute response time
        RT = (vars.EndRT - vars.StartRT);
        
        % Compute accuracy
        if thisTrialStim == 1   % 1 = warm stimulation
            if vars.Resp == 1   % 1 = warm prediction
                thisTrialCorrect = 1;
            elseif vars.Resp == 0 % 0 = cold prediction
                thisTrialCorrect = 0;
            end
        elseif thisTrialStim == 0 % 0 = cold stimulation
            if vars.Resp == 1     % 1 = warm prediction
                thisTrialCorrect = 0;
            elseif vars.Resp == 0 % 0 = cold prediction
                thisTrialCorrect = 1; 
            end
        elseif thisTrialStim == 2 % 2 = TGI stimulation´
            thisTrialCorrect = 2; 
        end
        
        % Write trial result to file
        Results.predResp(thisTrial) = vars.Resp;
        Results.predRT(thisTrial) = RT;
        Results.predAcc(thisTrial) = thisTrialCorrect;
        
        %% Prediction Display - display what the participant chose for a fixed time, to increase anticipation
        % wait after response defined in vars.PredD
        % edit the getResponse function to modify the prediction display
        
        %% Define stimulation temperatures
        if thisTrialStim == 1                    % Warm
            %warmCounter = warmCounter + 1;
            allWarmTemps = [vars.stim_warm vars.stim_warm vars.stim_warm vars.skinT vars.skinT; vars.skinT vars.stim_warm vars.stim_warm vars.stim_warm vars.skinT; vars.skinT vars.skinT vars.stim_warm vars.stim_warm vars.stim_warm];
            warmTempsIndex = randi(size(allWarmTemps,1)); % choose randomly among a predefined set of warm temperatures
            Results.targetT(thisTrial,1:5) = allWarmTemps(warmTempsIndex,:); % Define stimulation temperature for this trial
        elseif thisTrialStim == 0                % Cold
            %coldCounter = coldCounter + 1;
            allColdTemps = [vars.stim_cold vars.stim_cold vars.skinT vars.skinT vars.skinT; vars.skinT vars.stim_cold vars.stim_cold vars.skinT vars.skinT; vars.skinT vars.skinT vars.stim_cold vars.stim_cold vars.skinT; vars.skinT vars.skinT vars.skinT vars.stim_cold vars.stim_cold];
            coldTempsIndex = randi(size(allColdTemps,1)); % choose randomly among a predefined set of cold temperatures
            Results.targetT(thisTrial,1:5) =  allColdTemps(coldTempsIndex,:);
        elseif thisTrialStim == 2                % TGI
            %tgiCounter = tgiCounter + 1;
            Results.targetT(thisTrial,1:5) = [vars.stim_warm vars.stim_cold vars.stim_warm vars.stim_cold vars.stim_warm]; % unique stimulus configuration
        end
            
        %% Stimulation ON
        if ~vars.stimFlag 
            disp('debugging without stimulation')
        else
            disp(['Stimulation: ', tempPresented])
            Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
            [~, stimOn] = Screen('Flip',scr.win);
            Results.SOT_stimOn(thisTrial) = stimOn - Results.SessionStartT;
            
            % Stimulate
            [~, Results, ~, vars] = TrialPulse_TGI(keys, Results, scr, vars, thisTrial);
        end
        

        %% Stimulation OFF
        if vars.stimFlag
            TcsAbortStimulation(vars.ser) % stop stimulation
            Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
            [~, stimOff] = Screen('Flip',scr.win);
            Results.SOT_stimOff(thisTrial) = stimOff - Results.SessionStartT;
        end
        
        %% Skip ratings and turn stimulation OFF 
        if ~vars.ratings_yn(thisTrial) % no ratings trial
            % update results as Nan, because no ratings were required in
            % this trial
            Results.SOT_vasOn(thisTrial,1:vars.nRatings) = NaN;
            Results.SOT_vasOff(thisTrial,1:vars.nRatings) = NaN;
            Results.vasResp(thisTrial,1:vars.nRatings) = NaN;
            Results.vasRT(thisTrial,1:vars.nRatings) = NaN;
            Results.trialSuccess(thisTrial,1) = vars.ValidTrial(1);
        end

        %% Ratings
        Results.vasYN = vars.ratings_yn(thisTrial);
        if vars.ratings_yn(thisTrial) % ratings only for a subset of trials
            
            for thisRating = 1:vars.nRatings
                if useEyeLink
                    % EyeLink:  conf rating
                    startStimText = ['Trial ' num2str(thisTrial) ' confidence screen on'];
                    Eyelink('message', startStimText);
                end

                % Fetch the participant's ratings
                [vars] = getVasRatings(keys, scr, vars, thisRating);
                Results.SOT_vasOn(thisTrial,thisRating)     = vars.vasOnset - Results.SessionStartT;
                Results.SOT_vasOff(thisTrial,thisRating)    = vars.vasOffset - Results.SessionStartT;
               
                % If this trial was successfull, move on...
                %if(vars.ValidTrial(2)), WaitSecs(0.2); end
                % Write trial result to file
                Results.vasResp(thisTrial,thisRating) = vars.vasResp;
                Results.vasRT(thisTrial,thisRating)   = vars.vasRT;

                % Was this a successfull trial? (both emotion and confidence rating valid)
                % 1-success, 0-fail
                Results.trialSuccess(thisTrial,1) = logical(sum(vars.ValidTrial) == 4);
            end
                
        end
        
        Results.stimDuration = Results.SOT_stimOff(thisTrial) - Results.SOT_stimOn(thisTrial); 
        
        %% change in location
        if any(vars.changelocation == thisTrial)   
            % Show message "move arm/thermode"
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            locN = num2str(vars.locationNumber(thisTrial));
            DrawFormattedText(scr.win, uint8([vars.InstructionMove, ' ', locN]), 'center', 'center', scr.TextColour);
            [~, ~] = Screen('Flip', scr.win);
            WaitSecs(3);
        end
        %% ITI / fixation point 
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        if vars.fixCrossFlag
            scr = drawFixation(scr);
        end
        [~, StartIti] = Screen('Flip', scr.win);
        
        Results.SOT_iti(thisTrial) = StartIti - Results.SessionStartT;
        Results.itiDuration = vars.ITI(thisTrial);
        if useEyeLink
            % EyeLink:  ITI
            startStimText = ['Trial ' num2str(thisTrial) ' Fixation point start'];
            Eyelink('message', startStimText);
        end
        
        % Present the gray screen for ITI duration
        while (GetSecs - StartFixPoint) <= vars.ITI(thisTrial)
            
            if keys.KeyCode(keys.Escape)==1
                % Save, mark the run
                vars.RunSuccessfull = 0;
                vars.Aborted = 1;
                experimentEnd(keys, Results, scr, vars);
                return
            end
        end
        
        [~, ~, keys.KeyCode] = KbCheck;
        WaitSecs(0.001);
          
        %% Clean up
        Results.trialValidity(thisTrial,1:4) = vars.ValidTrial;
        Results.trialDuration(thisTrial) = GetSecs - Results.SOT_trial(thisTrial);
        vars.predResp = NaN;            % reset cold/warm prediction
        thisTrialCorrect = NaN;         % reset prediction accuracy
        vars.vasResp = NaN;             % reset vas responses
        vars.ValidTrial = zeros(1,4);   % reset trial validity

        if useEyeLink
            % EyeLink:  trial end
            startStimText = ['Trial ' num2str(thisTrial) ' end'];
            Eyelink('message', startStimText);
        end
        
        %% save data at every trial
        save(strcat(vars.OutputFolder, vars.UniqueFileName), 'Results', 'vars', 'scr', 'keys');
        
        %% Finish cleaning up after the trial
         % If the trial was missed, repeat it or go on...
        if vars.RepeatMissedTrials
            % if this was a valid trial, advance one. Else, repeat it.
            if vars.ValidTrial(1)            
                thisTrial = thisTrial + 1;
            else
                disp('Invalid response. Repeating trial.');
                % Repeat the trial...
            end
        else
            % Advance one trial (always in MR)
            thisTrial = thisTrial + 1;
        end
        
        %% Should we have a break here?
        if any(vars.breaks == thisTrial)
            switch vars.breaksSelfPace
                case 0 % fixed timing (no self pace)
                    % Gray screen - Take a short break
                    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                    DrawFormattedText(scr.win, uint8(vars.InstructionPause{1}), 'center', 'center', scr.TextColour);
                    [~, breakStartsNow] = Screen('Flip', scr.win);

                    % wait to continue by pressing the space bar

                    % wait for vars.breakT seconds
                    while (GetSecs - breakStartsNow) <= vars.breakT
                        % Draw time remaining on the screen
                        breakRemaining = vars.breakT - (GetSecs - breakStartsNow);
                        breakRemainingString = [num2str(round(breakRemaining)), ' seconds'];
                        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                        DrawFormattedText(scr.win, uint8(vars.InstructionPause), 'center', 'center', scr.TextColour);
                        DrawFormattedText(scr.win, uint8(breakRemainingString), 'center', ((scr.winRect(4)/2)+200), scr.TextColour);
                        [~, ~] = Screen('Flip', scr.win);
                        WaitSecs(1); 
                    end

                    % Get ready to continue...
                    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                    DrawFormattedText(scr.win, 'The experiment will now continue, get ready.', 'center', 'center', scr.TextColour);
                    [~, ~] = Screen('Flip', scr.win);
                    WaitSecs(3);
                    
               case 1 % self pace
                    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
                    DrawFormattedText(scr.win, uint8([vars.InstructionPause{2}]), 'center', 'center', scr.TextColour, 60);
                    [~, ~] = Screen('Flip', scr.win);

                    new_line;
                    disp(['Press SPACE to continue. Restarting from trial ' num2str(thisTrial)]); new_line;

                    while keys.KeyCode(keys.Space) == 0 % Wait for trigger
                        [~, ~, keys.KeyCode] = KbCheck;
                        WaitSecs(0.001);
                    end
            end
        end
    end%thisTrial
        
    vars.RunSuccessfull = 1;
    Results.SessionEndT = GetSecs - Results.SessionStartT; % session end

    % Save, mark the run
    experimentEnd(keys, Results, scr, vars);
    
    %% EyeLink: experiment end
    if useEyeLink
        ELshutdown(vars)
    end
    
    %% Draw plux trigger -- END
    if vars.pluxSynch
        Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
        [~, ~] = Screen('Flip', scr.win);
        WaitSecs(2)
    end
    [~, ~] = Screen('Flip', scr.win);
    
    %%
    ShowCursor;
    
catch ME% Error. Clean up...
    
    % Save, mark the run
    
    vars.RunSuccessfull = 0;
    vars.Error = 1;
    vars.Aborted = 0;
    experimentEnd(keys, Results, scr, vars);
    rethrow(ME)
end
