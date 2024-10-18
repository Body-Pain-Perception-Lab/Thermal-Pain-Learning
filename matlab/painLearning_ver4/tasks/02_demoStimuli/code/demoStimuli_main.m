function demoStimuli_main(vars, scr)
%
% Project: Demonstration of the stimuli uses in the pain learning task
%
% Francesca Fardo 
% Last edit: 03/07/2021

%% Load stimulation and task parameters
demoStimuli_loadParams;

%% Define Results struct
uniqueFilename = strcat(vars.OutputFolder, vars.UniqueFileName, '.mat');
if ~exist(uniqueFilename)
    DummyDouble = ones(vars.NTrialsTotal,1).*NaN;
    %DummyString = strings(vars.NTrialsTotal,1);
    Results = struct('SubID',           {DummyDouble}, ...
                     'trialN',          {DummyDouble}, ...
                     'baseline',        {DummyDouble}, ...
                     'targetT',         {DummyDouble}, ...
                     'SessionStartT',   {DummyDouble}, ...
                     'SessionEndT',     {DummyDouble});
else
    vars.confirmedSubjN = input('Subject already exists. Do you want to continue anyway (yes = 1, no = 0)?    ');
    if vars.confirmedSubjN
        load(uniqueFilename,'Results');
        vars.startTrialN = input('Define the trial number to restart from?   ');
        vars.date_time = datestr(now,'ddmmyyyy_HHMMSS');
        vars.DataFileName = strcat(vars.exptName, '_',vars.subIDstring, '_', vars.date_time);    % name of data file to write to
    else
        return
    end
end
%% Keyboard & keys configuration
[keys] = keyConfig();

% Reseed the random-number generator
SetupRand;

%% Prepare to start
 try
    %% Check if window is already open (if not, open screen window) 
      if ~exist('scr','var')
        if ~isfield(scr, 'win')
            % Diplay configuration
            [scr] = displayConfig(scr);
            scr.bkColor = scr.BackgroundGray;
            AssertOpenGL;
            [scr.win, scr.winRect] = PsychImaging('OpenWindow', scr.screenID, scr.BackgroundGray); %,[0 0 1920 1080] mr screen dim
            PsychColorCorrection('SetEncodingGamma', scr.win, 1/scr.GammaGuess);
            
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
            scr.dist = scr.ViewDist;
            scr.width  = scr.MonitorWidth;
            scr.resolution = scr.winRect(3);                    % number of pixels of display in horizontal direction
        end
      end
    
    %% Dummy calls to prevent delays
    vars.ValidTrial = zeros(1,2);
    vars.RunSuccessfull = 0;
    vars.Aborted = 0;
    vars.Error = 0;
    WaitSecs(0.1);
    GetSecs;
    vars.Resp = 888;
    vars.ConfResp = 888;
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
    
    %% Start session
    Results.SessionStartT = GetSecs;            % session start = trigger 1 + dummy vols
    if useEyeLink
        Eyelink('message','STARTEXP');
    end
    
    %% Run through trials
    WaitSecs(0.500);            % pause before experiment start
    thisTrial = vars.startTrialN; % trial counter (user defined)
    endOfExpt = 0;

    while endOfExpt ~= 1       % General stop flag for the loop
        thisTrialStim  = input('Which stimulation, 0 = cold, 1 = warm, 2 = tgi?     ');
        vars.ratings_yn(thisTrial) = 1;%input ('Ratings, yes = 1, no = 0?     ');
        if thisTrialStim == 0
            vars.stim_cold = input('Define cold temperature:    '); 
            vars.stimulation{thisTrial} = 'cold';
            %coldCounter = coldCounter + 1;
            allColdTemps = [vars.stim_cold vars.stim_cold vars.skinT vars.skinT vars.skinT; vars.skinT vars.stim_cold vars.stim_cold vars.skinT vars.skinT; vars.skinT vars.skinT vars.stim_cold vars.stim_cold vars.skinT; vars.skinT vars.skinT vars.skinT vars.stim_cold vars.stim_cold];
            coldTempsIndex = randi(size(allColdTemps,1)); % choose randomly among a predefined set of cold temperatures
            Results.targetT(thisTrial,1:5) =  allColdTemps(coldTempsIndex,:);
        elseif thisTrialStim == 1 % || vars.tempQual == 2
            vars.stim_warm = input('Define warm temperature:    ');
            vars.stimulation{thisTrial} = 'warm';
            %warmCounter = warmCounter + 1;
            allWarmTemps = [vars.stim_warm vars.stim_warm vars.stim_warm vars.skinT vars.skinT; vars.skinT vars.stim_warm vars.stim_warm vars.stim_warm vars.skinT; vars.skinT vars.skinT vars.stim_warm vars.stim_warm vars.stim_warm];
            warmTempsIndex = randi(size(allWarmTemps,1)); % choose randomly among a predefined set of warm temperatures
            Results.targetT(thisTrial,1:5) = allWarmTemps(warmTempsIndex,:); % Define stimulation temperature for this trial
        elseif thisTrialStim == 2
            vars.stim_cold = input('Define cold temperature:    '); 
            vars.stim_warm = input('Define warm temperature:    ');
            vars.stimulation{thisTrial} = 'tgi';
            %tgiCounter = tgiCounter + 1;
            Results.targetT(thisTrial,1:5) = [vars.stim_warm vars.stim_cold vars.stim_warm vars.stim_cold vars.stim_warm]; % unique stimulus configuration
        end
        
        %% Show instructions
        if thisTrial == 1
            % instruction stimuli
            showInstruction(scr, keys, vars.InstructionDemo{1})
            
            % instruction rating
            showInstruction(scr, keys, vars.InstructionDemo{2})
        end
            
        if useEyeLink
            % EyeLink:  this trial
            startStimText = ['Trial ' num2str(thisTrial) ' starts now'];
            Eyelink('message', startStimText); % Send message
        end

        %% Trial starts: draw fixation point
        Results.SOT_trial(thisTrial) = GetSecs - Results.SessionStartT; % trial starts
        disp(['Trial # ', num2str(thisTrial), ' ', vars.stimulation{thisTrial}]);
        
        % Draw Fixation
        [~, ~] = Screen('Flip', scr.win);            % clear screen
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        scr = drawFixation(scr); % fixation point
        [~, ~] = Screen('Flip', scr.win);
        
        %% Jitter
        vars.jitter(thisTrial) = randInRange(3,1,1); % time between the beginning of the trial and the beginning of the stimulation
        WaitSecs(vars.jitter(thisTrial));
            
        %% Stimulation ON
        if ~vars.stimFlag 
            disp('debugging without stimulation')
        else
            disp(['Stimulation: ', vars.stimulation{thisTrial}])
            Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
            [~, stimOn] = Screen('Flip',scr.win);
            Results.SOT_stimOn(thisTrial) = stimOn - Results.SessionStartT;
            
%             if vars.pluxSynch % not needed during demo
%                 % if were in the first pluxDurationSec seconds, draw the rectangle
%                 if thisTrialStim == 0     &&((GetSecs - StimOn) <= scr.pluxDurSec(2)) % Cold
%                     Screen('FillRect', scr.win, scr.pluxBlack, scr.pluxRect);
%                 else % Warm or TGI
%                     Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
%                 end
%             end
%             
            % Stimulate
            [~, Results, ~, vars] = TrialPulse_TGI(keys, Results, scr, vars, thisTrial);
        end
         
        %% Stimulation Off
        Results.SOT_stimOff(thisTrial) = GetSecs;
        
        if useEyeLink
            % EyeLink:  stimulation off
            startStimText = ['Trial ' num2str(thisTrial) ' stim off'];
            Eyelink('message', startStimText); % Send message
        end
        
%% Skip ratings and turn stimulation OFF 
        if ~vars.ratings_yn(thisTrial) % no ratings trial
            % stop stimulation
            if vars.stimFlag 
                TcsAbortStimulation(vars.ser) % stop stimulation
            end
            Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
            [~, stimOff] = Screen('Flip',scr.win);
            Results.SOT_stimOff(thisTrial) = stimOff - Results.SessionStartT;
            % update results as Nan, because no ratings were required in
            % this trial
%             Results.SOT_vas(thisTrial,1:vars.nRatings) = NaN;
%             Results.SOT_vasOn(thisTrial,1:vars.nRatings) = NaN;
%             Results.SOT_vasOff(thisTrial,1:vars.nRatings) = NaN;
            Results.vasResp(thisTrial,1:vars.nRatings) = NaN;
            Results.vasRT(thisTrial,1:vars.nRatings) = NaN;
            %Results.trialSuccess(thisTrial,1:vars.nRatings) = NaN;
            Results.trialSuccess(thisTrial,1:vars.nRatings) = logical(sum(vars.ValidTrial) == 1);
        end

        %% Ratings and then turn stimulation OFF
        if vars.ratings_yn(thisTrial) % ratings only for a subset of trials
            for thisRating = 1:vars.nRatings
                if useEyeLink
                    % EyeLink:  conf rating
                    startStimText = ['Trial ' num2str(thisTrial) ' confidence screen on'];
                    Eyelink('message', startStimText);
                end

                % Fetch the participant's ratings
                [vars] = getVasRatings(keys, scr, vars, thisRating);
                %Results.SOT_vasTime(thisTrial,thisRating)   = vars.vasTime - Results.SessionStartT;
                %Results.SOT_vasOn(thisTrial,thisRating)     = vars.vasOnset - Results.SessionStartT;
                %Results.SOT_vasOff(thisTrial,thisRating)    = vars.vasOffset - Results.SessionStartT;
               
                % If this trial was successfull, move on...
                %if(vars.ValidTrial(2)), WaitSecs(0.2); end
                % Write trial result to file
                Results.vasResp(thisTrial,thisRating) = vars.vasResp;
                Results.vasRT(thisTrial,thisRating)   = vars.vasRT;

                % Was this a successfull trial? (both emotion and confidence rating valid)
                % 1-success, 0-fail
                Results.trialSuccess(thisTrial,thisRating) = logical(sum(vars.ValidTrial) == 2);
            end
            
            % Stimulation OFF
            if vars.stimFlag
                TcsAbortStimulation(vars.ser) % stop stimulation
                Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
                [~, stimOff] = Screen('Flip',scr.win);
                Results.SOT_stimOff(thisTrial) = stimOff - Results.SessionStartT;
            end
        end

        %% Update Results
        Results.SubID(thisTrial,1)        = vars.subNo;
        Results.trialN(thisTrial,1)       = thisTrial;
        Results.baseline(thisTrial,1)     = vars.skinT;
        
        %% save data at every trial
        %save(strcat(vars.OutputFolder, vars.UniqueFileName), 'Results', 'vars', 'scr', 'keys' );
        save(strcat(vars.OutputFolder, vars.UniqueFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );

        %% ITI / prepare for next trial
        % Draw Fixation again
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        %scr = drawFixation(scr); % fixation point
        [~, StartITI] = Screen('Flip', scr.win);
        
        Results.SOT_ITI(thisTrial) = GetSecs;
        if useEyeLink
            % EyeLink:  ITI
            startStimText = ['Trial ' num2str(thisTrial) ' ITI start'];
            Eyelink('message', startStimText); % Send message
        end
        
        % Define ITI
        vars.ITI(thisTrial) = 0;
        
        % Present the gray screen for ITI duration
        while (GetSecs - StartITI) <= vars.ITI(thisTrial)
            
            if keys.KeyCode(keys.Escape)==1
                % Save, mark the run
                vars.RunSuccessfull = 0;
                vars.Aborted = 1;
                experimentEnd(vars, scr, keys, Results);
                return
            end
        end
        
        [~, ~, keys.KeyCode] = KbCheck;
        WaitSecs(0.001);
        
        Results.TrialDuration(thisTrial) = GetSecs - Results.SOT_ITI(thisTrial);
        
        % Reset Texture, ValidTrial, Resp
        vars.ValidTrial = zeros(1,2);
        
        if useEyeLink
            % EyeLink:  trial end
            startStimText = ['Trial ' num2str(thisTrial) ' ends now'];
            Eyelink('message', startStimText);          % Send message
        end
        
        %% Continue to next trial or time to stop? (max # trials reached)
        endOfExpt = input('End of demonstration, yes = 1 or no = 0?     ');
        if ~endOfExpt
            % Advance one trial
            thisTrial = thisTrial + 1;
        end
        
    end % end trial

    Results.SessionEndT = GetSecs  - Results.SessionStartT;
    vars.RunSuccessfull = 1;
    
    % Save, mark the run
    experimentEnd(vars, scr, keys, Results);
    
    %% EyeLink: experiment end
    if useEyeLink
        ELshutdown(vars)
    end
    
    ShowCursor;
    
catch ME% Error. Clean up...
    
    % Save, mark the run
    vars.RunSuccessfull = 0;
    vars.Error = 1;
    experimentEnd(vars, scr, keys, Results);
    rethrow(ME)
end
