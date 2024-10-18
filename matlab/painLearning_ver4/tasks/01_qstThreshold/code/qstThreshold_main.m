function qstThreshold_main(vars, scr)
%
% Project: Implementation of QST thermal thresholds
%
% Niia Nikolova & Francesca Fardo 
% Last edit: 03/07/2021

%% Load stimulation and task parameters
qstThreshold_loadParams;

%% Define Results struct
%%
uniqueFilename = strcat(vars.OutputFolder, vars.UniqueFileName, '.mat');
if ~exist(uniqueFilename)
    DummyDouble = ones(vars.NTrialsTotal,1).*NaN;
    %DummyString = strings(vars.NTrialsTotal,1);
    Results = struct('SubID',           {DummyDouble}, ...
                     'thresholdN',      {DummyDouble}, ...
                     'trialN',          {DummyDouble}, ...
                     'baseline',        {DummyDouble}, ...
                     'targetT',         {DummyDouble}, ...
                     'SOT_trial',       {DummyDouble}, ...
                     'SOT_jitter',      {DummyDouble}, ...
                     'SOT_stimOn',      {DummyDouble}, ...
                     'SOT_stimOff',     {DummyDouble}, ...
                     'SOT_ITI',         {DummyDouble}, ...
                     'TrialDuration',   {DummyDouble}, ...
                     'tcsData',         {DummyDouble}, ...
                     'SessionStartT',   {DummyDouble}, ...
                     'SessionEndT',     {DummyDouble});
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
    
    %% Plux parameters
    scr.pluxDur     = [2;25];
    scr.hz          = Screen('NominalFrameRate', scr.win); 
    scr.pluxDurSec  =  scr.pluxDur / scr.hz;
    
    %% Dummy calls to prevent delays
    vars.ValidTrial = zeros(1,2);
    vars.RunSuccessfull = 0;
    vars.Aborted = 0;
    vars.Error = 0;
    WaitSecs(0.1);
    GetSecs;
    vars.Resp = NaN;
    vars.ConfResp = NaN;
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
    WaitSecs(0.5);            % pause before experiment start
    thisTrial = vars.startTrialN; % trial counter (user defined)
    endOfExpt = 0;

    while endOfExpt ~= 1       % General stop flag for the loop
        Results.SOT_trial(thisTrial) = GetSecs - Results.SessionStartT; % trial starts
        %% show instructions
        if any(vars.showInstruction == thisTrial) 
            whichInstruction = find(vars.showInstruction == thisTrial);
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, uint8([vars.InstructionTask{whichInstruction}]), 'center', 'center', scr.TextColour, 60); %#ok<FNDSB>
            [~, ~] = Screen('Flip', scr.win);

            new_line;
            disp(['Press SPACE to start ' vars.qstMeasureName{thisTrial}]); new_line;

            while keys.KeyCode(keys.Space) == 0 % Wait for trigger
                 [~, ~, keys.KeyCode] = KbCheck;
                WaitSecs(0.001);
                
                if keys.KeyCode(keys.Escape)==1 % if ESC, quit the experiment
                    % Save, mark the run
                    vars.RunSuccessfull = 0;
                    vars.Aborted = 1;
                    experimentEnd(vars, scr, keys, Results)
                    return
                end
            end
        end
        
        if useEyeLink
            % EyeLink:  this trial
            startStimText = ['Trial ' num2str(thisTrial) ' starts now'];
            Eyelink('message', startStimText); % Send message
        end

        %% Trial starts: draw fixation point
        disp(['Trial # ', num2str(thisTrial), ' ', vars.qstMeasureName{thisTrial}]);
        
        % Draw Fixation
        [~, ~] = Screen('Flip', scr.win);            % clear screen
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        scr = drawFixation(scr); % fixation point
        [~, ~] = Screen('Flip', scr.win);
        
        %% Jitter
        WaitSecs(vars.jitter(thisTrial));
        
        %% Stimulation On
        Results.SOT_stimOn(thisTrial) = GetSecs;
        
        if useEyeLink
            % EyeLink:  fixation on
            startStimText = ['Trial ' num2str(thisTrial) ' stim on'];
            Eyelink('message', startStimText); % Send message
        end
        
        if ~vars.stimFlag 
            disp('debugging without stimulation')
            vars.thresholdT = NaN;
            vars.thresholdM = NaN;
        else
            [Results, vars] = TrialQstThresholds(Results, vars, thisTrial); % Stimulate
        end  
        %% Stimulation Off
        Results.SOT_stimOff(thisTrial) = GetSecs;
        
        if useEyeLink
            % EyeLink:  stimulation off
            startStimText = ['Trial ' num2str(thisTrial) ' stim off'];
            Eyelink('message', startStimText); % Send message
        end
        
        %% Provide feedback about button press
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        DrawFormattedText(scr.win, uint8(vars.FeedbackBP), 'center', 'center', scr.TextColour);
        [~, ~] = Screen('Flip', scr.win);
        WaitSecs(vars.feedbackBPtime(thisTrial));
        
        %% Update Results
        Results.SubID(thisTrial,1)        = vars.subNo;
        Results.thresholdN(thisTrial,1)   = vars.thresholdN(thisTrial);
        Results.trialN(thisTrial,1)       = thisTrial;
        Results.baseline(thisTrial,1)     = vars.baseline(thisTrial);
        Results.targetT(thisTrial,1)      = vars.targetThreshold(thisTrial);
        Results.SOT_jitter(thisTrial,1)   = vars.jitter(thisTrial);   
        Results.SOT_ITI(thisTrial,1)      = vars.ITI(thisTrial);
        Results.ThresholdT(thisTrial,1:5) = vars.thresholdT; % all 5 zones
        Results.ThresholdM(thisTrial,1)   = vars.thresholdM; % average of several zones
        
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
        
        Results.TrialDuration(thisTrial) = GetSecs - Results.SOT_jitter(thisTrial);
        
        % Reset Texture, ValidTrial, Resp
        vars.ValidTrial = zeros(1,2);
        
        if useEyeLink
            % EyeLink:  trial end
            startStimText = ['Trial ' num2str(thisTrial) ' ends now'];
            Eyelink('message', startStimText);          % Send message
        end
        
        %% Continue to next trial or time to stop? (max # trials reached)
        if (thisTrial == vars.NTrialsTotal)
            endOfExpt = 1;
        else
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
    
    %% Show summary of results
    CDT = mean(Results.ThresholdM(2:4));
    WDT = mean(Results.ThresholdM(6:8));
    CPT = mean(Results.ThresholdM(10:12));
    HPT = mean(Results.ThresholdM(14:16));
    TGI = mean(Results.ThresholdM(18:20));
    disp(['Threshold estimate for cold detetction: ', num2str(CDT)]);
    disp(['Threshold estimate for warm detection: ', num2str(WDT)])
    disp(['Threshold estimate for cold pain: ', num2str(CPT)]);
    disp(['Threshold estimate for heat pain: ', num2str(HPT)]);
    disp(['Threshold estimate for tgi: ', num2str(round(TGI,2)), ' and ', num2str(round(30-(TGI-30),2))]);
    
catch ME% Error. Clean up...
    
    % Save, mark the run
    vars.RunSuccessfull = 0;
    vars.Error = 1;
    experimentEnd(vars, scr, keys, Results);
    rethrow(ME)
end
