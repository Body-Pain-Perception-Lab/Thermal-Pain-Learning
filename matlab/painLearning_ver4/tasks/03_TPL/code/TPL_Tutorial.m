function TPL_Tutorial(scr, keys, vars)
%cwtTutorial(scr, keys, vars)
%
% Runs a tutorial for CWT task
%
%   Input:
%       scr       screen parameters structure
%       keys      key names structure
%       vars      general vars (set by loadParams.m)
%
% Niia Nikolova 25/06/2021


% 1. CWT instructions

% 2. CWT example trial, no confidence

% 3. CWT example trial, confidence

% 4. Learning task explanation

% 5. Prediction trial explanation

% 6. CWT, some example trials

%global language
%vars.language = language;

%% Set variables & instructions
nTrialsTPL      = 10;        % Number of tutorial trials to run
tutorialCues    = round(rand(1,nTrialsTPL));
tutorialStims   = round(rand(1,nTrialsTPL));
tutorialRatings = [0 0 0 0 1 0 0 1 0 1];
tutorialMoveArm = [0 0 0 0 0 0 0 1 0 0];

% Instructions
if vars.language == 1       % English
    instr.A     = '--- Learning task --- \n \n Please position the thermode on location 2. \n \n  On each trial, you will first hear a sound (a tone). Then you will be asked to guess whether the upcoming stimulation will be cold or warm (by pressing the LEFT / RIGHT mouse buttons). Afterwards you will feel a stimulus. \n \n \n \n Press SPACE to continue.';
    instr.B     = 'After some trials, you will rate the intensity of you feelings by clicking on a slider scale, as in the previous demo. You will have max 5 seconds to complete each rating. \n \n \n \n Press SPACE to continue.';
    instr.C     = 'In addition to this, there is a learning component to the task. There is a relationship between the tones and the stimuli in such a way that a given tone predicts whether the following stimulus will be cold or warm. \n \n \n \n Press SPACE to continue. ';
    instr.D     = 'For example, the higher tone might start out predicting a cold stimulation, while the lower tone predicts a warm stimulation. Crucially, these predictive associations will change over the course of the session. So the higher tone may eventually go on to predict warm stimuli, then again cold stimuli, and so on. Note that although the tones predict the thermal quality with some certainty, this is not 100%. This means that there may be some trials that do not go with the current relationship. \n \n \n \n Press SPACE to continue.';
    instr.E     = 'We would like you to try to learn what the associations are at any given time. \n \n Now you will do a few practice trials. \n \n \n \n Press SPACE to start. ';
    %instr.F     = 'Now you will do a few practice trials. \n \n \n \n Press SPACE to continue.';
    instr.G     = 'You have completed the tutorial and will now go on to the main experiment. This will take about 60 minutes. \n \n \n \n You will have several opportunities to take breaks. \n \n \n \n Press SPACE to continue.';
    instr.H     = 'Get ready…';
    
    instr.move = 'Please move the thermode to location 3.';
    
    instr.feedbackA = 'Response: Cold';
    instr.feedbackH = 'Response: Warm';
    instr.feedbackS = 'Too slow!';
    
elseif vars.language == 2       % Danish
    
    instr.A     = '--- 2. Læringsopgave  --- \n \n Placer venligst sensoren på lokation 2. \n \n I hver runde vil du først høre en lyd (en tone). \n \n Så vil du blive bedt om at gætte om den kommende stimulus vil være kold eller varm (ved at trykke på den VENSTRE / HØJRE mussetast). Efterfølgende vil du mærke en stimulus. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte.';
    instr.B     = 'Efter nogle runder skal du angive hvor intens følelsen var ved at trykke på en skala ligesom i den foregående demo. Du har maks 5 sekunder til at svare. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte.';
    instr.C     = 'Derudover er der en læringskomponent i denne opgave. Der er et forhold imellem tonerne og stimulien på en sådan måde at en given tone forudsiger om den efterfølgende stimulus vil være kold eller varm. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte. ';
    instr.D     = 'F.eks, den højere tone kan starte med at forudsige en kold stimulus, mens den lavere tone forudsiger en varm stimulus. Det er vigtigt at vide at disse forudsigende associationer vil ændre sig igennem denne opgave. Så den højere tone kommer måske til at forudse en varm stimulus og så en kold stimulus, og så videre. Bemærk at selvom tonerne forudser stimuliens karakter med nogen sikkerhed, så er det ikke 100%. Dette betyder at der vil være nogle runder som ikke passer med det nuværende forhold. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte';
    instr.E     = 'Vi vil gerne have at du prøver at lærer associationerne på det givne tidspunkt. \n \n Nu vil der være et par prøve runder. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at starte ';
    %instr.F     = 'Du vil nu gennemføre nogle få øve runder. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte.';
    instr.G     = 'Du har nu gennemført øvelsen og vil nu fortsætte til hoveddelen af eksperimentet. Dette tager omkring 60 minutter. \n \n \n \n Du vil have flere muligheder for at tage pauser. \n \n \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte.';
    instr.H     = 'Gør dig klar…';
    
    instr.move = 'Placer venligst sensoren på lokation 3';
    
    instr.feedbackA = 'Svar: Kold';
    instr.feedbackH = 'Svar: Varm';
    instr.feedbackS = 'For langsomt!';

end

global tutorialAbort

try
    
    pause(0.200);
    [~, ~, keys.KeyCode] = KbCheck;
    
    
    clear thisTrial
    
    %% 1. Example trials
    for thisTrial = 1 : (nTrialsTPL) % add two for introducing response & confidence
        
        if thisTrial == 1
            % General task instruction, example trial & response
             showInstruction(scr, keys, instr.A);
            % Ratings
             showInstruction(scr, keys, instr.B);
            % Learning component
             showInstruction(scr, keys, instr.C);
             showInstruction(scr, keys, instr.D);
             showInstruction(scr, keys, instr.E);
        end
        
        %% Cue
        thisCue = tutorialCues(thisTrial);
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

%                     % Draw plux trigger -- CUE
%                     if vars.pluxSynch
%                         % if were in the first pluxDurationSec seconds, draw the rectangle
%                         if thisCue == 0     &&((GetSecs - CueOn) <= scr.pluxDurSec(1)) 
%                             Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
%                         elseif thisCue == 1 &&((GetSecs - CueOn) <= scr.pluxDurSec(1))
%                             Screen('FillRect', scr.win, scr.pluxBlack, scr.pluxRect);
%                         end
%                     end

                    % Flip screen
                    Screen('Flip', scr.win);

                    % KbCheck for Esc key
                    if keys.KeyCode(keys.Escape)==1
                        % set tutorialAbort to 1
                        tutorialAbort = 1;
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
                [~, ~] = Screen('Flip', scr.win);
                %Results.SOT_cueOn(thisTrial) = CueOn - Results.SessionStartT;
                playBeep(cueFreq) % auditory cue
                
        end
        
        % cue off
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        [~, ~] = Screen('Flip', scr.win);
        %Results.SOT_cueOff(thisTrial) = CueOff - Results.SessionStartT;
        
        
        %% ISI (between cue off and display of prediction question)
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        if vars.fixCrossFlag
            scr = drawFixation(scr);
        end
        [~, StartJitter1] = Screen('Flip', scr.win);
        
        % Present the gray screen for ITI duration
        while (GetSecs - StartJitter1) <= vars.jitter1(thisTrial)
            
            if keys.KeyCode(keys.Escape)==1
                % set tutorialAbort to 1
                tutorialAbort = 1;
                return
            end
        end
        
        [~, ~, keys.KeyCode] = KbCheck;
        WaitSecs(0.001);
        
        %% Prediction
        % Will the outcome be C or W?
        thisTrialStim = tutorialStims(thisTrial);         % 1 Warm, 0 Cold, 2 TGI
        if thisTrialStim == 0
            tempPresented = 'cold';
        elseif thisTrialStim == 1
            tempPresented = 'warm';
        elseif thisTrialStim == 2
            tempPresented = 'tgi';
        end
       
        % Cold (L arrow) or Warm (R arrow)?
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        DrawFormattedText(scr.win, [vars.InstructionQ], 'center', 'center', scr.TextColour);
        
        [~, vars.StartRT] = Screen('Flip', scr.win);
        
        % Fetch the participant's response, via keyboard or mouse
        [vars] = getResponse(keys, scr, vars);
        
%         if vars.abortFlag               % Esc was pressed
%             % set tutorialAbort to 1
%             tutorialAbort = 1;
%             return
%         end
%         
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
        disp(['Stimulation: ', tempPresented])
        Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
        [~, ~] = Screen('Flip',scr.win);
        if ~vars.stimFlag 
            disp('debugging without stimulation')
        else

            %Results.SOT_stimOn(thisTrial) = stimOn - Results.SessionStartT;
            
%             if vars.pluxSynch
%                 % if were in the first pluxDurationSec seconds, draw the rectangle
%                 if thisTrialStim == 0    % &&((GetSecs - StimOn) <= scr.pluxDurSec(2)) % Cold
%                     Screen('FillRect', scr.win, scr.pluxBlack, scr.pluxRect);
%                 else % Warm or TGI
%                     Screen('FillRect', scr.win, scr.pluxWhite, scr.pluxRect);
%                 end
%             end
            
            % Stimulate
            [~, Results, ~, vars] = TrialPulse_TGI(keys, Results, scr, vars, thisTrial);
        end
        
        %% Ratings and then turn stimulation OFF
        if tutorialRatings(thisTrial) % ratings only for a subset of trials
            for thisRating = 1:vars.nRatings

                % Fetch the participant's ratings
                [vars] = getVasRatings(keys, scr, vars, thisRating);
                 
            end
        end
            
        %% Stimulation OFF
        if vars.stimFlag
            TcsAbortStimulation(vars.ser) % stop stimulation
            Screen('FillRect',scr.win, scr.BackgroundGray, scr.winRect);
            [~, ~] = Screen('Flip',scr.win);
        end
        
        %% Move arm 
        if tutorialMoveArm(thisTrial) % ratings only for a subset of trials
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, uint8(instr.move), 'center', 'center', scr.TextColour);
            [~, ~] = Screen('Flip', scr.win);
            WaitSecs(3);
        end

        %% ITI / fixation point 
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        if vars.fixCrossFlag
            scr = drawFixation(scr);
        end
        [~, StartFixPoint] = Screen('Flip', scr.win);
        
        % Present the gray screen for ITI duration
        while (GetSecs - StartFixPoint) <= vars.ITI(thisTrial)
            
            if keys.KeyCode(keys.Escape)==1
                % set tutorialAbort to 1
                tutorialAbort = 1;
                return
            end
        end
        
        [~, ~, keys.KeyCode] = KbCheck;
        WaitSecs(0.001);
        
        % Clean up
        vars.Resp = NaN;            % reset H A resp
        
    end
    
    %% Tutoial complete..
     showInstruction(scr, keys, instr.G);
    
catch ME
    rethrow(ME)
    
    
end

end
