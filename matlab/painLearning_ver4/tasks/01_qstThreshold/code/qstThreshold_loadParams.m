%% Define parameters
%
% Project: Implementation of QST thermal thresholds
%
% Francesca Fardo 
% Last edit: 03/07/2021

%% Key flags
vars.inputDevice    = 1;   % Response method for button presses 1 - mouse, %2 - qstLab button press device 
useEyeLink          = 0;  

%% Fixed stimulation parameters
vars.timepoint          = 0;   % timer 
vars.tcsData            = [];  % variable to store stimulation data   
%vars.skinT              = 32;  % baseline temperature
vars.speed_ramp         = 1;   % rate of temperature change from baseline to target
vars.return_ramp        = 1;   % rate of temperature change from target to baseline
vars.durations          = 32;  % max stimulation duration
vars.bp                 = 1;   % 1 = button press required
vars.save_file          = 1;   % 1 = save tcs data

%% Stimulation parameters
vars.nQstMeasures       = 5;
vars.nTrials            = 4;
vars.NTrialsTotal       = vars.nQstMeasures*vars.nTrials;
vars.qstMeasureName     = repelem({'cold detection threshold', ...
                           'warm detection threshold', ...
                           'cold pain threshold', ...
                           'heat pain threshold', ...
                           'tgi threshold'},1,vars.nTrials);
vars.baseline           = repelem([32 32 32 32 30],1,vars.nTrials);
vars.targetThreshold    = [repelem(zeros(1,5),vars.nTrials,1); ...
                           repelem(repmat(50,1,5),vars.nTrials,1); ...
                           repelem(zeros(1,5),vars.nTrials,1); ...
                           repelem(repmat(50,1,5),vars.nTrials,1); ...
                           repelem([50 0 50 0 50],vars.nTrials,1)];
vars.thresholdN         = repelem(1:vars.nQstMeasures,1,vars.nTrials);
vars.showInstruction    = sort([1:vars.nTrials:vars.NTrialsTotal 2:vars.nTrials:vars.NTrialsTotal]);

%% Task parameters
vars.showInstruction    = sort([1:vars.nTrials:vars.NTrialsTotal 2:vars.nTrials:vars.NTrialsTotal]); % present instructions for specific trials
vars.jitter             = randInRange(3,1,vars.NTrialsTotal); % time between the beginning of the trial and the beginning of the stimulation
vars.feedbackBPtime     = ones(1,vars.NTrialsTotal); % this determines how long the feedback "button press detected" is shown on the screen
vars.ITI                = 5 - (vars.jitter + vars.feedbackBPtime);
vars.ITI(1:4:12)        = 0; % remove ITI after each demo trial (as the next trial requires intructions and button press to start)
vars.ITI(end)           = 0; % remove last ITI (to move to the next task faster) 

%% Instructions
textSize = 35;
switch vars.language
    case 1 % English

        vars.InstructionPause = 'Take a short break...'; 
        vars.InstructionNext = 'Please press SPACE to continue.';
        vars.InstructionEnd = 'You have completed the first task!';
        vars.FeedbackBP1 = 'Button press detected';
        vars.FeedbackBP = vars.FeedbackBP1;
        vars.FeedbackBP2 = 'Trial ended. Threshold is outside the tested temperature range.';

        vars.InstructionTask =  {'Cold detection (demo) \n \n Please position the thermode to location 1. \n \n First we will test your ability to perceive cold sensations. Please press the left mouse button immediately as soon as you perceive a change in temperature to cool/cooler for the first time. \n \n \n \n  Press SPACE to start.', ...
                     'Cold detection \n \n Please press the left mouse button immediately as soon as you perceive a change in temperature to cool/cooler for the first time. \n \n This will be reapeated 3 times. \n \n \n \n  Press SPACE to continue.', ...
                     'Warm detection (demo) \n \n Now we will test your ability to perceive warm sensations. Please press the left mouse button immediately as soon as you perceive a change in temperature to warm/warmer for the fist time. \n \n \n \n  Press SPACE to start.', ...
                     'Warm detection \n \n Please press the left mouse button immediately as soon as you perceive a change in temperature to warm/warmer for the fist time. \n \n This will be repeated 3 times. \n \n \n \n  Press SPACE to continue.', ...
                     'Pain threshold 1 (demo) \n \n Now we will test when you perceive painfully cold sensations. Your skin will be slowly cooled. At some point in time you will feel a second sensation on top of the cold sensation. The impression of cold will change its quality towards an additional impression of a burning, stinging, drilling or aching sensation. Please press the left mouse button as soon as you perceive such a change. Please do not wait until the sensation has become unbearably painful. \n \n \n \n  Press SPACE to continue.', ... 
                     'Pain threshold 1 \n \n  Please press the left mouse button as soon as the cold sensation changes its quality to an additional burning, stinging, drilling or aching sensation. \n \n This will be repeated 3 times. \n \n \n \n  Press SPACE to continue.', ... 
                     'Pain threshold 2 (demo) \n \n Now we will test when you perceive painfully warm or heat sensations. Your skin will be slowly warmed. At some point in time you will feel a second sensation on top of the warm or heat sensation. The impression of warm/heat will change its quality towards an additional impression of a burning, stinging, drilling or aching sensation. Please press the left mouse button as soon as you perceive such a change. Please do not wait until the sensation has become unbearably painful. \n \n \n \n  Press SPACE to continue.',...
                     'Pain threshold 2 \n \n Please press the left mouse button as soon as the warm or heat sensation changes its quality to an additional burning, stinging, drilling or aching sensation. \n \n This will be repeated 3 times. \n \n \n \n Press SPACE to continue.',...
                     'Pain threshold 3 (demo) \n \n Now we will test when you perceive painful sensations again. Your skin will be slowly cooled or warmed. At some point in time you will feel a second sensation on top of a cold, warm or heat sensation. The impression of cold, warmth or heat will change its quality towards an additional impression of a burning, stinging, drilling or aching sensation. Please press the left mouse button as soon as you perceive such a change. Please do not wait until the sensation has become unbearably painful. \n \n \n \n  Press SPACE to continue.',...
                     'Pain threshold 3 \n \n Please press the left mouse button as soon as the cold or heat sensation changes its quality to an additional burning, stinging, drilling or aching sensation. \n \n This will be repeated 3 times. \n \n \n \n Press SPACE to continue.'};        
    case 2 % Danish

        vars.InstructionPause = 'Tag en kort pause...'; % change this one
        vars.InstructionNext = 'Tryk MELLEMRUM for at fortsætte.'; % change this one
        vars.InstructionEnd = 'Opgaven er nu slut. Tak!';
        vars.FeedbackBP = 'Der blev trykket på knappen.';
        vars.FeedbackBP1 = vars.FeedbackBP;
        vars.FeedbackBP2 = 'Runden er slut. Grænseværdien er udenfor den testede temperaturramme.';
 
        vars.InstructionTask =  {'Kuldeopfattelse (demo) \n \n Placer venligst sensoren på lokation 1.  \n \n Først vil vi teste din evne til at føle kulde. Vær venlig at klikke på den venstre museknap ligeså snart du føler temperaturen skifte til kold eller koldere for første gang. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ...
                     'Kuldeopfattelse \n \n Vær venlig at klikke på den venstre museknap ligeså snart du føler temperaturen skifte til kold eller koldere for første gang. \n \n Dette måles i alt 3 gange. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ...
                     'Varmeopfattelse (demo) \n \n Vi vil nu teste din evne til at føle varme. Vær venlig at klikke på den venstre museknap ligeså snart du føler temperaturen skifte til varm eller varmere for første gang. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ...
                     'Varmeopfattelse \n \n Vær venlig at klikke på den venstre museknap ligeså snart du føler temperaturen skifte til varm eller varmere for første gang.  \n \n Dette måles i alt 3 gange. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ...
                     'Smerte-grænseværdi 1 (demo) \n \n Vi vil nu teste, hvornår du føler kulden som smertefuld mod din hud. Sensoren vil langsomt køle din hud. Kulden vil på et tidspunkt skifte fra blot at føles kold til yderligere at føles som brændende, sviende, borende eller smertende. Så snart du mærker et sådant skift, skal du klikke på den venstre musetast. Du skal ikke vente med at klikke, til smerten bliver uudholdelig. Herefter vil sensoren varmes op igen, indtil den når udgangstemperaturen. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ... 
                     'Smerte-grænseværdi 1 \n \n Vær venlig at klikke på den venstre museknap, så snart kuldefornemmelsen skifter til yderligere at føles brændende, sviende, borende eller smertende. \n \n Dette måles i alt 3 gange. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ... 
                     'Smerte-grænseværdi 2 (demo) \n \n Nu vil vi ligeledes teste, hvornår du føler varmen som smertefuld mod din hud. Sensoren vil langsomt varme din hud. Varmen vil på et tidspunkt skifte fra blot at føles varm til også at føles som brændende, sviende, borende eller smertende. Så snart du mærker et sådant skift, skal du klikke på den venstre musetast. Du skal ikke vente med at bruge din afbryder, til smerten bliver uudholdelig. Herefter vil sensoren køle ned igen, indtil den når udgangstemperaturen. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.',...
                     'Smerte-grænseværdi 2\n \n Vær venlig at klikke på den venstre musetast, så snart varmefornemmelsen skifter til yderligere at føles brændende, sviende, borende eller smertende. \n \n Dette måles i alt 3 gange. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.',...
                     'Smerte-grænseværdi 3 (demo) \n \n Vi vil nu teste hvornår du bemærker en smerte fornemmelse igen. Din hud vil langsomt blive kølet eller varmet. På et tidspunkt vil du bemærke en anden fornemmelse ovenpå en kold eller varm fornemmelse. Varme eller kuldefornemmelsen vil ændre sig mod en anden fornemmelse som er brændende, sviende, borende eller smertende. Klik på den venstre musetast lige så snart du bemærker sådan en ændring. Du skal ikke vente med at bruge din afbryder, til smerten bliver uudholdelig. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.', ...
                     'Smerte-grænseværdi 3 \n \n Tryk venligst på den venstre musetast så snart den kolde eller varme fornemmelse ændre sig til en brændende, sviende, borende eller smertende fornemmelse. \n \n Dette måles i alt 3 gange. \n \n \n \n  Tryk på MELLEMRUMSTASTEN for at fortsætte.'};
end