%% Define parameters
%
% Project: Demonstration of the stimuli uses in the pain learning task
%
% Francesca Fardo 
% Last edit: 03/07/2021

%% Key flags
vars.expN           = 2;
vars.inputDevice    = 1;   % Response method for button presses 1 - mouse, %2 - qstLab button press device 
useEyeLink          = 0;  

%% Fixed Stimulation parameters
vars.timepoint = 0;         % timer 
vars.tcsData = [];          % variable to store stimulation data   
vars.skinT = 32;            % baseline temperature
vars.speed_ramp = 20;       % rate of temperature change from baseline to target
vars.return_ramp = vars.speed_ramp/2;      % rate of temperature change from target to baseline
vars.bp = 0;                % 1 = button press required
vars.save_file = 1;         % 1 = save tcs data

%% Task timing and parameters
vars.InputDevice        = 2;    % mouse (2)
vars.NTrialsTotal       = 1;    % minimum number of stimuli
vars.RecT               = 5;    % record tcs data for 5 sec, then auditory cue
vars.StimT              = 35;   % stimulate for max 30 sec (but abort stimulation earlier
vars.ConfT              = 10;   % present each rating scale for max X seconds
vars.showInstruction    = 1;    % show instruction only before the first trial
vars.nRatings           = 3;    % number of VAS ratings per trial

%% Instructions
textSize = 35;
switch vars.language
    case 1 % English   
        vars.InstructionDemo{1} = 'Stimuli Demo \n \n Please position the thermode to location 2. \n \n First we will demonstrate the cold and warm stimuli we will use during the main task. Most of the stimuli will be mild cold or warm temperatures. In some trial, the stimuli may be felt more intensely. In all trials the stimulation temperatures are safe, without risk of skin damage. \n \n Press space bar to continue.';
        vars.InstructionDemo{2} = 'Ratings Demo \n \n Sometimes you will be asked to rate how you perceived the stimuli. The following instructions will appear on the screen one at a time, always in the same order: \n \n (1) Please rate the most intense cold sensation you felt? \n \n (2) Please rate the most intense warm sensation you felt \n \n (3)  Please rate the most intense burning sensation you felt \n \n Ready to start with the stimuli & ratings demo?';
        vars.InstructionDemo{3} = 'Please move the indicator along the line and confirm with a left click, as fast and accurately as possible. \n \n If you do not perceive the sensation that is described in the question, make sure to select the extreme left position (rating = 0/not at all).';
        value = {'cold', 'warm', 'burning', 'painful'};
        valueBold = {'COLD', 'WARM', 'BURNING', 'PAINFUL'};
        for ii = 1:4
            vars.InstructionRating{ii} = ['Please rate the most intense ' valueBold{ii} '  you felt.'];
        end
        vars.ConfEndPoins = {'Not at all', 'Extremely'};
        
    case 2 % Danish
        vars.InstructionDemo{1} = 'Stimuli Demo \n \n Placer venligst sensoren på lokation 2. \n \n Vi vil først demonstrere den kolde og varme stimuli, som vi vil bruge i hovedeopgaven. De fleste stimuli vil have milde kolde eller milde varme temperaturer. I nogle runder vil stimulien muligvis føles mere intens. Temperaturerne er ikke farlige i nogle af runderne. Der er ingen risiko for skade på din hud. \n \n Tryk på MELLEMRUMSTASTEN for at fortsætte.';
        vars.InstructionDemo{2} = 'Bedømmelsesdemo \n \n Nogle gange vil du blive spurgt, hvordan stimulien blev fornemmet. De følgende instruktioner vil komme til syne på skærmen et af gangen. Instruktionerne kommer altid i samme rækkefølge. \n \n (1) Bedøm venligst den mest intense kolde fornemmelse du bemærkede \n \n (2) Bedøm venligst den mest intense varme fornemmelse du bemærkede \n \n (3)  Bedøm venligst den mest intense brændende fornemmelse du bemærkede \n \n Er du klar til at starte med stimuli og bedømmelsesdemoen?';
        vars.InstructionDemo{3} = 'Du bedes venligst rykke indikatorlinjen langs skalen og bekræfte med en venstre click, så hurtigt som muligt. \n \n Hvis du ikke oplever den følelse beskrevet i spørgsmålet, skal du vælge positionen helt til venstre (bedømmelse = 0/slet ikke).';
        value = {'kolde', 'varme', 'brændende', 'smertefuld'};
        valueBold = {'KOLD', 'VARM', 'BRÆNDENDE', 'SMERTEFULD'};
        for ii = 1:4
            vars.InstructionRating{ii} = ['Bedøm venligst den mest intense ' value{ii} ' fornemmelse du bemærkede. \n \n' valueBold{ii} '?'];
        end
        vars.ConfEndPoins = {'Overhovedet ikke', 'Ekstremt'};
end
 
