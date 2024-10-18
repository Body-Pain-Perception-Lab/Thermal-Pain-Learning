function experimentEnd(keys, Results, scr, vars)
%function experimentEnd(keys, Results, scr, vars)
%
% Project: CWT task
%
% End of experiment routine. Shows a message to let the user know if the
% run has been aborted or crashed, saves results, and cleans up

%
% Niia Nikolova
% Last edit: 20/07/2020

if isfield(vars,'Aborted') || isfield(vars,'Error')
    if vars.Aborted
        if isfield(scr, 'win')          % if a window is open, display a brief message
            % Abort screen
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, 'Experiment aborted. Exiting...', 'center', 'center', scr.TextColour);
            [~, ~] = Screen('Flip', scr.win);
            %         WaitSecs(0.5);
        end
        
        %     ListenChar(0);
        ShowCursor;
        %     sca;
        disp('Experiment aborted by user!');
        
        % Save, mark the run
        vars.DataFileName = ['Aborted_', vars.DataFileName];
        save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );        
        %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
        disp(['Run was aborted. Results were saved as: ', vars.DataFileName]);
        
        % and as .csv
        csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
        struct2csv(Results, csvName);
        
    elseif vars.Error
        if isfield(scr, 'win')          % if a window is open, display a brief message
            % Abort screen
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, 'Error. Exiting... ', 'center', 'center', scr.TextColour);
            [~, ~] = Screen('Flip', scr.win);
            WaitSecs(0.5);
        end
        
        % Error, save file
        vars.DataFileName = ['Error_',vars.DataFileName];
        %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
        save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );
        % and as .csv
        csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
        struct2csv(Results, csvName);
        
        disp(['Run crashed. Results were saved as: ', vars.DataFileName]);
        disp(' ** Error!! ***')
        %     ListenChar(0);
        ShowCursor;
        %     sca;
        
        % Output the error message that describes the error:
        psychrethrow(psychlasterror);
        rethrow(ME)
    end
end

if vars.RunSuccessfull  % Successfull run
    % Show end screen and clean up
    Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
    DrawFormattedText(scr.win, vars.InstructionEnd, 'center', 'center', scr.TextColour);
    [~, ~] = Screen('Flip', scr.win);
    WaitSecs(6);
    
    % Save the data
    %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
    save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );
    disp(['Run complete. Results were saved as: ', vars.DataFileName]);
    % and as .csv
    csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
    struct2csv(Results, csvName);                                       %<----- PsiAdaptive: NOT SAVING .csv due to PF objects in Results struct#####
    
end

%% close COM port
if vars.stimFlag 
    TcsAbortStimulation(vars.ser)
    TcsSetBaseLine(vars.ser,32)
    TcsCloseCom(vars.ser)
end
% rmpath(genpath('code'));
% ListenChar(0);          % turn on keypresses -> command window
% sca;
ShowCursor;
% fclose('all');
% Priority(0);
end