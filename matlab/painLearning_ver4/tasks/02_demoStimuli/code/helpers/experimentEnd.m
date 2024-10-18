function experimentEnd(vars, scr, keys, Results)
%
% Project: Implementation of QST thermal thresholds
%
% End of experiment routine. Shows a message to let the user know if the
% run has been aborted or crashed, saves results, and cleans up

%
% Niia Nikolova & Francesca Fardo 
% Last edit: 02/07/2021

if isfield(vars,'Aborted') || isfield(vars,'Error')
    if vars.Aborted
        % Abort screen
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        DrawFormattedText(scr.win, 'Run was stopped by experimenter...', 'center', 'center', scr.TextColour);
        [~, ~] = Screen('Flip', scr.win);
        WaitSecs(0.5);
        ShowCursor;
        %     sca;
        disp('Experiment aborted by user!');
        
        % Save, mark the run
        vars.DataFileName = ['Aborted_', vars.DataFileName];
        %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
        save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );
        disp(['Run was aborted. Results were saved as: ', vars.DataFileName]);
        
        % Display message, close win and return to project directory
        WaitSecs(3);
        sca;
        cd(vars.projdir);
        % and as .csv
        %csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
        %struct2csv(Results, csvName);
        
    elseif vars.Error
        % Abort screen
        Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
        DrawFormattedText(scr.win, 'Run crashed ...', 'center', 'center', scr.TextColour);
        [~, ~] = Screen('Flip', scr.win);
        
        ShowCursor;
        %     sca;
        disp(' ** Error!! ***')
        disp(['Run crashed. Results were saved as: ', vars.DataFileName]);
        
        % Error
        vars.DataFileName = ['Error_',vars.DataFileName];
        save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );
        %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
        % and as .csv
        %csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
        %struct2csv(Results, csvName);
        
        % Output the error message that describes the error:
        psychrethrow(psychlasterror);
        
        % Display message, close win and return to project directory
        WaitSecs(3)
        sca;
        cd(vars.projdir);
    end
end

if vars.RunSuccessfull      % Successfull run
    % Show end screen and clean up
    %Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
    %DrawFormattedText(scr.win, vars.InstructionEnd, 'center', 'center', scr.TextColour);
    %[~, ~] = Screen('Flip', scr.win);
    %     WaitSecs(3);
    ShowCursor;
    %     sca;
    
    % Save the data
    %save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys' );
    save(strcat(vars.OutputFolder, vars.DataFileName), 'Results', 'vars', 'scr', 'keys', '-regexp', ['^(?!', 'vars.ser' , '$).'] );
    disp(['Run complete. Results were saved as: ', vars.DataFileName]);
    % Create dataframe for analysis
    %mkdf(vars,Results) 
%     copy2VMPaux(subID)
    
    % and as .csv
    %csvName = strcat(vars.OutputFolder, vars.DataFileName, '.csv');
    %struct2csv(Results, csvName);
    
end


%% close COM port
if vars.stimFlag 
    TcsAbortStimulation(vars.ser)
    TcsSetBaseLine(vars.ser,32)
    TcsCloseCom(vars.ser)
end
end