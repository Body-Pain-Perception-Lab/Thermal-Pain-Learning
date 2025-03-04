function [vars] = getResponse(keys, scr, vars)
%function [vars] = getResponse(keys, scr, vars)
%
% Get the participants response - either keyboard or mouse
%
% Project: CWT task, for fMRI.
%
% Input:
%   keys (struct)
%   scr (struct)
%   vars (struct)
%
%
% Output:
%   vars (struct)
%
%
% Niia Nikolova
% Last edit: 07/07/2020

feedbackString = 'O';
% loop until valid key is pressed or RespT is reached
while ((GetSecs - vars.StartRT) <= vars.RespT)
    
    switch vars.InputDevice
        
        case 1 % Keyboard response
            
            % KbCheck for response
            if keys.KeyCode(keys.Left)==1         % Angry
                % update results
                vars.Resp = 0;
                vars.ValidTrial(1) = 1;
                
            elseif keys.KeyCode(keys.Right)==1    % Happy
                % update results
                vars.Resp = 1;
                vars.ValidTrial(1) = 1;
                
            elseif keys.KeyCode(keys.Escape)==1
                vars.abortFlag = 1;
                
            else
                % ? DrawText: Please press a valid key...
            end
            
            [~, vars.EndRT, keys.KeyCode] = KbCheck;
            WaitSecs(0.001);
            
            
        case 2 % Mouse
            
            [~,~,buttons] = GetMouse;
            while (~any(buttons)) && ((GetSecs - vars.StartRT) <= vars.RespT) % wait for press & response time
                [~,~,buttons] = GetMouse; % L [1 0 0], R [0 0 1]
            end
            
            if buttons == [1 0 0] % Left - angry
                % update results
                vars.Resp = 0;
                vars.ValidTrial(1) = 1;
                
            elseif buttons == [0 0 1] % Right - happy
                % update results
                vars.Resp = 1;
                vars.ValidTrial(1) = 1;
                
            else
                vars.Resp = NaN;
                vars.ValidTrial(1) = 0;
                
            end
            vars.EndRT = GetSecs;
    end
    
    %% Brief feedback
    if vars.Resp == 1% happy
        emotString = 'Happy';
        feedbackXPos = ((scr.winRect(3)/2)+150);
    elseif vars.Resp == 0
        emotString = 'Angry';
        feedbackXPos = ((scr.winRect(3)/2)-250);
    else
        emotString = '';
        feedbackXPos = ((scr.winRect(3)/2));
    end
    
    % fixed timing - wait for response interval to pass
    if vars.fixedTiming
        if ~isnan(vars.Resp) && (vars.ValidTrial(1))    % valid trial
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, [vars.InstructionQ], 'center', 'center', scr.TextColour);
            DrawFormattedText(scr.win, feedbackString, feedbackXPos, ((scr.winRect(4)/2)+150), scr.AccentColour);
            [~, ~] = Screen('Flip', scr.win);
            
            outputString = ['Response recorded: ', emotString];
        else
            outputString = 'No response recorded';
        end
        
    else    % Variable timing
        if ~isnan(vars.Resp) && (vars.ValidTrial(1))
            Screen('FillRect', scr.win, scr.BackgroundGray, scr.winRect);
            DrawFormattedText(scr.win, [vars.InstructionQ], 'center', 'center', scr.TextColour);
            DrawFormattedText(scr.win, feedbackString, feedbackXPos, ((scr.winRect(4)/2)+150), scr.AccentColour);
            [~, ~] = Screen('Flip', scr.win);
            WaitSecs(0.3);
            
            outputString = ['Response recorded: ', emotString];
        else
            outputString = 'No response recorded';
        end
        
        WaitSecs(0.2);
        break;
    end
    
end

disp(outputString);

end
