% CHOICEROTOR Display some text on screen and a rotating bar.
%   [DISPLAYINFO, KEYINFO] = CHOICEBAR(SCRINFO, MESSAGE [, TEXTCOLOR]) 
%   displays on screen SCRINFO the text MESSAGE and the slider bar.
%   Waits for user response after setting the slider.
%   Returns the structures DISPLAYINFO, KEYINFO which contain all the
%   temporal information regarding the display of text on screen and the
%   user response. All timestamps are expressed in seconds.
%
function [displayInfo, answerInfo] = ChoiceRotor(scrInfo, message, textColor, sliderColor, nticks, slidertitle, titlepos, arrowlabel, ticktext, feedback)

% Disable all keys except arrows, spacebar and escape
allkeys = 1:256;
enabledkeys = KbName({'Backspace', 'Space', 'LeftArrow', 'RightArrow'});
disabledkeys = allkeys(~ismember(allkeys, enabledkeys));
olddisabledkeys = DisableKeysForKbCheck(disabledkeys);

% If not specified, default textColor is foregroundColor
if ~exist('textColor','var'); textColor = []; end
if isempty(textColor); textColor = scrInfo.textColor; end

% If not specified, ticktext is empty
if ~exist('ticktext','var'); ticktext = []; end

% If not specified, there is no feedback
if ~exist('feedback','var'); feedback = []; end

% Clear screen
% Screen('TextFont', scrInfo.win, 'Times Roman');
% Screen('TextSize', scrInfo.win, 24);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
Screen('Flip', scrInfo.win);

% Display introduction message
DrawFormattedText(scrInfo.win, message, 'center', scrInfo.height/6, textColor);

% Slider element
sliderel = floor(scrInfo.width*5/6/nticks);
sliderSize = sliderel*nticks;

% Display slider
sliderpos = [floor((scrInfo.width-sliderSize)/2), scrInfo.height*3/6, sliderSize+floor((scrInfo.width-sliderSize)/2), scrInfo.height*13/24];
midsliderpos = 0.5*[sliderpos(1) + sliderpos(3), sliderpos(2) + sliderpos(4)];

% Randomize initial slider position (from 1 to nticks)
currentselection = floor((nticks+1)/2);
sliderindex = floor(rand()*nticks)+1;
answerInfo.startPosition = sliderindex;
oldsliderindex = sliderindex;
DrawRotor();

% Display screen
[displayInfo.VBLTimestamp displayInfo.StimulusOnsetTime ...
        displayInfo.FlipTimestamp displayInfo.Missed displayInfo.Beampos] = ...
        Screen('Flip', scrInfo.win, 0, 1, 2);

answerInfo.startTime = GetSecs();
    
% Start loop
moving = 0;
waitmoving = 10;
idleTime = 0.02;
while 1
    [keyIsDown, secs, keyCode] = KbCheck();
    if keyIsDown
        % Increase slider
        if keyCode(KbName('RightArrow'))
            if moving < 1 || moving > waitmoving
                oldsliderindex = sliderindex;
                sliderindex = mod(sliderindex, nticks)+1;
                DrawRotor();
            end
            moving = max(1, moving + 1);
        % Move arrow left
        elseif keyCode(KbName('LeftArrow'))
            if moving < -waitmoving || moving > -1
                oldsliderindex = sliderindex;
                sliderindex = mod(sliderindex-2, nticks)+1;
                DrawRotor();
            end
            moving = min(-1, moving - 1);
        elseif keyCode(KbName('Backspace')) || keyCode(KbName('Space'))
            break;
        end
        Screen('Flip', scrInfo.win, 0, 1, 2);
    else
        moving = 0;
    end
    
    % You might want to wait some time in order not to overload the cpu
    WaitSecs('YieldSecs', idleTime); 
end

% The time of the keypress is the end
answerInfo.keyCode = keyCode;
answerInfo.endTime = secs;
answerInfo.endPosition = currentSlider(currentselection);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN FEEDBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(feedback) && ~answerInfo.keyCode(KbName('Backspace'))
    WaitRandomSecs(scrInfo.randomWait);
    
    feedbackColor = [255 255 0];
    
    % Draw feedback arrow
    feedbackindex = mod(feedback - sliderindex, nticks) + 1;
    xpos = sliderpos(1) + round((feedbackindex-0.5)*sliderel);
    Screen('DrawLine', scrInfo.win, feedbackColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
    if ~isempty(arrowlabel)
        DrawFormattedText(scrInfo.win, arrowlabel, xpos - 6, sliderpos(4)+arrowHeight + 4, feedbackColor);         
    end
    
    % Write feedback (only short text)
    cr = strfind(ticktext{answerInfo.endPosition}, '\n');
    if isempty(cr); youranswer = ticktext{answerInfo.endPosition};
    else youranswer = ticktext{answerInfo.endPosition}(1:(cr(1)-1));
    end

    cr = strfind(ticktext{feedback}, '\n');
    if isempty(cr); correctanswer = ticktext{feedback};
    else correctanswer = ticktext{feedback}(1:(cr(1)-1));
    end
    
    if answerInfo.endPosition == feedback
        feedbacktext = ['Correct answer! (' correctanswer ')\n\nPress SPACE to continue.'];
    else
        feedbacktext = ['Your answer: ' youranswer '\n\n'];
        feedbacktext = [feedbacktext 'Correct answer: ' correctanswer '\n\nPress SPACE to continue.'];
    end
    
    textypos = sliderpos(4)+arrowHeight + scrInfo.height/8;
    eraserect = [1, textypos, scrInfo.width, textypos+90];
    Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);

    DrawFormattedText(scrInfo.win, feedbacktext, 'center', sliderpos(4)+arrowHeight + scrInfo.height/8);
    
    % Show the feedback
    Screen('Flip', scrInfo.win);
    
    while 1
        [keyIsDown, secs, keyCode] = KbCheck();
        if keyIsDown
            if keyCode(KbName('Backspace')) || keyCode(KbName('Space'))
                break;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END FEEDBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wait some random time and clear screen
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
WaitRandomSecs(scrInfo.randomWait);
Screen('Flip', scrInfo.win);

% Re-enable old keys
DisableKeysForKbCheck(olddisabledkeys);

return;

    % Draw the arrow
    function drawarrow(x, oldx, label)
        arrowHeight = scrInfo.height/24;
        
        if isempty(label)
            % Delete the previous arrow
            oldxpos = sliderpos(1) + round((oldx-0.5)*sliderel);
            xpos = sliderpos(1) + round((x-0.5)*sliderel);
            Screen('DrawLine', scrInfo.win, scrInfo.backgroundColor, oldxpos, sliderpos(4)+1, oldxpos, sliderpos(4)+arrowHeight, 1);
            Screen('DrawLine', scrInfo.win, scrInfo.foregroundColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
        else
            oldxpos = sliderpos(1) + round((oldx-0.5)*sliderel);
            eraserect = [oldxpos-30, sliderpos(4)+1, oldxpos + 30, sliderpos(4)+100];
            Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);
            xpos = sliderpos(1) + round((x-0.5)*sliderel);
            Screen('DrawLine', scrInfo.win, scrInfo.foregroundColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
            DrawFormattedText(scrInfo.win, label, xpos - 6, sliderpos(4)+arrowHeight + 4);         
        end

        % Draw the associated text
        if ~isempty(ticktext)
            textypos = sliderpos(4)+arrowHeight + scrInfo.height/8;
            eraserect = [1, textypos, scrInfo.width, textypos+90];
            Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);            
            DrawFormattedText(scrInfo.win, ticktext{x}, 'center', sliderpos(4)+arrowHeight + scrInfo.height/8);                        
        end
        
    end    
    
    % Draw the moving title
    function drawtitle()
        textypos = sliderpos(2) - 28;
        eraserect = [1, textypos, scrInfo.width, textypos + 24];
        Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);            
        titleindex = mod(titlepos - sliderindex, nticks) + 1;
        xpos = sliderpos(1) + round((titleindex-1)*sliderel);
        
        DrawFormattedText(scrInfo.win, slidertitle, xpos, textypos);
    end

    % Draw the slider
    function DrawRotor()
        Screen('FrameRect', scrInfo.win, scrInfo.textColor, sliderpos + [-1, -1, 1, 1]);
        elementbase = [sliderpos(1), sliderpos(2), sliderpos(1) + sliderel - 1, sliderpos(4)];
                
        for i = 1:nticks
            Screen('FillRect', scrInfo.win, sliderColor(currentSlider(i), :), elementbase + sliderel*(i-1)*[1 0 1 0]);
        end
        
        % Draw the sliding title
        drawtitle();
        
        % Draw the arrow
        drawarrow(currentselection, currentselection, arrowlabel);
        
        % Draw the associated text
        if ~isempty(ticktext)
            textypos = sliderpos(4) +arrowHeight + scrInfo.height/8;
            eraserect = [1, textypos, scrInfo.width, textypos+90];
            Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);
            
            DrawFormattedText(scrInfo.win, ticktext{currentSlider(currentselection)}, 'center', sliderpos(4)+arrowHeight+ scrInfo.height/8);                        
        end

        
    end

    function jndex = currentSlider(index)
        jndex = mod(sliderindex + index - 2, nticks) + 1;
    end

end

% Fill a rectangle with a specific gradient (from left to right)
function FillGradient(screenwin, rect, colleft, colright)
for x = rect(1):rect(3)
    alpha = 1 - (x - rect(1))/(rect(3) - rect(1));
    gradcolor = alpha*colleft + (1-alpha)*colright;
    Screen('DrawLine', screenwin, gradcolor, x, rect(2), x, rect(4), 1); 
end
end


