% CHOICEBAR Display some text on screen and a slider bar.
%   [DISPLAYINFO, KEYINFO] = CHOICEBAR(SCRINFO, MESSAGE [, TEXTCOLOR]) 
%   displays on screen SCRINFO the text MESSAGE and the slider bar.
%   Waits for user response after setting the slider.
%   Returns the structures DISPLAYINFO, KEYINFO which contain all the
%   temporal information regarding the display of text on screen and the
%   user response. All timestamps are expressed in seconds.
%
function [displayInfo, answerInfo] = ChoiceBar(scrInfo, message, textColor, nticks, slidertitle, arrowlabel, ticktext, feedback)

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

% Display slider
sliderpos = [scrInfo.width/12, scrInfo.height*3/6, scrInfo.width*11/12, scrInfo.height*13/24];
midsliderpos = 0.5*[sliderpos(1) + sliderpos(3), sliderpos(2) + sliderpos(4)];

% Fill the slider
slidergradient = [255 255 0; 50 50 50; 0 255 255];
FillGradient(scrInfo.win, [sliderpos(1)+1, sliderpos(2)+1, midsliderpos(1)-1, sliderpos(4)-1], slidergradient(1, :), slidergradient(2, :));
FillGradient(scrInfo.win, [midsliderpos(1)+1, sliderpos(2)+1, sliderpos(3)-1, sliderpos(4)-1], slidergradient(2, :), slidergradient(3, :));
Screen('FrameRect', scrInfo.win, scrInfo.textColor, sliderpos);

% Display title of the slider
DrawFormattedText(scrInfo.win, slidertitle, 'center', sliderpos(2) - 36);

% Randomize initial arrow position (from 1 to nticks)
arrow = floor(rand()*nticks)+1;
answerInfo.startPosition = arrow;

oldarrow = arrow;
drawarrow(arrow, oldarrow, arrowlabel);

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
        % Move arrow right
        if keyCode(KbName('RightArrow'))
            if moving < 1 || moving > waitmoving
                oldarrow = arrow;
                arrow = min(arrow+1, nticks);
                drawarrow(arrow, oldarrow, arrowlabel);
            end
            moving = max(1, moving + 1);
        % Move arrow left
        elseif keyCode(KbName('LeftArrow'))
            if moving < -waitmoving || moving > -1
                oldarrow = arrow;
                arrow = max(arrow-1, 1);
                drawarrow(arrow, oldarrow, arrowlabel);
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
answerInfo.endPosition = arrow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN FEEDBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(feedback) && ~answerInfo.keyCode(KbName('Backspace'))
    WaitRandomSecs(scrInfo.randomWait);
    
    feedbackColor = [255 255 0];
    
    % Draw feedback arrow
    xpos = ((feedback-1) / (nticks-1)) * (sliderpos(3) - sliderpos(1)) + sliderpos(1);
    Screen('DrawLine', scrInfo.win, feedbackColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
    if ~isempty(arrowlabel)
        DrawFormattedText(scrInfo.win, arrowlabel, xpos - 4, sliderpos(4)+arrowHeight + 4, feedbackColor);         
    end
    
    % Write feedback
    if arrow == feedback
        feedbacktext = ['Correct answer! (' ticktext{feedback} ')\n\nPress SPACE to continue.'];
    else
        feedbacktext = ['Your answer: ' ticktext{arrow} '\n\n'];
        feedbacktext = [feedbacktext 'Correct answer: ' ticktext{feedback} '\n\nPress SPACE to continue.'];
    end
    
    textypos = sliderpos(4)+arrowHeight + scrInfo.height/8;
    eraserect = [1, textypos, scrInfo.width, textypos+30];
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
            oldxpos = ((oldx-1) / (nticks-1)) * (sliderpos(3) - sliderpos(1)) + sliderpos(1);
            xpos = ((x-1) / (nticks-1)) * (sliderpos(3) - sliderpos(1)) + sliderpos(1);
            Screen('DrawLine', scrInfo.win, scrInfo.backgroundColor, oldxpos, sliderpos(4)+1, oldxpos, sliderpos(4)+arrowHeight, 1);
            Screen('DrawLine', scrInfo.win, scrInfo.foregroundColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
        else
            oldxpos = ((oldx-1) / (nticks-1)) * (sliderpos(3) - sliderpos(1)) + sliderpos(1);
            eraserect = [oldxpos-30, sliderpos(4)+1, oldxpos + 30, sliderpos(4)+100];
            Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);
            xpos = ((x-1) / (nticks-1)) * (sliderpos(3) - sliderpos(1)) + sliderpos(1);
            Screen('DrawLine', scrInfo.win, scrInfo.foregroundColor, xpos, sliderpos(4)+1, xpos, sliderpos(4)+arrowHeight, 1);
            DrawFormattedText(scrInfo.win, label, xpos - 4, sliderpos(4)+arrowHeight + 4);         
        end

        % Draw the associated text
        if ~isempty(ticktext)
            textypos = sliderpos(4)+arrowHeight + scrInfo.height/8;
            eraserect = [1, textypos, scrInfo.width, textypos+30];
            Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);
            
            DrawFormattedText(scrInfo.win, ticktext{x}, 'center', sliderpos(4)+arrowHeight + scrInfo.height/8);                        
        end
        
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