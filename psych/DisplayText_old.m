% DISPLAYTEXT Display some text on screen and wait for a keypress.
%   [DISPLAYINFO, KEYINFO] = DISPLAYTEXT(SCRINFO, ENABLEDKEYS, MESSAGE [, TEXTCOLOR]) 
%   displays on screen SCRINFO the text MESSAGE and waits for user response
%   through ENABLEDKEYS. ENABLEDKEYS is an array of allowed keys, and it
%   covers keyboard and other devices. The keymap is:
%   1:256 keyboard. (Other devices to be implemented.)
%   Returns the structures DISPLAYINFO, KEYINFO which contain all the
%   temporal information regarding the display of text on screen and the
%   user response. All timestamps are expressed in seconds.
%
%   [DISPLAYINFO, KEYINFO] = DISPLAYTEXT(SCRINFO, ENABLEDKEYS, MESSAGE, 
%   TEXTCOLOR) draws text in color TEXTCOLOR. By default,
%   SCRINFO.FOREGROUNDCOLOR is used.
%
function [displayInfo, keyInfo] = DisplayText(scrInfo, enabledkeys, message, textColor, timeLimit)

% Do nothing if message is not provided
if isempty(message)
   displayInfo = [];
   keyInfo.keyCode = [];
   keyInfo.responseTime = 0;
   return;
end

% Disable not enabled keys
if isempty(enabledkeys)
    olddisabledkeys = DisableKeysForKbCheck([]);
else
    allkeys = 1:256;
    disabledkeys = allkeys(~ismember(allkeys, enabledkeys));
    olddisabledkeys = DisableKeysForKbCheck(disabledkeys);
end

% If not specified, default textColor is foregroundColor
if ~exist('textColor','var'); textColor = []; end
if isempty(textColor); textColor = scrInfo.textColor; end    

% If not specified, default timeLimit is infinity
if ~exist('timeLimit','var'); timeLimit = []; end
if isempty(timeLimit); timeLimit = Inf; end    

% Clear screen
% Screen('TextFont', scrInfo.win, 'Times Roman');
% Screen('TextSize', scrInfo.win, 24);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
Screen('Flip', scrInfo.win);

% Display text and prepare blank screen
DrawFormattedText(scrInfo.win, message, 'center', 'center', textColor);
[displayInfo.VBLTimestamp displayInfo.StimulusOnsetTime ...
        displayInfo.FlipTimestamp displayInfo.Missed displayInfo.Beampos] = ...
        Screen('Flip', scrInfo.win);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);

% Wait for keypress

startTime = GetSecs();

% [keyInfo.downTime, keyInfo.keyCode, keyInfo.deltaSecs] = KbWait();
while 1
    elapsedTime = GetSecs() - startTime;
    [keyIsDown, keyInfo.downTime, keyInfo.keyCode] = KbCheck();
    if keyIsDown || elapsedTime > timeLimit; break; end
end

if elapsedTime <= timeLimit
    keyInfo.deltaSecs = elapsedTime;
    % Wait till key is released
    while 1
        [keyIsDown, keyInfo.upTime] = KbCheck();
        if ~keyIsDown; break; end
    end
    keyInfo.timeout = 0;
else
    keyInfo.deltaSecs = Inf;    
    keyInfo.timeout = 1;
end

% Wait some random time and clear screen
WaitRandomSecs(scrInfo.randomWait);
Screen('Flip', scrInfo.win);

% Re-enable old keys
DisableKeysForKbCheck(olddisabledkeys);

end