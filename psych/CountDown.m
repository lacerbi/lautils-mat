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
function [displayInfo, keyInfo] = CountDown(scrInfo, enabledkeys, message, duration, textColor)

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

% Clear screen
% Screen('TextFont', scrInfo.win, 'Times Roman');
% Screen('TextSize', scrInfo.win, 24);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
Screen('Flip', scrInfo.win);

% Display first text
DrawFormattedText(scrInfo.win, message{1}, 'center', scrInfo.height/6, textColor);
[displayInfo.VBLTimestamp displayInfo.StimulusOnsetTime ...
        displayInfo.FlipTimestamp displayInfo.Missed displayInfo.Beampos] = ...
        Screen('Flip', scrInfo.win, 0, 1, 2);
    
% Start countdown
timeStart = GetSecs();
timeToGo = duration;

eraserect = round([scrInfo.width/4 scrInfo.height*4/9 scrInfo.width*3/4 scrInfo.height*5/9]);

while timeToGo > 0
    Screen('FillRect', scrInfo.win, scrInfo.backgroundColor, eraserect);    
    DrawFormattedText(scrInfo.win, num2str(round(timeToGo)), 'center', 'center', textColor);
    Screen('Flip', scrInfo.win, 0, 1, 2);
    timeToGo = duration - GetSecs() + timeStart;
end
    
% Display second text
DrawFormattedText(scrInfo.win, message{2}, 'center', scrInfo.height*3/4, textColor);
[displayInfo.VBLTimestamp displayInfo.StimulusOnsetTime ...
        displayInfo.FlipTimestamp displayInfo.Missed displayInfo.Beampos] = ...
        Screen('Flip', scrInfo.win, 0, 1, 2);
    
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
    
% Wait for keypress
[keyInfo.downTime, keyInfo.keyCode, keyInfo.deltaSecs] = KbWait();
while 1
    [keyIsDown, keyInfo.upTime] = KbCheck();
    if ~keyIsDown; break; end     
end

% Wait some random time and clear screen
WaitRandomSecs(scrInfo.randomWait);
Screen('Flip', scrInfo.win);

% Re-enable old keys
DisableKeysForKbCheck(olddisabledkeys);

end