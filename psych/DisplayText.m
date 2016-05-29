function [responseKey,RT,displayInfo,keyInfo]=DisplayText(scrInfo, enabledkeys, msg, textColor, timeLimit)
% DISPLAYTEXT Display some text on screen and wait for a keypress.
%   [RESPONSEKEY,RT] = DISPLAYTEXT(SCRINFO,ENABLEDKEYS,MSG) displays 
%   on screen SCRINFO the text MSG and waits for user response
%   through ENABLEDKEYS. ENABLEDKEYS is either an array of allowed keys, or a struct of key
%   names. Returns the first pressed key RESPONSEKEY and the reaction time
%   RT in seconds (time from screen display to key press).
%
%   [RESPONSEKEY,RT,DISPLAYINFO,KEYINFO] = DISPLAYTEXT(...) also returns
%   the structures DISPLAYINFO, KEYINFO which contain detailed information 
%   regarding the display of text on screen and the user response. All 
%   timestamps are expressed in seconds.
%
%   [...]=DISPLAYTEXT(SCRINFO,ENABLEDKEYS,MSG,TEXTCOLOR) 
%   draws text in TEXTCOLOR. Default color is SCRINFO.FOREGROUNDCOLOR.
%
%   [...]=DISPLAYTEXT(SCRINFO,ENABLEDKEYS,MSG,TEXTCOLOR,TIMELIMIT) sets a
%   time limit of TIMELIMIT seconds for the response.
%
%   Note: RESPONSEKEY returns the code of one single key. If simultaneous
%   multiple keys are needed, data are in KEYINFO.KEYCODE.

% Do nothing if message is not provided
if isempty(msg)
   displayInfo = [];
   keyInfo.keyCode = [];
   keyInfo.responseTime = 0;
   responseKey = NaN;
   RT = NaN;
   return;
end

% If not specified, default textColor is foregroundColor
if nargin<4 || isempty(textColor); textColor = scrInfo.textColor; end    

% If not specified, default timeLimit is infinity
if nargin<5 || isempty(timeLimit); timeLimit = Inf; end    

% Get list of enabled keys
if isempty(enabledkeys) % All keys are enabled
    enabledkeyslist = 1:256;
elseif isnumeric(enabledkeys)
    enabledkeyslist = enabledkeys(:)';
else
    enabledkeyslist = [];
    for i = 1:length(enabledkeys)
        enabledkeyslist = [enabledkeyslist, KbName(enabledkeys{i})];
    end    
end

% Disable not enabled keys
allkeys = 1:256;
disabledkeys = allkeys(~ismember(allkeys, enabledkeyslist));
olddisabledkeys = DisableKeysForKbCheck(disabledkeys);

enabledkeyslist

% Clear screen
% Screen('TextFont', scrInfo.win, 'Times Roman');
% Screen('TextSize', scrInfo.win, 24);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
Screen('Flip', scrInfo.win);

% Display text and prepare blank screen
DrawFormattedText(scrInfo.win, msg, 'center', 'center', textColor);
[displayInfo.VBLTimestamp displayInfo.StimulusOnsetTime ...
        displayInfo.FlipTimestamp displayInfo.Missed displayInfo.Beampos] = ...
        Screen('Flip', scrInfo.win);
Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);

% Wait for keypress

startTime = GetSecs();

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

% Response key (only one key is reported)
responseKey = find(keyInfo.keyCode, 1);

% Reaction time
RT = keyInfo.downTime - displayInfo.StimulusOnsetTime;

% Wait some random time and clear screen
WaitRandomSecs(scrInfo.randomWait);
Screen('Flip', scrInfo.win);

% Re-enable old keys
DisableKeysForKbCheck(olddisabledkeys);

end