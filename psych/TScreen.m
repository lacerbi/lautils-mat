% TSCREEN Initialize or finalize the screen for the time experiments.
%   SCRINFO = INITSCREEN() returns a structure that contains all
%   information about the screen. The fields of the structure are
%   scrInfo.olddebuglevel, scrInfo.oldSupressAllWarnings
%   scrInfo.screens, scrInfo.scrNumber
%   scrInfo.win, scrInfo.rect, scrInfo.center
%   scrInfo.monitorFlipInterval, scrInfo.black, scrInfo.white
%   scrInfo.backgroundColor, scrInfo.foregroundColor
%   scrInfo.randomWait
function scrInfo = TScreen(action, scrInfo, varargin)

if strcmpi(action, 'open')
    
    % The first variable argument is the colormap
    if length(varargin) > 0
        scrcolormap = varargin{1};
    else
        scrcolormap = [0 0 0; 125 125 125; 200 200 200; 255 255 0; 255 0 255];        
    end    
    
    % Fix some bugs in the video driver using data for the 
    % ATI Radeon HD 2400 PRO graphic card plus
    % DELL M782p monitor. See also 'help Beampositionqueries'
    vtotal = 527;
    Screen('Preference', 'ConserveVRAM', 4096);
    Screen('Preference', 'VBLEndlineOverride', vtotal);

    scrInfo.oldVisualDebugLevel=Screen('Preference', 'VisualDebugLevel', 3);
    scrInfo.oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    scrInfo.screens=Screen('Screens');
    scrInfo.scrNumber=max(scrInfo.screens);

    [scrInfo.win,scrInfo.rect] = Screen('OpenWindow',scrInfo.scrNumber);
    scrInfo.center = round([scrInfo.rect(3:4)/2 scrInfo.rect(3:4)/2]);
    scrInfo.width = scrInfo.rect(3) - scrInfo.rect(1);
    scrInfo.height = scrInfo.rect(4) - scrInfo.rect(2);
    scrInfo.unitSize = round(scrInfo.width/1280*12);

    % Check the frame rate of the screen
    scrInfo.monitorFlipInterval=Screen('GetFlipInterval', scrInfo.win, 1, 0.0001);
    
    % The average lag is half the flip interval
    scrInfo.monitorAverageLag = scrInfo.monitorFlipInterval/2;

    % Find the colours; we impose a dark gray background and light gray text
    scrInfo.black = BlackIndex(scrInfo.win);
    scrInfo.white = WhiteIndex(scrInfo.win);    
    
    scrInfo.backgroundColor = scrcolormap(1, :);
    scrInfo.foregroundColor = scrcolormap(2, :);
    scrInfo.textColor = scrcolormap(3, :);
    
    scrInfo.stimulusColor = scrcolormap(4, :);
    scrInfo.spuriousStimulusColor = scrcolormap(5, :);

    Screen('FillRect', scrInfo.win, scrInfo.backgroundColor);
    
    HideCursor;

    % Default Text format
    Screen('TextFont', scrInfo.win, 'Times Roman');
    Screen('TextSize', scrInfo.win, 2*scrInfo.unitSize);
    Screen('TextColor', scrInfo.win, scrInfo.textColor);
    
    % Specify time interval for random wait
    scrInfo.randomWait = [0.5 1];
    
% Close the screen with data SCRINFO
elseif strcmpi(action, 'close')
    Screen('CloseAll');
    ShowCursor;
    Screen('Preference', 'VisualDebugLevel', scrInfo.oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', scrInfo.oldSupressAllWarnings);
    scrInfo = [];
    
% Unknown action
else
    disp(['Unknown option ''', action, ''' for command TScreen. Digit ''help TScreen'' to see the allowed options.']);
    return
end
    
end