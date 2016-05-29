function varargout = MouseTapper(cmd, tapper, varargin)
% MouseTapper - Use mouse as a tapping device.
% 
%
% Commands and their syntax:
% --------------------------
%
% tapper = MouseTapper('Open', polhemus [, activeButton=1]);
% - Open connection with mouse tapper. Initialize it, return a struct 
%   'tapper'. You'll have to pass 'tapper' to all following functions to 
%   access the tapping device. If you provide a non-empty tapper handle, 
%   Open will preserve the existing data.
%   If you want to open a new connection from scratch, pass 'tapper' as
%   an empty structure. The optional argument 'activeButton' specifies the
%   active tapping mouse button.
%
%
% tapper = MouseTapper('Close', tapper);
% - Close connection to mouse tapping device with handle 'tapper'.
%
%
% tapper = MouseTapper('Start', tapper);
% - Start reading session from the mouse tapper 'tapper'. All the buffers 
%   are reset to zero. 
%
%
% tapper = MouseTapper('Stop', tapper);
% - Stop reading data from open connection 'tapper'. Notice that the
%   buffers are not deleted, so that they can be saved or analyzed.
%
%
% [tapper, evt] = MouseTapper('GetEvent', tapper [, waitForEvent=0]);
% - Check mouse for an event (button press or button release), place it in
% struct 'evt' and in tapper buffer.
% If no event is happening and the optional 'waitForEvent' is set to
% 1, then the function will wait until at least one valid event becomes
% available and return that event. Otherwise it will return an empty struct,
% ie., evt = [] to signal that no events are happening.
%
% The following subfields are available in 'evt' if 'evt' is non-empty:
%
% evt.state = a value which encodes the new status of the tapper.
% 1 == tap, 2 == untap.
%
% evt.time  = Psychtoolbox GetSecs() timestamp of the time when the event
% was received from the device. The accuracy depends on the properties of
% the mouse and system load.
%
% evt.trouble = If zero, then evt is probably valid and good. If non-zero,
% then the timestamp is likely screwed and useless, as are probably all
% following timestamps and events. (Not supported yet.)
%
%
% updatedTimetable = PolhemusTapper('UpdateTimeTable', polhemus [, timetable]);
% - If 'timetable' is empty, use buffered polhemus data to generate a phase
% space dataset AND calibrate the tapping device. If 'timetable' is not
% empty, use 'polhemus' buffered data to update the phase space data.
%
%
% tapPrediction = PolhemusTapper('GetTapPrediction', polhemus);
% - Return expected time of next tap (in seconds from now) given the 
% current state and position, as recorded from the last GetEvent.
% If no prediction is achieved, return NaN.

if nargin < 1
    help MouseTapper;
    return;
end

if nargin < 2
    tapper = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GETEVENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(cmd, 'GetEvent')    
    checkInputArgument('GetEvent');
        
    if length(varargin) < 1
        waitEvent = [];
    else
        waitEvent = varargin{1};
    end
    
    if isempty(waitEvent); waitEvent = 0; end
    
    % Start with empty event 'evt':
    evt = [];
    
    % Reset timestamp trouble flag:
    tTrouble = 0;
    
    % Repeat loop till you find an event (only if waitEvent is TRUE)    
    while 1        
        % Read mouse state
        % (ignore x, y positions, consider only left mouse button)
        timeNow = GetSecs();
        [x,y,buttons] = GetMouse();
        timeNow = 0.5*(GetSecs() + timeNow);

        % If button WAS down and now is up, UNTAP
        if tapper.tappingStatus == 1 && buttons(tapper.activeButton) == 0
            evt.time = timeNow;
            evt.state = 2; % UNTAP
            evt.trouble = tTrouble;
            tapper.tappingStatus = 2;
            display('Untap');

        % If button WAS up and now is down, TAP    
        elseif tapper.tappingStatus == 2 && buttons(tapper.activeButton) == 1
            evt.time = timeNow;
            evt.state = 1; % TAP
            evt.trouble = tTrouble;
            tapper.tappingStatus = 1;
            display('Tap');        
        end

        % If an event happened, bufferize it
        if ~isempty(evt); 
            tapper.bufEvents{length(tapper.bufEvents)+1} = evt;
            break;
        elseif waitEvent == 0; 
            break;
        end
        
        % You might want to wait some time in order not to overload the cpu
        WaitSecs('YieldSecs', 0.00025); 
    end
   
    % Return evt if any:
    varargout{1} = tapper;
    varargout{2} = evt;
    
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open connection to MouseTapper:
if strcmpi(cmd, 'Open')
    if length(varargin) < 1; activeButton = [];
    else activeButton = varargin{1}; end
        
    if isempty(activeButton); activeButton = 1; end
    
    % Check if a non-empty tapper handle is provided
    if isempty(tapper); 
        notapper = 1; 
        clear tapper;
    else
        checkInputArgument('Open');
        notapper = 0;
    end
        
    tapper.activeButton = activeButton;
    tapper.recording = 0;
    tapper.tappingStatus = 2; % Starting position is UNTAP

    % Initialize tapper temporary collected data
    if notapper
        tapper.bufEvents = []; 
        % Perform a test?        
        fprintf('MouseTapper: Connection open!\n');
    end
    
    % Preheat GetSecs:
    GetSecs;
    
    varargout{1} = tapper;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close connection to tapping device:
if strcmpi(cmd, 'Close')
    checkInputArgument('Close');
    
    % It does nothing...
    
    varargout{1} = polhemus;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start tapping session:
if strcmpi(cmd, 'Start')
    checkInputArgument('Start');
   
   tapper.tappingStatus = 2; % Starting position is UP
   tapper.bufEvents = [];
   tapper.recording = 1;
                     
   varargout{1} = polhemus;     
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stop tapping session:
if strcmpi(cmd, 'Stop')
    checkInputArgument('Stop');
    
    tapper.recording = 0;
    varargout{1} = tapper;
        
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESTART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restart tapping session (flush data, set tappingStatus to two):
if strcmpi(cmd, 'Restart')
    checkInputArgument('Restart');
   
   tapper.tappingStatus = 2; % Starting position is UP
   tapper.bufEvents = [];
   tapper.recording = 1;
                     
   varargout{1} = polhemus;     
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVEBUFFER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the buffer (a short version of the tapper struct)
if strcmpi(cmd, 'SaveBuffer')
    checkInputArgument('SaveBuffer');
    ttemp.bufEvents = tapper.bufEvents;
    
    varargout{1} = ttemp;    
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATETIMETABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create or update the time table starting from the given buffered data
if strcmpi(cmd, 'UpdateTimeTable')
    checkInputArgument('UpdateTimeTable');
    
    % The optional argument should be a Polhemus timetable struct
    if size(varargin, 1) > 0;
        timetable = varargin{1};
    else
        timetable = [];
    end
    
    evtList = polhemus.bufEvents;
    
    % Update timetable data (proper events timing)
    [tapList, timetable] = NewTapList(evtList, dataList, timetable);
            
    % Update phase space data for all proper events
    for i = 1:size(tapList, 1)
        evtIndexes = tapList(i, :);
        evtTimes = [evtList{evtIndexes(1)}.time, evtList{evtIndexes(2)}.time, evtList{evtIndexes(3)}.time, evtList{evtIndexes(4)}.time];            
        tapTime = evtTimes(4);
        % evtTimesToTap = tapTime - [tapTime, evtTimes(1), evtTimes(2), evtTimes(3)];

        % Update phase space data
        for j = (evtList{evtIndexes(1)}.index+2):evtList{evtIndexes(4)}.index
            ypos = polhemus.bufData(j, 2);
            tpos = tapTime - polhemus.bufTimes(j);
            % Avoid points exactly at the tapping time (they are
            % not informative, and the calculation of the derivative
            % might be erroneous too!)
            if tpos == 0; continue; end
            % Calculate velocity (the function also check that timing
            % is reasonably spaced)
            vpos = getVelocity(j);
            phasespace = [ypos vpos tpos];
            if any(isnan(phasespace)); continue; end
            timetable.phasespaceData = [timetable.phasespaceData; phasespace];
        end
                
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN CONSTRUCTION OF PHASE SPACE LOOKUP TABLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tablesize = [100 100];
    for i1 = 1:tablesize(1);
        for i2 = 1:tablesize(2);
            timetable.phasespaceLookupTable{i1}{i2} = [];
        end
    end
    % Not all phase space history is checked -- up to last 10k entries
    lookuptablemaxsize = 10000;
    lookupmax = size(timetable.phasespaceData, 1);
    lookupmin = max(lookupmax - lookuptablemaxsize, 1);
    
    % Build timetable only if there are some data
    if lookupmax > 0
        timetable.phasespacemins = min(timetable.phasespaceData(lookupmin:lookupmax, 1:2));
        timetable.phasespacesteps = (max(timetable.phasespaceData(lookupmin:lookupmax, 1:2)) - timetable.phasespacemins)./(tablesize - 1);            

        % Fill table
        for j = lookupmin:lookupmax
           xindex = floor((timetable.phasespaceData(j, 1) - timetable.phasespacemins(1))/timetable.phasespacesteps(1)) + 1;
           vindex = floor((timetable.phasespaceData(j, 2) - timetable.phasespacemins(2))/timetable.phasespacesteps(2)) + 1;
           timetable.phasespaceLookupTable{xindex}{vindex} = [timetable.phasespaceLookupTable{xindex}{vindex}; timetable.phasespaceData(j, :)];                
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END CONSTRUCTION OF PHASE SPACE LOOKUP TABLE UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout{1} = timetable;
    
    return;
end

% Compute the prediction error over the current buffered data
if strcmpi(cmd, 'ComputePredictionError')
    if isempty(polhemus)
        error('PolhemusTapper: ComputePredictionError: No "handle" for device provided!');
    end
    
    if length(polhemus) ~= 1
        error('PolhemusTapper: ComputePredictionError: Passed argument is not a valid Polhemus "handle"!');
    end
        
    % The optional argument should be a Polhemus timetable struct
    if size(varargin, 1) < 1;
        error('PolhemusTapper: ComputePredictionError: A valid Polhemus timetable must be provided!');
    else
        timetable = varargin{1};
    end
    
    evtList = polhemus.bufEvents;
    dataList = polhemus.bufData;
    timesList = polhemus.bufTimes;
    
    % Update timetable data (proper events timing)
    [tapList, timetable] = NewTapList(evtList, dataList, timetable);
        
    predictionErrorList = [];
                
    % Compute prediction error
    for i = 1:size(tapList, 1)
        evtIndexes = tapList(i, :);
        evtTimes = [evtList{evtIndexes(1)}.time, evtList{evtIndexes(2)}.time, evtList{evtIndexes(3)}.time, evtList{evtIndexes(4)}.time];            
        tapTime = evtTimes(4);                        

        timetable.evtTimesToTap(timetable.indexEvt(1), 1) = 0;
        timetable.evtTimesToTap(timetable.indexEvt(2), 2) = tapTime - evtTimes(1);
        timetable.evtTimesToTap(timetable.indexEvt(3), 3) = tapTime - evtTimes(2);
        timetable.evtTimesToTap(timetable.indexEvt(4), 4) = tapTime - evtTimes(3);

        % Increase index in a discrete ring
        timetable.indexEvt = timetable.indexEvt + 1;
        if timetable.indexEvt > size(timetable.evtTimesToTap, 1); timetable.indexEvt = ones(1, 4); end

        % Update tapping errors            
        currentTappingState = 2;
        for j = evtList{evtIndexes(1)}.index:evtList{evtIndexes(4)}.index
            if j >= evtList{evtIndexes(3)}.index; currentTappingState = 4;
            elseif j >= evtList{evtIndexes(2)}.index; currentTappingState = 3;
            end

            ypos = polhemus.bufData(j, 2);
            vpos = getVelocity(j);
            elapsedTimeFromLastEvent = polhemus.bufTimes(j) - evtTimes(currentTappingState-1);
            tapActual = tapTime - polhemus.bufTimes(j);

            tapPrediction = getTapPrediction(ypos, vpos, elapsedTimeFromLastEvent, currentTappingState);
            predictionError = tapPrediction - tapActual;                
            if ~isnan(predictionError)
                predictionErrorList = [predictionErrorList; tapActual predictionError];
            end

        end
    end
    
    % Calculate binned response (10 ms bins)
    nBins = 100;
    binSize = 0.010;
    for k = 1:nBins; peBins{k} = []; end
    for k = 1:size(predictionErrorList, 1);
        binIndex = min(floor(predictionErrorList(k, 1)/binSize) + 1, nBins);
        peBins{binIndex} = [peBins{binIndex} predictionErrorList(k, 2)];
    end    
        
    % Return prediction error list:
    varargout{1} = predictionErrorList;
    varargout{2} = peBins;
    
    return;
end

% Get the prediction from the current Polhemus status
if strcmpi(cmd, 'GetTapPrediction')
    if isempty(polhemus)
        error('PolhemusTapper: GetTapPrediction: No "handle" for device provided!');
    end
    
    if length(polhemus) ~= 1
        error('PolhemusTapper: GetTapPrediction: Passed argument is not a valid Polhemus "handle"!');
    end
    
    % The optional argument should be a Polhemus timetable struct
    if size(varargin, 1) < 1;
        error('PolhemusTapper: ComputePredictionError: A valid Polhemus timetable must be provided!');
    else
        timetable = varargin{1};
    end
    
    tapPrediction = NaN;
    if polhemus.nextBuf > 1 && polhemus.tappingStatus > 1
        ypos = polhemus.bufData(polhemus.nextBuf-1, 2);
        vpos = getVelocity(polhemus.nextBuf-1);
        elapsedTimeFromLastEvent = GetSecs() - polhemus.bufEvents{length(polhemus.bufEvents)}.time;
        tapPrediction = getTapPrediction(ypos, vpos, elapsedTimeFromLastEvent, polhemus.tappingStatus);        
    end
        
    varargout{1} = tapPrediction;
    return;
end

% Invalid command!
error('PolhemusTapper: Invalid or unknown command specified!');

    % GETTAPPREDICTION Calculate the predicted time of the next tap (from
    % now, in seconds) given an y position, a velocity, the time from the 
    % last event and the current tapping status.
    %
    function taptime = getTapPrediction(ypos, vpos, timeFromLastEvent, tappingStatus)
        % yindex = round((ypos - polhemus.timetable.positionBase)/polhemus.timetable.positionStep);
        % yindex = max(1, min(yindex, size(polhemus.timetable.guessedTappingTimes, 1)));
        % taptime = polhemus.timetable.guessedTappingTimes(yindex, tappingStatus);
        % avgevtTimesToTap = nanmean(polhemus.timetable.evtTimesToTap(:, tappingStatus)) - timeFromLastEvent;
        % taptime = polhemus.timetable.positionWeight(tappingStatus)*taptime + (1-polhemus.timetable.positionWeight(tappingStatus))*avgevtTimesToTap;
                    
        % BEGIN PHASE SPACE PREDICTION
        if isempty(timetable.phasespaceLookupTable)
            taptime = NaN;
            return;
        end
        
        % Calculate y index in phase space lookup table
        lookupRadius = 1;
        yindex = floor((ypos - timetable.phasespacemins(1))/timetable.phasespacesteps(1)) + 1;
        vindex = floor((vpos - timetable.phasespacemins(2))/timetable.phasespacesteps(2)) + 1;
        minbounds = [max(yindex - lookupRadius, 1), max(vindex - lookupRadius, 1)];
        maxbounds = [min(yindex + lookupRadius, size(timetable.phasespaceLookupTable, 2)), min(vindex + lookupRadius, size(timetable.phasespaceLookupTable{1}, 2))];
        
        taptime = 0; % Estimated tapping time
        nf = 0; % Normalization factor
        nfdist = timetable.phasespacesteps.^2; % Normalization factor for distances (squared)        
        for yloop = minbounds(1):maxbounds(1)
            for vloop = minbounds(2):maxbounds(2)
                phasespacecell = timetable.phasespaceLookupTable{yloop}{vloop};
                for kloop = 1:size(phasespacecell, 1);
                    normdistance = ([ypos - phasespacecell(kloop, 1), vpos - phasespacecell(kloop, 2)].^2)./nfdist;
                    weight = exp(-sum(normdistance));
                    taptime = taptime + weight*phasespacecell(kloop, 3);
                    nf = nf + weight;
                end
            end
        end
        taptime = taptime / nf;
        
        %if isnan(taptime)
        %    yindex
        %    vindex
        %    minbounds
        %    maxbounds
        %   nf 
        %end
        
        % END PHASE SPACE PREDICTION
    end

% End of main function:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MACROS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CHECKINPUTARGUMENT A macro for checking the validity of the tapper handle.
    function checkInputArgument(cmdString)
        if isempty(tapper)
            error(['MouseTapper: ' cmdString, ': No "handle" for mouse tapper provided!']);
        end

        if length(tapper) ~= 1
            error(['MouseTapper: ' cmdString, ': Passed argument is not a valid tapper "handle"!']);
        end
    end

end

% NEWTIMETABLE Timetable struct constructor. Provide a downposition and a
% upposition (6-vectors) for the tapping device.
function newtt = NewTimeTable(downPosition, upPosition)
   timetableSize = [20 50];
   newtt.downPosition = downPosition;
   newtt.upPosition = upPosition;
   newtt.positionBase = downPosition(2);
   newtt.positionStep = (upPosition(2) + 0.5 - downPosition(2))/(timetableSize(2)-1);
   %newtt.positionTimesToTap = NaN*ones(timetableSize(1), timetableSize(2), 2);
   newtt.evtTimesToTap = NaN*ones(timetableSize(1), 4);
   newtt.guessedTappingTimes = zeros(timetableSize(2), 4);
   newtt.positionWeight = [0.9 0.9 0.9 0.9];
   newtt.indexEvt = ones(1, 4);
   newtt.indexPosition = ones(timetableSize(2), 2);
   newtt.phasespaceData = [];
   newtt.phasespaceLookupTable = [];
   newtt.predictionError = [];
end

% NEWTAPLIST Build a list of valid taps from an event list and data list 
% (and optionally a timetable)
function [tapList, timetable] = NewTapList(evtList, dataList, timetable)

    tapList = [];
    
    % Check existence of timetable, otherwise generate one
    if ~exist('timetable','var'); timetable = []; end    
    if isempty(timetable)
        % Determine average tapping up and down positions
        tapDown = []; 
        tapUp = [];
        for ievt = 1:size(evtList, 2)
            if evtList{ievt}.state == 1; 
                tapDown = [tapDown; dataList(evtList{ievt}.index, :)];
            elseif evtList{ievt}.state == 4; 
                tapUp = [tapUp; dataList(evtList{ievt}.index, :)];
            end
        end
        downPosition = mean(tapDown);
        upPosition = mean(tapUp);
        timetable = NewTimeTable(downPosition, upPosition);   
        clear tapDown tapUp;
    end    

    % Update timetable
    for i = 1:(length(evtList) - 3)
        % Identify a full tapping event
        if evtList{i}.state == 2 && evtList{i+1}.state == 3 && evtList{i+2}.state == 4 && evtList{i+3}.state == 1
            evtIndexes = [i, i+1, i+2, i+3];
            evtTimes = [evtList{evtIndexes(1)}.time, evtList{evtIndexes(2)}.time, evtList{evtIndexes(3)}.time, evtList{evtIndexes(4)}.time];            
            tapTime = evtTimes(4);
                                    
            evtTimesToTap = tapTime - [tapTime, evtTimes(1), evtTimes(2), evtTimes(3)];
            
            % Check for anomaly (if there are at least four datapoints)
            if ~isnan(timetable.evtTimesToTap(4, 1))
                meanTimesToTap = nanmean(timetable.evtTimesToTap, 1);
                stdTimesToTap = max(nanstd(timetable.evtTimesToTap, 1), 0.020);
                distanceTimesToTap = abs(evtTimesToTap - meanTimesToTap);
                % If any measurement is more than FIVE std away, do not use
                % tap -- it is anomalous
                if any(distanceTimesToTap - 5*stdTimesToTap > 0)
                    continue;
                end
            end
            
            tapList = [tapList; evtIndexes];
            
            timetable.evtTimesToTap(timetable.indexEvt(1), 1) = 0;
            timetable.evtTimesToTap(timetable.indexEvt(2), 2) = tapTime - evtTimes(1);
            timetable.evtTimesToTap(timetable.indexEvt(3), 3) = tapTime - evtTimes(2);
            timetable.evtTimesToTap(timetable.indexEvt(4), 4) = tapTime - evtTimes(3);

            % Increase index in a discrete ring
            timetable.indexEvt = timetable.indexEvt + 1;
            if timetable.indexEvt > size(timetable.evtTimesToTap, 1); timetable.indexEvt = ones(1, 4); end

        end
    end
    
    varargout{1} = tapList;
    varargout{2} = timetable;    
end

