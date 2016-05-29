function varargout = PolhemusTapper(cmd, polhemus, varargin)
% PolhemusTapper - Access Polhemus Liberty motion tracker as a tapping device.
% 
%
% Commands and their syntax:
% --------------------------
%
% polhemus = PolhemusTapper('Open', polhemus, activeSensorNumber);
% - Open connection with Polhemus motion tracker, considering the
% 'activeSensorNumber'-th sensor as the tapping device. Initialize it, 
% return a struct 'polhemus'. You'll have to pass 'polhemus' to all 
% following functions to access the tapping device. If you provide
% a non-empty polhemus handle, Open will preserve the existing data.
% If you want to open a new connection from scratch, pass 'polhemus' as
% an empty structure.
%
%
% polhemus = PolhemusTapper('Close', polhemus);
% - Close connection to Polhemus device with handle 'polhemus'.
%
%
% [polhemus, notInPosition] = 
% PolhemusTapper('Start', polhemus [, WaitTillInPosition=0, bufSize=5000]);
% - Start reading data from open connection 'polhemus'. If the flag
% 'WaitTillInPosition' is 1, the device waits till the sensor is in the
% 'resting' position on the tapping pad. 'bufSize' is the starting size of 
% the buffer matrix. All the buffers are reset to zero. 
% Returns the error flag 'notInPosition' which is 1 if the sensor is not in
% the starting position.
%
%
% polhemus = PolhemusTapper('Stop', polhemus);
% - Stop reading data from open connection 'polhemus'. Notice that the
% buffers are not deleted, so that they can be saved or analyzed.
%
%
% [polhemus, timetable] = PolhemusTapper('Calibrate', polhemus);
% - Active calibration of the tapping device. Return the 'polhemus' struct
% and a initialized 'timetable' struct.
%
%
% [polhemus, evt] = PolhemusTapper('GetEvent', polhemus [, ForEvent=0]);
% - Retrieve next queued event received from the device in the struct 'evt'.
% If no new events are available and the optional 'waitForEvent' is set to
% 1, then the function will wait until at least one valid event becomes
% available and return that event. Otherwise it will return an empty struct,
% ie., evt = [] to signal that no new events are available.
%
% The following subfields are available in 'evt' if 'evt' is non-empty:
%
% evt.state = a value which encodes the new status of the tapper.
% 1 == tap, 2 == untap, 3 == top, 4 == untop, 5 == start.
%
% evt.index = index in the buffer array of polhemus corresponding to the
% event.
%
% evt.time  = Psychtoolbox GetSecs() timestamp of the time when the event
% was received from the device. The accuracy depends on the properties of
% the Polhemus Libery device and system load.
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
% tapPrediction = PolhemusTapper('GetTapPrediction', polhemus, polhemustt [, options]);
% - Return expected time of next tap (in seconds from now) given the 
% current state and position, as recorded from the last GetEvent.
% If no prediction is achieved, return NaN. A vector with 'options' might
% be provided if the prediction is offline.
%
%
%
% Note: Propenko library (MatLab interface for the Polhemus Liberty device) 
% installed on the SMLC Multisensory Perception LAB machine on day 18/01/11. 
% Since it is a shareware product, it expires the 18/02/11
% (but you can simply change the system date to make it run).

if nargin < 1
    help PolhemusTapper;
    return;
end

if nargin < 2
    polhemus = [];
end

if strcmpi(cmd, 'GetEvent')
    checkInputArgument('GetEvent');
    
    if length(varargin) < 1
        waitEvent = [];
    else
        waitEvent = varargin{1};
    end
    
    if isempty(waitEvent)
        waitEvent = 0;
    end
    
    % Start with empty event 'evt':
    evt = [];
    
    % Reset timestamp trouble flag:
    tTrouble = 0;
    skipped = 0;
    
    % Private function for reading new samples
    readNewSamples();
        
    % If no new data and everything has been already checked, return
    if polhemus.nextBuf <= polhemus.index
        varargout{1} = polhemus;
        varargout{2} = evt;        
        return;
    end 
    
    switch polhemus.tappingStatus
        case 0 % Finger nowhere, wait to enter down position
            for i = polhemus.index:(polhemus.nextBuf-1)
               if insideBottomBox(polhemus.bufData(i, :))
                    evt.time = polhemus.bufTimes(i);
                    evt.state = 5; % START
                    evt.trouble = tTrouble;
                    %polhemus.index = i+1;
                    %polhemus.tappingStatus = 2;
                    polhemus.index = i;
                    polhemus.tappingStatus = 1;
                    polhemus.lastTapPosition = polhemus.bufData(i, :);                        
                    display('Start');                    
                    break;
               end
            end
        
        
        case 1 % Finger down
            for i = polhemus.index:(polhemus.nextBuf-1)
                % Check if ascension has started
                totUpSteps = 3;
                if (i > totUpSteps)
                    upSteps = 0;
                    for j = 1:totUpSteps
                        if (polhemus.bufData(i-j, 2) < polhemus.bufData(i-j+1, 2)); 
                            upSteps = upSteps + 1; 
                        end
                    end                    
                    raiseThres = polhemus.downThres(2);
                    % raiseThres = polhemus.lastTapPosition(2) + 0.2*(polhemus.downThres(2)-polhemus.lastTapPosition(2));
                    
                    if upSteps == totUpSteps && polhemus.bufData(i, 2) > raiseThres
                        evt.time = polhemus.bufTimes(i-totUpSteps);
                        evt.index = i-totUpSteps;
                        evt.state = 2; % UNTAP
                        evt.trouble = tTrouble;
                        polhemus.index = i-totUpSteps+1;
                        polhemus.tappingStatus = 2;
                        display('Untap');
                        break;
                    end                    
                end
            end
                
        case 2 % Finger going up            
            % Scan all measurements
            for i = polhemus.index:(polhemus.nextBuf-1)
                % Check if sensor is inside top box
                if insideTopBox(polhemus.bufData(i, :))
                    evt.time = polhemus.bufTimes(i);
                    evt.index = i;
                    evt.state = 3; % TOP
                    evt.trouble = tTrouble;
                    polhemus.index = i+1;
                    polhemus.tappingStatus = 3;
                    display('Top');
                    break;
                end
            end
            
        case 3 % Finger up
            for i = polhemus.index:(polhemus.nextBuf-1)
                % Check for inversion
                if (i > 2 && polhemus.bufData(i, 2) > polhemus.upThres(2))
                    % Check for inversion
                    if (polhemus.bufData(i-2, 2) < polhemus.bufData(i-1, 2) && polhemus.bufData(i-1, 2) > polhemus.bufData(i, 2))
                        evt.time = polhemus.bufTimes(i-1);
                        evt.index = i-1;                        
                        evt.state = 4; % UNTOP                        
                        evt.trouble = tTrouble;
                        polhemus.index = i;
                        polhemus.tappingStatus = 4;
                        % polhemus.lastTapPosition = polhemus.bufData(i-1, :); ???                       
                        display('Untop');
                        break;                        
                    end                    
                end
            end            
            
        case 4 % Finger up going down
            for i = polhemus.index:(polhemus.nextBuf-1)
                % Check if position is below threshold
                if (i > 2 && polhemus.bufData(i, 2) < polhemus.downThres(2))
                    % Check for inversion
                    if (polhemus.bufData(i-2, 2) > polhemus.bufData(i-1, 2) && polhemus.bufData(i-1, 2) < polhemus.bufData(i, 2))
                        evt.time = polhemus.bufTimes(i-1);
                        evt.index = i-1;
                        evt.state = 1; % TAP                        
                        evt.trouble = tTrouble;
                        polhemus.index = i;
                        polhemus.tappingStatus = 1;
                        polhemus.lastTapPosition = polhemus.bufData(i-1, :);                        
                        display('Tap');
                        break;                        
                    end                    
                end
            end            
        
    end

    % No top-event found
    if isempty(evt); 
        polhemus.index = polhemus.nextBuf;
    else
        % If an event happened, bufferize it
        polhemus.bufEvents{length(polhemus.bufEvents)+1} = evt; 
    end
   
    % Return evt if any:
    varargout{1} = polhemus;
    varargout{2} = evt;
    
    return;
end

% Open connection to Polhemus:
if strcmpi(cmd, 'Open')
    if length(varargin) < 1; activeSensorNumber = [];
    else activeSensorNumber = varargin{1}; end
        
    if ~isnumeric(activeSensorNumber)
        error('PolhemusTapper: Open: Need to specify active sensor number!\n');
    end
    
    % Check if a non-empty polhemus handle is provided
    if isempty(polhemus); 
        nopolhemus = 1; 
        clear polhemus;
    else
        checkInputArgument('Open');
        nopolhemus = 0;
    end
    
    polhemus.h = actxserver('Prokopenko.Liberty');
    polhemus.h.connect;
    
    sensor_map = zeros(1, 32);
    sensor_map(activeSensorNumber) = 1;
    polhemus.h.sensor_map = sensor_map;
    polhemus.recording = 0;    
    polhemus.tappingStatus = 0; % Starting position is unknown

    % Initialize Polhemus temporary collected data
    if nopolhemus
        polhemus.bufTimes = [];
        polhemus.bufData = [];
        polhemus.bufEvents = []; 
        polhemus.bufSize = 0;
        polhemus.index = 0;
        polhemus.lastTapPosition = zeros(1, 6);
        polhemus.upThres = zeros(1, 6);
        polhemus.downThres = zeros(1, 6);
        % Perform a test?        
        fprintf('PolhemusTapper: Connection open!\n');
    end
    
    % Preheat GetSecs:
    GetSecs;
    
    varargout{1} = polhemus;
    return;
end

% Close connection to Polhemus tracker:
if strcmpi(cmd, 'Close')
    checkInputArgument('Close');
    
    if ~isempty(polhemus.h)    
        polhemus.h.disconnect;
        polhemus.h = [];
    end
    
    varargout{1} = polhemus;
    return;
end

% Stop tapping session:
if strcmpi(cmd, 'Stop')
    checkInputArgument('Stop');
    
    polhemus.recording = 0;
    polhemus.h.flush;
    polhemus.h.stop;

    varargout{1} = polhemus;
        
    return;
end


% Start tapping session:
if strcmpi(cmd, 'Start')
    checkInputArgument('Start');

    % Wait finger in resting position
    if size(varargin, 1) > 0; waitPosition = varargin{1};
    else waitPosition = 0; end
    
    % Standard buffer size
    if size(varargin, 2) > 1; bufSize = 5000;
    else bufSize = varargin{1}; end
    
   polhemus.tappingStatus = 0; % Starting position is unknown
   clear polhemus.bufTimes polhemus.bufData; % Clear buffers
   
   % Create and initialize data buffers
   polhemus.bufSize = bufSize;
   polhemus.nextBuf = 1;
   polhemus.index = 1;
   polhemus.bufTimes = zeros(bufSize, 1);
   polhemus.bufData = zeros(bufSize, 6);   
   polhemus.bufEvents = []; 
   polhemus.lastTapPosition = zeros(1, 6);
   polhemus.recording = 1;
       
   polhemus.h.start; % Start recording
   polhemus.h.flush;
           
   if waitPosition; waitPosition = WaitTillInPosition(); end
   
   varargout{1} = polhemus;
   varargout{2} = waitPosition;   
     
   return;
end

% Restart tapping session (flush data, set tappingStatus to zero):
if strcmpi(cmd, 'Restart')
    checkInputArgument('Restart');
        
    % Wait finger in resting position
    if size(varargin, 1) > 0; waitPosition = varargin{1};
    else waitPosition = 0; end
    
    polhemus.recording = 1;
    polhemus.tappingStatus = 0; % Starting position is unknown   
    polhemus.h.flush;
           
    if waitPosition; waitPosition = WaitTillInPosition(); end
   
    varargout{1} = polhemus;
    varargout{2} = waitPosition;   
         
    return;
end

% Calibrate Polhemus tapping device:
if strcmpi(cmd, 'Calibrate')
    checkInputArgument('Calibrate');
    
    if polhemus.recording
        error('PolhemusTapper: Calibrate: Cannot calibrate while recording data!');
    end
    
    polhemustt = [];
    downPosition = [];
    upPosition = [];

    % Check if calibration data are already provided
    if ~isempty(varargin)
        if isfield(varargin{1}, 'downPosition')
            polhemustt = varargin{1};
            downPosition = polhemustt.downPosition
            upPosition = polhemustt.upPosition
        elseif ~isempty(varargin{1})
            downPosition = varargin{1}(1, :)
            upPosition = varargin{1}(2, :)        
        end
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN CALIBRATION (DOWN AND UP POSITIONS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if isempty(downPosition)
        % Measure down position
        fprintf('Down: place your finger on the tapping pad and keep it in position. Press SPACE when you are ready.');
        KbWait();

       polhemus.h.start; % Start recording
       polhemus.h.flush;
       WaitSecs('YieldSecs', 1.0); 

       availableSample = polhemus.h.available_samples;
       if availableSample > 0
           sensorData = polhemus.h.read(availableSample);
           downPosition = mean(sensorData, 1)       
           % Check for trouble
       end
       polhemus.h.stop;
       clear sensorData;
    end

    if isempty(upPosition)
       % Measure up threshold
       fprintf('Up: raise your finger up to a comfortable angle and keep it in position. Press SPACE when you are ready.');
       KbWait();

       polhemus.h.start; % Start recording
       polhemus.h.flush;
       WaitSecs('YieldSecs', 1.0); 

       availableSample = polhemus.h.available_samples;
       if availableSample > 0
           sensorData = polhemus.h.read(availableSample);
           upPosition = mean(sensorData, 1)   
           % Check for trouble
       end
       polhemus.h.stop;
       clear sensorData;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END CALIBRATION (DOWN AND UP POSITIONS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   
   polhemus.downThres =  downPosition + (upPosition - downPosition)*0.25;
   polhemus.upThres =  downPosition + (upPosition - downPosition)*0.70;
       
   % Generate timetable
   if isempty(polhemustt)
    polhemustt = PolhemusNewTimeTable(downPosition, upPosition);
   end
    
   varargout{1} = polhemus;
   varargout{2} = polhemustt;
    return;
end

% Save the buffer (a short version of the tapper struct)
if strcmpi(cmd, 'SaveBuffer')
    checkInputArgument('SaveBuffer');
    ttemp.nextBuf = polhemus.nextBuf;
    ttemp.index = polhemus.index;
    ttemp.bufTimes = polhemus.bufTimes(1:(ttemp.nextBuf-1), :);
    ttemp.bufData = polhemus.bufData(1:(ttemp.nextBuf-1), :);   
    ttemp.bufEvents = polhemus.bufEvents;
    
    varargout{1} = ttemp;    
    return;
end

if strcmpi(cmd, 'Status')
    checkInputArgument('Status');
    
    % Return box status struct:
    % varargout{1} = box;
    
    return;
end

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
    dataList = polhemus.bufData;
    timesList = polhemus.bufTimes;
    
    % Update timetable data (proper events timing); it also updates
    % the 'evtTimesToTap' part of the timetable
    [tapList, timetable] = PolhemusNewTapList(evtList, dataList, timetable);
            
    % Update phase space data for all proper events
    for i = 1:size(tapList, 1)
        evtIndexes = tapList(i, :);
        evtTimes = [evtList{evtIndexes(1)}.time, evtList{evtIndexes(2)}.time, evtList{evtIndexes(3)}.time, evtList{evtIndexes(4)}.time];            
        tapTime = evtTimes(4);
        %evtTimesToTap = tapTime - [tapTime, evtTimes(1), evtTimes(2), evtTimes(3)];
        %timetable.evtTimesToTap = [timetable.evtTimesToTap; evtTimesToTap];        
        
        % Update times to tap in position matrices (ascension)
        %for j = evtList{i}.index:evtList{i+2}.index
        %    ypos = dataList(j, 2);
        %    yindex = round((ypos - timetable.positionBase)/timetable.positionStep);
        %    yindex = max(1, min(yindex, size(timetable.positionTimesToTap, 2)));
        %    timetable.positionTimesToTap(timetable.indexPosition(yindex, 1), yindex, 1) = tapTime - timesList(j);
        %    timetable.indexPosition(yindex, 1) = timetable.indexPosition(yindex, 1) + 1;
        %    if timetable.indexPosition(yindex, 1) > size(timetable.positionTimesToTap, 1); timetable.indexPosition(yindex, 1) = 1; end
        %end
        % Update times to tap in position matrices (descent)
        %for j = evtList{i+2}.index:evtList{i+3}.index
        %    ypos = dataList(j, 2);
        %    yindex = round((ypos - timetable.positionBase)/timetable.positionStep);
        %    yindex = max(1, min(yindex, size(timetable.positionTimesToTap, 2)));
        %    timetable.positionTimesToTap(timetable.indexPosition(yindex, 2), yindex, 2) = tapTime - timesList(j);
        %    timetable.indexPosition(yindex, 2) = timetable.indexPosition(yindex, 2) + 1;
        %    if timetable.indexPosition(yindex, 2) > size(timetable.positionTimesToTap, 1); timetable.indexPosition(yindex, 2) = 1; end
        %end

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
    % Not all phase space history is checked -- up to last tot entries
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
        
    % Update table of guessed times
    %timetable.guessedTappingTimes(:, 1) = NaN;
    %timetable.guessedTappingTimes(:, 2) = nanmean(timetable.positionTimesToTap(: ,:, 1), 1)';
    %timetable.guessedTappingTimes(:, 3) = nanmean(timetable.positionTimesToTap(: ,:, 1), 1)';
    %timetable.guessedTappingTimes(:, 4) = nanmean(timetable.positionTimesToTap(: ,:, 2), 1)';
    
    %for j = 1:4
    %    avgevtTimesToTap = nanmean(timetable.evtTimesToTap(:, j));
    %    timetable.guessedTappingTimes(:, j) = timetable.positionWeight(j)*timetable.guessedTappingTimes(:, j) + (1-timetable.positionWeight(j))*avgevtTimesToTap;
    %end
                
    varargout{1} = timetable;
    
    return;
end

% Compute the prediction error over the current buffered data
if strcmpi(cmd, 'ComputePredictionError')
    checkInputArgument('ComputePredictionError');
        
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
    [tapList, timetable] = PolhemusNewTapList(evtList, dataList, timetable);
        
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
    checkInputArgument('GetTapPrediction');
        
    % The optional argument should be a Polhemus timetable struct
    if size(varargin, 1) < 1;
        error('PolhemusTapper: ComputePredictionError: A valid Polhemus timetable must be provided!');
    else
        timetable = varargin{1};
    end
    
    tapPrediction = NaN;
    evtPrediction = NaN;
    phasespacePrediction = NaN;
    
    if length(varargin) > 1
        polhemus.nextBuf = varargin{2}(1)+1;
        elapsedTimeFromLastEvent = varargin{2}(2);
        tappingStatus = varargin{2}(3);
        ypos = polhemus.bufData(polhemus.nextBuf-1, 2);
        vpos = getVelocity(polhemus.nextBuf-1);
        [tapPrediction, evtPrediction, phasespacePrediction] = ...
        	getTapPrediction(ypos, vpos, elapsedTimeFromLastEvent, tappingStatus);       
    elseif polhemus.nextBuf > 1 && polhemus.tappingStatus > 1
        ypos = polhemus.bufData(polhemus.nextBuf-1, 2);
        vpos = getVelocity(polhemus.nextBuf-1);
        elapsedTimeFromLastEvent = GetSecs() - polhemus.bufEvents{length(polhemus.bufEvents)}.time;
        [tapPrediction, evtPrediction, phasespacePrediction] = ...
            getTapPrediction(ypos, vpos, elapsedTimeFromLastEvent, polhemus.tappingStatus);        
    end
        
    varargout{1} = tapPrediction;
    varargout{2} = evtPrediction;
    varargout{3} = phasespacePrediction;

    return;
end

% Invalid command!
error('PolhemusTapper: Invalid or unknown command specified!');

    % GETTAPPREDICTION Calculate the predicted time of the next tap (from
    % now, in seconds) given an y position, a velocity, the time from the 
    % last event and the current tapping status.
    %
    function [taptime, evtprediction, phasespaceprediction] = ...
            getTapPrediction(ypos, vpos, timeFromLastEvent, tappingStatus)
        % yindex = round((ypos - polhemus.timetable.positionBase)/polhemus.timetable.positionStep);
        % yindex = max(1, min(yindex, size(polhemus.timetable.guessedTappingTimes, 1)));
        % taptime = polhemus.timetable.guessedTappingTimes(yindex, tappingStatus);
        % avgevtTimesToTap = nanmean(polhemus.timetable.evtTimesToTap(:, tappingStatus)) - timeFromLastEvent;
        % taptime = polhemus.timetable.positionWeight(tappingStatus)*taptime + (1-polhemus.timetable.positionWeight(tappingStatus))*avgevtTimesToTap;
        
        phasespaceprediction = NaN;
        evtprediction = NaN;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EVENT-BASED PREDICTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(timetable.evtTimesToTap)
            taptime = NaN;
            return;
        end
        
        lookupmax = size(timetable.evtTimesToTap, 1);
        lookupmin = max(lookupmax - 10, 1);
        evtprediction = nanmean(timetable.evtTimesToTap(lookupmin:lookupmax, tappingStatus)) - timeFromLastEvent;
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BEGIN PHASE SPACE PREDICTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(timetable.phasespaceLookupTable)
            taptime = evtprediction;
            return;
        end
        
        % Calculate y index in phase space lookup table
        lookupRadius = 2;
        yindex = floor((ypos - timetable.phasespacemins(1))/timetable.phasespacesteps(1)) + 1;
        vindex = floor((vpos - timetable.phasespacemins(2))/timetable.phasespacesteps(2)) + 1;
        minbounds = [max(yindex - lookupRadius, 1), max(vindex - lookupRadius, 1)];
        maxbounds = [min(yindex + lookupRadius, size(timetable.phasespaceLookupTable, 2)), min(vindex + lookupRadius, size(timetable.phasespaceLookupTable{1}, 2))];
        
        phasespaceprediction = 0; % Estimated tapping time
        nf = 0; % Normalization factor
        nfdist = timetable.phasespacesteps.^2; % Normalization factor for distances (squared)        
        for yloop = minbounds(1):maxbounds(1)
            for vloop = minbounds(2):maxbounds(2)
                phasespacecell = timetable.phasespaceLookupTable{yloop}{vloop};
                for kloop = 1:size(phasespacecell, 1);
                    normdistance = ([ypos - phasespacecell(kloop, 1), vpos - phasespacecell(kloop, 2)].^2)./nfdist;
                    weight = exp(-sum(normdistance));
                    phasespaceprediction = phasespaceprediction + weight*phasespacecell(kloop, 3);
                    nf = nf + weight;
                end
            end
        end
        phasespaceprediction = phasespaceprediction / nf;
                
        %if isnan(taptime)
        %    yindex
        %    vindex
        %    minbounds
        %    maxbounds
        %   nf 
        %end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END PHASE SPACE PREDICTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isnan(phasespaceprediction)
            taptime = evtprediction;
            return;
        end
        
        taptime = phasespaceprediction;
        
    end

    % GETVELOCITY Calculate velocity from index position ipos. Returns
    % NaN if the velocity cannot be calculated meaningfully (e.g. error in
    % timing or insufficient buffering).
    function v = getVelocity(ipos)
        v = NaN;
        % Buffering needed
        if ipos < 3; return; end
        tdif = polhemus.bufTimes(ipos) - polhemus.bufTimes(ipos-2);
        % Each timestep should be about 4 ms; if not, the data are
        % unreliable and do not consider them
        if tdif < 0.006 || tdif > 0.010; return; end
        % Three-points derivative
        v = (1.5*polhemus.bufData(ipos, 2) - 2*polhemus.bufData(ipos-1, 2) + 0.5*polhemus.bufData(ipos-2, 2))/tdif;        
    end


    % READNEWSAMPLES Read data from stream buffer, save them into
    % 'polhemus' and return number of data points read.
    function availableSample = readNewSamples()
        startTime = GetSecs();
        availableSample = polhemus.h.available_samples;

        if availableSample > 0
            iNew = polhemus.nextBuf;
            sensorData = polhemus.h.read(availableSample);
            
            endTime = GetSecs();
            avgTime = 0.5*(startTime + endTime);

            % Check buffer capacity and double it if necessary
            if (iNew + availableSample - 1) > polhemus.bufSize
                polhemus.bufTimes = [polhemus.bufTimes; zeros(polhemus.bufSize, 1)];
                polhemus.bufData = [polhemus.bufData; zeros(polhemus.bufSize, 6)];
                polhemus.bufSize = size(polhemus.bufTimes, 1);
            end

            polhemus.bufTimes(iNew:(iNew + availableSample - 1), 1) = avgTime*ones(availableSample, 1);
            polhemus.bufData(iNew:(iNew+availableSample-1), :) = sensorData;
            polhemus.nextBuf = iNew + availableSample;

            % Check for trouble
        end        
    end

    % INSIDETOPBOX Boolean function, TRUE if y position is inside the top
    % (virtual) box in 3D space.
    function b = insideTopBox(position)
        if position(2) > polhemus.upThres(2); b = 1;
        else b = 0;
        end        
    end

    % INSIDEBOTTOMBOX Boolean function, TRUE if y position is inside the 
    % bottom (virtual) box in 3D space.
    function b = insideBottomBox(position)
        tapbox = [2.5 0 2.5; -2.5 -5 -2.5]; % Virtual box in 3D space
        if (all(position(1:3) < polhemus.downThres(1:3) + tapbox(1, :)) && all(position(1:3) > polhemus.downThres(1:3) + tapbox(2, :)))        
            b = 1;
        else b = 0;
        end        
    end

    % WAITTILLINPOSITION Stop execution of the program till the finger is 
    % in the resting position.
    function notInPosition = WaitTillInPosition()
        maxwaitingtime = 10;    % Wait till a maximum of 10 seconds
        tstart = GetSecs();
        notInPosition = 1;
            while (notInPosition)
                % Wait an instant in order not to fry the CPU
                WaitSecs('YieldSecs', 0.001);

                [polhemus, evt] = PolhemusTapper('GetEvent', polhemus);
                if ~isempty(evt) && evt.state == 5
                    notInPosition = 0;
                end

                tnow = GetSecs();                
                if (tnow - tstart > maxwaitingtime); return; end
            end
    end
% End of main function:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MACROS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % CHECKINPUTARGUMENT A macro for checking the validity of the tapper handle.
    function checkInputArgument(cmdString)
        if isempty(polhemus)
            error(['PolhemusTapper: ' cmdString, ': No "handle" for polhemus tapper device provided!']);
        end

        if length(polhemus) ~= 1
            error(['PolhemusTapper: ' cmdString, ': Passed argument is not a valid Polhemus "handle"!']);
        end
    end

end



