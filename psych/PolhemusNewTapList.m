% POLHEMUSNEWTAPLIST Build a list of valid taps from an event list and data list 
% (and optionally a timetable)
function [tapList, timetable] = PolhemusNewTapList(evtList, dataList, timetable)

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
        timetable = PolhemusNewTimeTable(downPosition, upPosition);   
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
            %if ~isnan(timetable.evtTimesToTap(4, 1))
                %meanTimesToTap = meannan(timetable.evtTimesToTap, 1);
                %stdTimesToTap = max(stdnan(timetable.evtTimesToTap, 1), 0.020);
                %distanceTimesToTap = abs(evtTimesToTap - meanTimesToTap);
                % If any measurement is more than FIVE std away, do not use
                % tap -- it is anomalous
                
                %if any(distanceTimesToTap - 5*stdTimesToTap > 0)
                %    display('anomalous tap');
                %    continue;
                %end
            %end
            
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
