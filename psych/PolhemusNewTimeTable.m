% POLHEMUSNEWTIMETABLE Timetable struct constructor. Provide a downposition and a
% upposition (6-vectors) for the tapping device.
function newtt = PolhemusNewTimeTable(downPosition, upPosition)
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
   newtt.evtTimesToTap = [];
   newtt.phasespaceLookupTable = [];
   newtt.predictionError = [];
   newtt.predictionWeight = [0.9 0.9 0.5 0.];
end