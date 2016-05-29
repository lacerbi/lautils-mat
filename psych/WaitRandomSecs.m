% WAITRANDOMSECS Wait for a random time interval.
%   WAKEUP=WAITRANDOMSECS(INTERVAL) waits between INTERVAL(1) and 
%   INTERVAL(2) seconds.
%
% See also WAITSECS.
function WaitRandomSecs(interval)
WaitSecs(interval(1) + rand()*(interval(2)-interval(1)));
end