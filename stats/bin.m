function [b,bb,x] = bin(y, x)
%HIST  Histogram.
%   B = BIN(Y,X) bins the elements of Y into containers with centers
%   specified by X and return an array of bins.  The first bin includes
%   data between -inf and the first center and the last bin
%   includes data between the last bin and inf. 

if size(y,2) > 1; y2 = y'; else y2 = y; end
if size(x,1) > 1; x2 = x'; else x2 = x; end

[x2 ord] = sort(x2); % Sort x

m = size(y2,1);
k = size(x2,2);

d = abs(y2*ones(1, k) - ones(m, 1)*x2);
[~, index] = min(d, [], 2);

for i = 1:k
    bb{ord(i)} = y(index == i);
end

b = ord(index);



