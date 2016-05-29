function S = nansumall(X)
%SUM Sum of all elements, ignoring NaN values.
%   S = NANSUMALL(X) is the sum of all the elements of the array X,
%   computed after removing NaN values.
S = nansum(X(:));