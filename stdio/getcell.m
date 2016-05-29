function E = getcell(X,n)
%GETCELL Return n-th element of cell array.
%
%   E = GETCELL(X) returns the first element of cell array X, or X itself
%   if X is not a cell array.
%
%   E = GETCELL(X,N) returns the N-th element of cell array X, or the empty
%   matrix if N > length(X). If X is not a cell array, returns X if N = 1, 
%   or the empty matrix otherwise.

if nargin < 2 || isempty(n); n = 1; end

if iscell(X) && n <= length(X)
    E = X{n};
elseif n == 1
    E = X;
else
    E = [];
end

end