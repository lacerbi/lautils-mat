function s = numarray2str(x,format,sep,left,right)
%NUMARRAY2STR Convert a numerical array into a bracketed string
%
%   S = NUMARRAY2STR(X) converts numerical array X into a string S of the
%   form '[X(1),X(2),...,X(end)]'.
%
%   S = NUMARRAY2STR(X,N) converts the elements of X into a string 
%   representation with a maximum N digits of precision.  The default 
%   number of digits is based on the magnitude of the elements of X.
%
%   S = NUMARRAY2STR(X,FORMAT) uses the format string FORMAT (see SPRINTF 
%   for details).
%
%   S = NUMARRAY2STR(A,FORMAT,SEP) uses char array SEP to separate distinct 
%   elements of X in S. The default value is comma ','.
%
%   S = NUMARRAY2STR(A,FORMAT,SEP,LEFT,RIGHT) uses char array LEFT to open
%   the bracket, and char array RIGHT to close the bracket. The default 
%   values are square brackets, respectively '[' and ']'.
%
%   See also NUM2STR, SPRINTF.

assert(isnumeric(x), 'X needs to be a numerical array.');

if nargin < 2; format = []; end
if nargin < 3; sep = ','; end
if nargin < 4; left = '['; end
if nargin < 5; right = ']'; end

if isempty(format); func = @(x_) num2str(x_);
else func = @(x_) num2str(x_,format); end

if isempty(x)
    s = ''; 
else
    s = left;
    for i = 1:numel(x)-1; s = [s,func(x(i)),sep]; end
    s = [s,func(x(end)),right];
end