function C = cellmatfun(fun,A,B)
%CELLMATFUN Apply a function to each cell/element of a (cell) array.
%   C = CELLMATFUN(FUN,A) applies the function specified by FUN to the
%   contents of each cell of cell array A, and returns the results in
%   the cell array C (or multiple cell arrays), where the (I,J,...)th cell 
%   contains the value FUN(A{I,J,...}, ...). If A is a numerical array,
%   C is simply the array FUN(A).
%
%   C = CELLMATFUN(FUN,A,B) applies the element-by-element binary 
%   operation specified by the function handle FUN to cell arrays A and B, 
%   with singleton expansion enabled. FUN can be any of the functions
%   accepted by BSXFUN. A and/or B can also be simple numerical arrays.
%   If both A and B are cell arrays, C is a cell array where the 
%   (I,J,...)th cell contains the value BSXFUN(FUN, A{I,J,...}, B{I,J,...}). 
%   If A is a cell array and B is a numerical array (or vice versa), C is 
%   a cell array  whose (I,J,...)th cell contains the value 
%   BSXFUN(FUN, A{I,J,...}, B) (or BSXFUN(FUN, A, B{I,J,...})).
%   If both A and B are numerical arrays, C is the array BSXFUN(FUN, A, B).
%
%   See also BSXFUN, CELLFUN, FUNCTION_HANDLE.

if nargin < 3
    if ~iscell(A)
        C = fun(A);
    else
        C = cellfun(fun, A, 'UniformOutput', false);
    end
else
    if ~iscell(A) && ~iscell(B)
        C = bsxfun(fun,A,B);
    elseif iscell(A) && ~iscell(B)
        C = cellfun(@(x) bsxfun(fun, x, B), A, 'UniformOutput', false);
    elseif ~iscell(A) && iscell(B)
        C = cellfun(@(x) bsxfun(fun, A, x), B, 'UniformOutput', false);
    else
        C = cellfun(@(x,y) bsxfun(fun, x, y), A, B, 'UniformOutput', false);
    end
end