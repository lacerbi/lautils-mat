function M = gridmean(DIM,F,GRID)
%GRIDMEAN computes the mean value of a grid variable along the DIM-th
% dimension over regular grid GRID and weighting function FUN.
%                 
%   M = GRIDMEAN(DIM,F,GRID) computes the mean M over grid GRID of FUN 
%   weighting by the grid values along dimension DIM. 
%   values of function FUN at points G. 
%   F is a D-dimensional weighting matrix. GRID is a structure of D vectors 
%   specifying grid values for each dimension.
%
%   GRIDMEAN works up to D = 5.
%
%   See also fmingrid, grideval, gridsum.

if nargin < 3
    error('GRIDMEAN requires at least a sum direction DIM, a weighting matrix FUN and a grid structure GRID. Digit ''help gridmean'' for instructions.');
end

M = gridsum(DIM, F, GRID)/numel(F);

end