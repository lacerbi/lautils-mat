function M = gridsum(DIM,F,GRID)
%GRIDSUM computes the sum value of a grid variable along the DIM-th
% dimension over regular grid GRID and weighting function FUN.
%                 
%   M = GRIDSUM(DIM,F,GRID) computes the sum M over grid GRID of FUN 
%   weighting by the grid values along dimension DIM. 
%   values of function FUN at points G. 
%   F is a D-dimensional weighting matrix. GRID is a structure of D vectors 
%   specifying grid values for each dimension.
%
%   GRIDSUM works up to D = 5.
%
%   See also fmingrid, grideval, gridmean.

if nargin < 3
    error('GRIDSUM requires at least a sum direction DIM, a weighting matrix FUN and a grid structure GRID. Digit ''help gridmean'' for instructions.');
end

% Number of dimensions, and length of each grid dimension
D = length(GRID);
for i = 1:D; DL(i) = length(GRID{i}); end

% Perform the grid search
index = ones(1, D);
X = zeros(1, D);
for i = 1:D; X(i) = GRID{i}(index(i)); end

M = 0;

cycleOn = 1;
while cycleOn    
    switch D
        case 1;
            f = F(index(1));
        case 2;
            f = F(index(1), index(2));
        case 3;
            f = F(index(1), index(2), index(3));
        case 4;
            f = F(index(1), index(2), index(3), index(4));            
        case 5;
            f = F(index(1), index(2), index(3), index(4), index(5));            
    end
    
    M = M + f*GRID{DIM}(index(DIM));
        
    % Increase the index -- the cycle is off if all dimensions have a carry
    cycleOn = 0;
    for i = 1:D
        index(i) = index(i) + 1;
        if index(i) > DL(i); 
            index(i) = 1;
            X(i) = GRID{i}(1);
        else
            X(i) = GRID{i}(index(i));
            cycleOn = 1;
            break;
        end        
    end
    
end

end