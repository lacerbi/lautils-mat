function [Z,GRID] = grideval(FUN,GRID)
%GRIDEVAL computes the values taken by a function FUN of several variables
% on a regular grid.
%                 
%   Z = GRIDEVAL(FUN,GRID) computes the values of function FUN at points G. 
%   FUN accepts input X of dimension D and returns a scalar function value 
%   F evaluated at X. GRID is a structure of D vectors specifying grid 
%   values for each dimension. Z is a D-dimensional matrix holding the
%   values of FUN at the point in GRID.
%
%   GRIDEVAL works up to D = 5. Performing a grid evaluation for greater
%   number of dimensions is arguably pointless.
%
%   See also fmingrid, gridmean.

if nargin < 2
    error('GRIDEVAL requires at least a function handle FUN and a grid structure GRID. Digit ''help grideval'' for instructions.');
end

% Number of dimensions, and length of each grid dimension
D = length(GRID);
for i = 1:D; DL(i) = length(GRID{i}); end

% Perform the grid search
index = ones(1, D);
X = zeros(1, D);
for i = 1:D; X(i) = GRID{i}(index(i)); end

switch D
    case 1;
        Z = zeros(length(GRID{1}), 1);
    case 2;
        Z = zeros(length(GRID{1}), length(GRID{2}));
    case 3;
        Z = zeros(length(GRID{1}), length(GRID{2}), length(GRID{3}));
    case 4;
        Z = zeros(length(GRID{1}), length(GRID{2}), length(GRID{3}), length(GRID{4}));
    case 5;
        Z = zeros(length(GRID{1}), length(GRID{2}), length(GRID{3}), length(GRID{4}), length(GRID{5}));
    otherwise
        error('GRIDEVAL for the moment supports only up to D = 5. Do you really need more?');            
end

% Create grid of zeros
G = Z;

cycleOn = 1;
while cycleOn
    f = FUN(X);
    
    switch D
        case 1;
            Z(index(1)) = f;
        case 2;
            Z(index(1), index(2)) = f;
        case 3;
            Z(index(1), index(2), index(3)) = f;
        case 4;
            Z(index(1), index(2), index(3), index(4)) = f;            
        case 5;
            Z(index(1), index(2), index(3), index(4), index(5)) = f;            
    end
        
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