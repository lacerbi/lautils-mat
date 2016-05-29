function [X,FVAL,EXITFLAG,OUTPUT] = myfmingrid(FUN,GRID,NM,LB,UB,options,fminfun)
%FMINGRID performs a grid search of the minimum of a function of several variables.
%              
%   FMINGRID minimizes a function F of several variables by first
%   performing a brute-force sweep over a provided grid. The best NM minima
%   are kept and a standard FMINCON minimization is performed starting from 
%   them.
%   
%   X = FMINGRID(FUN,GRID [, NM = 10]) finds a minimum X to the function FUN. 
%   FUN accepts input X of dimension D and returns a scalar function value 
%   F evaluated at X. GRID is a structure of D vectors specifying grid 
%   values for each dimension. NM is the number of saved points after the
%   grid search, where a standard minimization is performed. If NM is zero
%   do not perform additional minimizatin.
% 
%   X = FMINGRID(FUN,GRID, NM, LB, UB) as above, where the standard
%   minimization is subject to the lower bounds LB and upper bounds UB.
%
%   X = FMINGRID(FUN,GRID, NM, LB, UB, OPTIONS) as above, where the standard
%   minimization is passed the option list OPTIONS. 
%
%   X = FMINGRID(FUN,GRID, NM, LB, UB, OPTIONS, FMINFUN) use minimization
%   algorithm specified in FMINFUN; it can be either 'fminsearchbnd' or 
%   'fmincon' (by default, FMINGRID uses fmincon).
%
%   See also fminsearchbnd, fmincon.

% Default values
if ~exist('NM', 'var'); NM = [] ; end
if isempty(NM); NM = 10; end
if ~exist('LB', 'var'); LB = [] ; end
if ~exist('UB', 'var'); UB = [] ; end
if ~exist('options', 'var'); options = [] ; end
if ~exist('fminfun', 'var'); fminfun = [] ; end
if isempty(fminfun); fminfun = 'fmincon'; end

% Minimization function
if strcmpi(fminfun, 'fminsearchbnd'); altfun = 1;
elseif strcmpi(fminfun, 'fmincon'); altfun = 0;
else error('Unknown function specified in FMINFUN.');
end

if nargin < 2
    error('FMINGRID requires at least a function handle FUN and a grid structure GRID. Digit ''help fmingrid'' for instructions.');
end

verbose = 0;

% Look for verbosity parameter
if ~isempty(options)
    if strcmpi(options.Display, 'off'); verbose = 0; else verbose = 1; end
end

% Number of dimensions, and length of each grid dimension
D = length(GRID);
for i = 1:D; DL(i) = length(GRID{i}); end

% Initialize best values
xvalues = ones(max(NM, 1), D)*NaN;
fvalues = ones(max(NM, 1), 1)*Inf;

% Perform the grid search
index = ones(1, D);
X = zeros(1, D);

for i = 1:D; X(i) = GRID{i}(index(i)); end

cycleOn = 1;
while cycleOn
    f = FUN(X);
    
    % Check if the function value at the point is lower than the worst    
    [fmax, maxindex] = max(fvalues);    
    if f < fmax
        fvalues(maxindex) = f;
        xvalues(maxindex, :) = X;
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

% Trim the unused solutions
xnans = any(isnan(xvalues), 2);
xvalues(xnans, :) = [];
fvalues(xnans, :) = [];

exitflag = zeros(1, length(fvalues));
output = [];

% Minimize with a standard fmincon the found best values
if NM > 0
    for i = 1:length(fvalues)
        if verbose; display(xvalues(i, :)); end
        if altfun
            [xvalues(i, :) fvalues(i) exitflag(i)] = fminsearchbnd(FUN, xvalues(i, :), LB, UB, options);
        else
            [xvalues(i, :) fvalues(i) exitflag(i)] = fmincon(FUN, xvalues(i, :), [], [], [], [], LB, UB, [], options);
        end
        if verbose; display(xvalues(i, :)); display(fvalues(i)); end
    end
end

% Find the minimum among the minima
[FVAL ival] = min(fvalues);
X = xvalues(ival, :);
if ~isempty(exitflag); EXITFLAG = exitflag(ival); end
if ~isempty(output); OUTPUT = output{ival}; end

end