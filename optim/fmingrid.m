function [x,fval,exitflag,output] = fmingrid(funfcn,varargin)
%FMINGRID Multidimensional minimization via grid search.
%   X = FMINGRID(FUN,LB,UB,N) finds a global minimum X of the function FUN 
%   in a grid defined by lower bounds LB and upper bounds LB. FUN accepts a 
%   row vector X as input and returns a scalar function value F evaluated 
%   at X. LB and UB are row vectors of the same length as the input 
%   accepted by FUN. If N is a scalar, FMINGRID creates a grid with N 
%   points per dimension. If N is a vector, FMINGRID creates a grid with 
%   N(1) points along the 1st dimension, N(2) points along the 2nd 
%   dimension, and so on.
%
%   X = FMINGRID(FUN,GVEC) uses vectors contained in GVEC to build a grid.
%   GVEC is a cell array where GVEC{1} contains the vector of points for 
%   the 1st dimension of the grid, GVEC{2} contains the vector of points 
%   for the 2nd dimension of the grid, and so on. The vectors in GVEC need
%   not be equally spaced.
% 
%   X = FMINGRID(...,NOUT) returns NOUT grid points with the lowest values 
%   of FUN, ordered from the minimum to the maximum. Default value is 
%   NOUT = 1. Use NOUT = Inf to return all points.
%
%   X = FMINGRID(...,NOUT,1) performs a vectorized grid search (VECFLAG = 1).
%   For this, FUN needs to accept a matrix X as input, where each row of X 
%   is a vector to be evaluated, and returns a column vector F such that
%   F(n) = FUN(X(n,:)). A vectorized grid search can be orders of magnitude
%   faster than a non-vectorized search based on a for loop. By default the 
%   search is not vectorized (VECFLAG = 0).
%
%   X = FMINGRID(...,NOUT,VECFLAG,FCON) subjects the minimization to the 
%   constraints defined in FCON. The function FCON accepts a matrix X whose
%   rows are grid vectors and returns the vector C. The n-th element of C
%   is 1 if the n-th vector X(n,:) respects the constraints, 0 otherwise.
%   Points that do not satisfy the constraints are removed from the grid.
%
%   [X,FVAL] = FMINGRID(...) returns the value of the objective function,
%   described in FUN, at X. If NOUT > 1, FVAL is a colum vector, ordered
%   from the lowest to the highest reported value (same as X).
%
%   [X,FVAL,EXITFLAG] = FMINGRID(...) returns an EXITFLAG that describes 
%   the exit condition of FMINGRID. Currently EXITFLAG is always 0. This
%   output is used for compatibility with other minimization functions.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINGRID(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%        
%   FMINGRID is a very coarse minimization method. It is generally 
%   recommended to use more efficient methods.
%
%   See also FMINCON, FMINSEARCH, FMINUNC.

% Author: Luigi Acerbi
% e-mail: luigi.acerbi@gmail.com
% Release: 0.9
% Release date: 17/Jul/2015

if nargin < 2
    error('fmingrid:NotEnoughInputs','FMINGRID requires at least two input arguments.');
end

if isnumeric(varargin{1})   % Provide grid lower/upper bounds
    if nargin < 4
        error('fmingrid:NotEnoughInputs','FMINGRID requires at least four input arguments with lower/upper bounds.');
    end
    lb = varargin{1};
    ub = varargin{2};
    n = varargin{3};
    % Check for non-double inputs
    if ~isa(lb,'double') || ~isa(ub,'double')
      error('fmingrid:NonDoubleInput','FMINGRID accepts boundary inputs LB and UB only of data type double.')
    end

    % Default values
    if nargin < 5 || isempty(varargin{4}); nout = 1; else nout = varargin{4}; end
    if nargin < 6 || isempty(varargin{5}); vecflag = 0; else vecflag = varargin{5}; end
    if nargin < 7; fcon = []; else fcon = varargin{6}; end

    % Check that LB and UB are vectors of the same size
    if ~isvector(lb) || ~isvector(ub) || numel(lb) ~= numel(ub)
      error('fmingrid:BoundsDifferentSize','FMINGRID bounds LB and UB should be vectors with the same size.')    
    end
    
    % Number of dimensions
    D = length(lb);

    % Expand n to a vector if scalar
    if isscalar(n); n = n*ones(1,D); end

    % Create grid vectors
    for d = 1:D; gvec{d} = linspace(lb(d),ub(d),n(d)); end
else
    gvec = varargin{2};

    % Default values
    if nargin < 3 || isempty(varargin{2}); nout = 1; else nout = varargin{2}; end
    if nargin < 4 || isempty(varargin{3}); vecflag = 0; else vecflag = varargin{3}; end
    if nargin < 5; fcon = []; else fcon = varargin{4}; end
    
    % Number of dimensions
    D = length(gvec);
    
    % Convert all grid vectors to row vectors
    for d = 1:D; gvec{d} = gvec{d}(:)'; end
end

% This output is just to be compatible with other minimization functions
exitflag = 0;   

% Create list
xlist = cartvec(gvec{:})';

% Remove grid points that do not satisfy constraints
if ~isempty(fcon)
    if ~isa(fcon,'function_handle')
        error('fmingrid:ConstraintsNotAFunction','FMINGRID constraints FCON should be a function handle.');
    end
    xlist = xlist(logical(fcon(xlist)),:); 
end

ntot = size(xlist,1);
fvallist = NaN(ntot,1);
nout = min(ntot,nout); % Maximum length of output

if vecflag     % Vectorized function
    fvallist = funfcn(xlist);
    output.iterations = 1;
    output.funcCount = 1;
else
    for i = 1:ntot
        fvallist(i) = funfcn(xlist(i,:));
    end
    output.iterations = ntot;
    output.funcCount = ntot;
end

if nout == 1
    [fval,index] = min(fvallist);
    x = xlist(index,:);
else
    [fvallist,index] = sort(fvallist,1,'ascend');
    fval = fvallist(1:nout);
    x = xlist(index(1:nout),:);
end

if nargout > 3
    output.algorithm = 'Grid search';
    output.message = 'Optimization terminated: full grid explored.';
end

end % FMINGRID

%--------------------------------------------------------------------------
function y = cartvec(varargin)
%CARTVEC Cartesian product of row vectors.
%
% Inspired by Mark Beale's COMBVEC.

if nargin < 1
    y = [];
else
    y = varargin{1};
    for i=2:length(varargin)
        z = varargin{i}(:)';
        y = [block_copy(y,size(z,2)); interleave_copy(z,size(y,2))];
    end
end

    function b = block_copy(m,n)
        [mr,mc] = size(m);
        b = zeros(mr,mc*n);
        ind = 1:mc;
        for j=(0:(n-1))*mc; b(:,ind+j) = m; end
    end

    function b = interleave_copy(m,n)
        mc = size(m,2);
        b = zeros(n,mc);
        for k=0:(n-1); b(1+k,:) = m; end
        b = reshape(b,1,n*mc);
    end
end
