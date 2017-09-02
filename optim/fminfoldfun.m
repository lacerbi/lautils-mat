function varargout = fminfoldfun(x,varargin)
%FMINFOLDFUN Wrapper function for simultaneous minimization of folds.
%
%   FMINFOLDFUN([],FUN,MASK,TESTMASK) initialize FMINFOLD wrapper for objective 
%   function FUN. FUN must be either a string with the function name or a 
%   function handle for the objective function F(X). F takes a parameter 
%   vector X as input and returns a numeric array as output. Each element
%   of the output represents the function value at a given data point.
%   FUN can also return as a second output a matrix array of gradients per
%   data point (each column is the gradient for a data point).
%   MASK is a logical matrix where each column represents a trial and each 
%   row represents a fold (subset of data points). Typically the first row
%   of MASK is all TRUE (the full data set), whereas other rows represent
%   subsets of the data (e.g., training folds in K-fold cross-validation).
%   TESTMASK is a logical matrix similar to MASK.
%
%   FMINFOLDFUN([],FUN,MASK,TESTMASK,IFOLD,N) also sets the currently queried fold
%   to IFOLD (default IFOLD=1), and the cache size (number of stored function
%   values) to N (default N=1e4).
%
%   FMINFOLDFUN([],[],[],[],IFOLD) sets the currently queried fold to IFOLD.
%
%   FVAL = FMINFOLDFUN(X) returns the value of the logged function at X,
%   where X is a numeric array.
%
%   [FVAL,DF] = FMINFOLDFUN(X) also returns the gradient of the logged
%   function at X (summed over trials of the current fold). The logged
%   function must return the gradient as second output.
%
%   [...] = FMINFOLDFUN(X,ARG1,ARG2,...) passes a series of additional 
%   arguments to the call of the logged function.
%
%   FOLDLOG = FMINFOLDFUN() returns the function log structure.
%

%--------------------------------------------------------------------------
% To be used under the terms of the GNU General Public License 
% (http://www.gnu.org/copyleft/gpl.html).
%
%   Author (copyright): Luigi Acerbi, 2017
%   e-mail: luigi.acerbi@{gmail.com,nyu.edu}
%   URL: http://luigiacerbi.com
%   Version: Aug/23/2017
%   Code repository: https://github.com/lacerbi/bads
%--------------------------------------------------------------------------

persistent foldlog;     % Record log of function calls

%% Case 1: No arguments ==> return function log
if nargin < 1
    varargout = {foldlog};

    %if isstruct(x)  % Swapped variable order
    %    temp = probstruct;
    %    probstruct = x;
    %    x = temp;
    %end

%% Case 2: Empty X, no output, pass function name or handle ==> Re-initialize log struct
elseif isempty(x) && nargout == 0 && ~isempty(varargin{1}) && ...
        (ischar(varargin{1}) || isa(varargin{1},'function_handle'))
    
    % Reset log struct
    foldlog = [];
    
    % Read inputs
    fun = varargin{1};
    mask = varargin{2};
    testmask = varargin{3};
    if nargin > 4; iFold = varargin{4}; else, iFold = 1; end
    if nargin > 5; N = varargin{5}; else, N = []; end

    % Store function name and function handle
    if ischar(fun)
        foldlog.FuncName = fun;
        foldlog.FuncHandle = str2func(fun);
    else
        foldlog.FuncName = func2str(fun);
        foldlog.FuncHandle = fun;        
    end
    
    % Number of stored function values
    if isempty(N); N = 1e4; end
    
    % Size of function output arguments
    foldlog.OutSize = size(mask,2);
    
    % Fold mask
    if any(all(mask == 0,2))
        error('At least one row of MASK contains all zeroes.');
    end
    foldlog.Mask = mask;
    nfolds = size(mask,1);
        
    % Test mask (used for testing)
    if isempty(testmask)
        testmask = ~mask;   % By default, take complement
        idx = find(all(testmask,2));
        if ~isempty(idx); testmask(idx,:) = true; end
    end
    if any(size(mask) ~= size(testmask))
        error('MASK and TESTMASK should have the same size.');
    end
    foldlog.TestMask = testmask;
    
    % Current fold
    testifold(iFold,nfolds);
    foldlog.iFold = iFold;    
    
    foldlog.Clock = tic;                    % Time of initialization
    foldlog.FuncCount = zeros(1,nfolds);    % Function count per fold
    foldlog.X = [];
    foldlog.Fval = [];
    foldlog.FvalTest = [];
    
    % Best results for each fold
    foldlog.BestX = [];
    foldlog.BestFval = Inf(nfolds,1);
    foldlog.TestFvalAtBestX = Inf(nfolds,1);
    foldlog.NewSearchFlag = false(nfolds,1);
    
    foldlog.N = N;           % Size of X matrix, temporary
    foldlog.last = 0;        % Last entry
    
%% Case 3: Empty X, no output, non-empty fifth argument, set IFOLD
elseif isempty(x) && nargout == 0 && nargin == 5 ...
        && isempty(varargin{1}) && isempty(varargin{2}) && isempty(varargin{3})
    
    if isempty(foldlog)
        error('FMINFOLDFUN has not been initialized.');        
    end
    
    iFold = varargin{4};
    nfolds = size(foldlog.Mask,1);
    testifold(iFold,nfolds);    
    foldlog.iFold = iFold;
    
%% Case 4: Pass X, get output ==> Standard function call    
else

    gradient_flag = (nargout > 1);    % Is the gradient requested?
    
    % You need to initialize FMINFOLDFUN first
    if isempty(foldlog)
        error('FMINFOLDFUN has not been initialized.');
    end
    
    if isempty(foldlog.iFold)
        error('IFOLD not specified. Use FMINFOLD([],[],IFOLD) to select a fold to evaluate.');
    end

    % Call function
    func = foldlog.FuncHandle;

    % Check if need to pass additional arguments
    try        
        if nargin > 1
            if ~gradient_flag
                fval = func(x,varargin{:});
            else
                [fval,df] = func(x,varargin{:});    % Get gradient as well
            end
        else
            if ~gradient_flag
                fval = func(x);
            else
                [fval,df] = func(x);                % Get gradient as well
            end
        end
    catch funcError
        warning(['Error in executing the logged function ''' foldlog.FuncName ''' with input:']);
        x
        rethrow(funcError);
    end
    
    % Check that returned gradient has the right shape
    if gradient_flag
        if size(df,1) ~= numel(x) || size(df,2) ~= foldlog.OutSize
            error('The 2nd output argument of the logged function should be a matrix of the gradient per data point (each column is a different data point).');
        end
    end

    % Update function counter for the current fold
    foldlog.FuncCount(foldlog.iFold) = foldlog.FuncCount(foldlog.iFold) + 1;

    % Record log
    if isempty(foldlog.X)
        foldlog.X = NaN(foldlog.N,numel(x));
        foldlog.Fval = NaN(foldlog.N,size(foldlog.Mask,1));
        foldlog.FvalTest = NaN(foldlog.N,size(foldlog.TestMask,1));
    end
    foldlog.last = max(1,mod(foldlog.last+1, foldlog.N+1));
    foldlog.X(foldlog.last,:) = x;
    
    % Compute function value for each fold
    Fval_fold = sum(bsxfun(@times,foldlog.Mask,fval(:)'),2);
    
    % Compute function value for each test fold
    FvalTest_fold = sum(bsxfun(@times,foldlog.TestMask,fval(:)'),2);    
    
    % Store function value computed on each fold
    foldlog.Fval(foldlog.last,:) = Fval_fold;
    foldlog.FvalTest(foldlog.last,:) = FvalTest_fold;
    
    % Update best values for each fold
    newbestflag = Fval_fold < foldlog.BestFval;
    for idx = find(newbestflag)
        foldlog.BestX(idx,1:numel(x)) = repmat(x,[numel(idx),1]);
        foldlog.BestFval(idx) = Fval_fold(idx);
        foldlog.TestFvalAtBestX(idx) = FvalTest_fold(idx);
    end
    
    % Update flags for starting new searches
    foldlog.NewSearchFlag = foldlog.NewSearchFlag | newbestflag;
    foldlog.NewSearchFlag(foldlog.iFold) = false;
    
    % First output, function value of current fold
    varargout{1} = foldlog.Fval(foldlog.last,foldlog.iFold);
    
    % Second returned output is assumed to be the gradient
    % (summed over trials, only for the current fold)
    if gradient_flag
        varargout{2} = sum(bsxfun(@times,foldlog.Mask(foldlog.iFold,:),df),2);
    end
    
end

end

%--------------------------------------------------------------------------
function testifold(iFold,nfolds)
%TESTIFOLD Test validity of IFOLD input.

if ~isscalar(iFold) || isempty(iFold) || ...
        iFold < 1 || iFold > nfolds || round(iFold) ~= iFold
    error('IFOLD should be an integer between 1 and NFOLDS.');
end

end