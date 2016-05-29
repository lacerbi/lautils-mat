function varargout = fminfold_wrap(x,varargin)
%FMINFOLD_WRAP Auxiliary wrapper function for FMINFOLD (do not call directly).
%
%   FMINFOLD_WRAP([],FUN,MASK,IFOLD) initialize wrapper for function FUN. 
%   FUN must be either a string with the function name or a function handle 
%   for function F(X). F takes a numeric array as input and returns a 
%   numeric array as output (one element per data point). MASK is a logical
%   matrix where each column represents a trial and each row represents a
%   fold. IFOLD is the number of the currently queried fold.
%
%   FMINFOLD_WRAP([],FUN,MASK,IFOLD,N) stores up to N function calls 
%   (default N=1e4).
%
%   FVAL = FMINFOLD_WRAP(X) returns the value of the logged function at X,
%   where X is a numeric array.
%
%   [FVAL,DF] = FMINFOLD_WRAP(X) also returns the gradient of the logged
%   function at X (summed over trials of the current fold). The logged
%   function must return the gradient as second output.
%
%   [...] = FMINFOLD_WRAP(X,ARG1,ARG2,...) passes a series of additional 
%   arguments to the call of the logged function.
%
%   FOLDLOG = FMINFOLD_WRAP() returns the function log structure.
%
% 
%   Author: Luigi Acerbi
%   Version: Mar/03/2016
%
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
elseif isempty(x) && nargout == 0 && ...
        (ischar(varargin{1}) || isa(varargin{1},'function_handle'))
    foldlog = [];
    fun = varargin{1};
    mask = varargin{2};
    iFold = varargin{3};
    if nargin > 4; N = varargin{4}; else N = []; end

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
    foldlog.Mask = mask;
    
    % Current fold
    foldlog.iFold = iFold;
    
    foldlog.Clock = tic;     % Time of initialization
    foldlog.FuncCount = 0;
    foldlog.X = [];
    foldlog.Y = [];
    
    foldlog.N = N;           % Size of X matrix, temporary
    foldlog.last = 0;        % Last entry
    
%% Case 3: Pass X, get output ==> Standard function call    
else

    % You need to initialize FMINFOLD_WRAP first
    if isempty(foldlog)
        error('FMINFOLD_WRAP has not been initialized.');
    end

    % Call function
    func = foldlog.FuncHandle;

    % if isfield(probstruct,'octstruct'); x = optimct(x,probstruct.octstruct,1); end

    % Check if need to pass probstruct
    try
        
        if nargin > 1
            if nargout == 1
                fval = func(x,varargin{:});
            else
                [fval,df] = func(x,varargin{:});    % Get gradient as well
            end            
        else
            if nargout == 1
                fval = func(x);
            else
                [fval,df] = func(x);                % Get gradient as well
            end
        end
    catch
        warning(['Error in executing the logged function ''' foldlog.FuncName ''' with input:']);
        x
        error('Execution interrupted.');
        % fval = NaN(1,foldlog.OutSize);
    end

    foldlog.FuncCount = foldlog.FuncCount + 1;

    % Record log
    if isempty(foldlog.X)
        foldlog.X = NaN(foldlog.N,numel(x));
        foldlog.Y = NaN(foldlog.N,size(foldlog.Mask,1));
    end

    foldlog.last = max(1,mod(foldlog.last+1, foldlog.N+1));
    foldlog.X(foldlog.last,:) = x;
    
    foldlog.Y(foldlog.last,:) = sum(bsxfun(@times,foldlog.Mask,fval(:)'),2);
    
    varargout{1} = foldlog.Y(foldlog.last,foldlog.iFold);
    
    % Second returned output is assumed to be the gradient 
    % (summed over trials, only for the current fold)
    if nargout > 1
        varargout{2} = sum(bsxfun(@times,foldlog.Mask(foldlog.iFold,:)',df),1);
    end
    
end

end