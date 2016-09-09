function [dy,y0] = fgrad(fun,x0,method,varargin)
%FGRAD Estimation of function gradient via finite differences
%
%   DY = FGRAD(FUN,X0) estimates the gradient of function FUN at X0 via 
%   finite differences. FUN is a function handle.  FUN accepts input X and 
%   returns either a scalar or a row vector function value F evaluated at 
%   X. If F is a row vector, the gradient is computed separately for 
%   each element of F. X0 can be a scalar or a row vector. 
%
%   DY is a row vector of partial derivatives of FUN evaluated at X0. If
%   FUN has a vector output, DY is a matrix and each row of DY corresponds 
%   to the gradient for a different output.
%
%   DY = FGRAD(FUN,X0,METHOD) specifies alternate methods.
%   The default is forward differences. Use an empty matrix [] to specify
%   the default. Available methods are:
%
%     'forward'     - forward differences       (D+1 points)
%     'backward'    - backward differences      (D+1 points)
%     'central'     - central differences       (2*D points)
%     'five-points' - five-point stencil        (4*D points)
%   
%   D is the number of input dimensions (length of X0). Methods that use 
%   more points to evaluate the finite differences are more precise.
%
%   DY = FGRAD(...,'Vectorized') assumes that FUN can operate on a 
%   vectorized input X. X is a matrix where each row corresponds to a 
%   separate input. The output F of FUN should be also a vector array in
%   which the i-th row contains the value of FUN evaluated on the i-th row 
%   of X.
%
%   DY = FGRAD(...,'Step',H) specifies the step used to compute finite
%   differences. The default is SQRT(EPS).
%
%   [DY,Y0] = FGRAD(...) also returns the value of FUN evaluted at X0.
%
%   See also DIFF.

% Copyright (c) 2015 Luigi Acerbi.

if nargin < 3 || isempty(method); method = 'forward'; end

% Default options
Vectorized = false;
Transpose = false;
h = sqrt(eps)*(abs(x0)+sqrt(eps));

% Parse additional options
idx = 1;
while idx <= numel(varargin)
    switch lower(varargin{idx})
        case 'step'; h = varargin{idx+1}; idx = idx + 2;
        case 'vectorized'; Vectorized = true; idx = idx + 1;
        case 'transpose'; Transpose = true; x0 = x0'; idx = idx + 1;
        otherwise
            error('Unknown options %s.', varargin{idx});
    end
end

if ~Transpose && size(x0,1) > 1
    for i = 1:size(x0,1)
        dy(i,:) = fgrad(fun,x0(i,:),method,varargin{:});
    end
    return;
end

D = size(x0,2);
if isscalar(h); h = h*ones(1,D); end
K = diag(h);

% Evaluate f(x0) only if y0 is requested
EvaluateFx0 = nargout > 1;

% Prepare vectors for finite differences methods
switch lower(method)
    case 'central'
        xi = bsxfun(@plus,x0,0.5*[K;-K]);
    case 'forward'
        xi = bsxfun(@plus,x0,K);
        EvaluateFx0 = 1;
    case 'backward'
        xi = bsxfun(@minus,x0,K);
        EvaluateFx0 = 1;
    case 'five-points'
        xi = bsxfun(@plus,x0,[2*K;K;-K;-2*K]/12);
end

% Add x0 to the evaluation list only if f(x0) is requested
if EvaluateFx0; xi = [xi; x0]; end

ni = size(xi,1);

% Evaluate function at each point in the evaluation list
if Vectorized || ni == 1
    yi = fun(xi);        
else
    if Transpose
        yi = fun(xi(1,:)');        
    else
        yi = fun(xi(1,:));
    end
    nout = size(yi,2);
    yi = [yi; zeros(ni-1,nout)];
    if Transpose
        for i = 2:ni; yi(i,:) = fun(xi(i,:)'); end        
    else
        for i = 2:ni; yi(i,:) = fun(xi(i,:)); end
    end
end

if EvaluateFx0; y0 = yi(end,:); end

% Compute gradient
switch method
    case 'central'
        dy = bsxfun(@rdivide, yi(1:D,:) - yi(D+(1:D),:), h(:));
    case 'forward'
        dy = bsxfun(@rdivide, yi(1:D,:) - y0, h(:));
    case 'backward'
        dy = bsxfun(@rdivide, y0 - yi(1:D,:), h(:));
    case 'five-points'
        dy = bsxfun(@rdivide, -yi(1:D,:) + 8*yi(D+(1:D),:) - 8*yi(2*D+(1:D),:) + yi(3*D+(1:D),:), h(:));
end

if EvaluateFx0; y0 = y0'; end
dy = dy';

end