function [H,z0] = fhess(fun,x0,method,varargin)
%FHESS Estimation of function Hessian via finite differences
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

if nargin < 3 || isempty(method); method = 'central'; end

% Default options
Vectorized = false;
h = eps^(1/4);

% Parse additional options
idx = 1;
while idx <= numel(varargin)
    switch lower(varargin{idx})
        case 'step'; h = varargin{idx+1}; idx = idx + 2;
        case 'vectorized'; Vectorized = true; idx = idx + 1;
        otherwise
            error('Unknown options %s.', varargin{idx});
    end
end

D = size(x0,2);

% Prepare vectors for finite differences methods
switch lower(method)
    case 'central'
        xi = [x0; bsxfun(@plus,x0,0.5*h*[eye(D);-eye(D)])];
    case 'five-points'
        xi = [x0; bsxfun(@plus,x0,h/12*[2*eye(D);eye(D);-eye(D);-2*eye(D)])];
end
offset = size(xi,1);

lowtrmat = triu(ones(D),1);
lowtrmat(logical(lowtrmat)) = 0:D*(D-1)/2-1; 
lowtrmat = lowtrmat';

ximix = zeros(D*(D-1), D);
for i = 2:D
    for j = 1:i-1
        ximix(lowtrmat(i,j)*2+(1:2),[i j]) = 0.5*h*[1 1; -1 -1];
    end
end
ximix = bsxfun(@plus,x0,ximix);

xi = [xi; ximix];
ni = size(xi,1);

% Evaluate function at each point in the evaluation list
if Vectorized || ni == 1
    zi = fun(xi);        
else
    zi = fun(xi(1,:));
    nout = size(zi,2);
    zi = [zi; zeros(ni-1,nout)];
    for i = 2:ni; zi(i,:) = fun(xi(i,:)); end        
end

z0 = zi(1,:);

% Compute Hessian
H = zeros(D,D);
for i = 1:D
    % Second derivative
    switch method
        case 'central'
            baseidx1 = i+1;
            H(i, i) = (zi(baseidx1) -2*z0 + zi(baseidx1+D))/(0.25*h^2);
    end
    
    % Mixed derivatives
    for j = 1:i-1
        baseidx2 = j+1;
        basemix = offset+lowtrmat(i,j)*2+1;
        H(i, j) = -1/(0.5*h^2)*(zi(baseidx1)+zi(baseidx1+D)+ ...
            zi(baseidx2)+zi(baseidx2+D)-2*z0-zi(basemix)-zi(basemix+1));
    end
    
end

% Symmetrize Hessian
H = H + tril(H,-1)';

if 0
    % Compute gradient
    switch method
        case 'central'
            dy = (yi(1:D,:) - yi(D+(1:D),:))./h;
        case 'forward'
            dy = (yi(1:D,:) - y0)./h;
        case 'backward'
            dy = -(yi(1:D,:) - y0)./h;
        case 'five-points'
            dy = (-yi(1:D,:) + 8*yi(D+(1:D),:) - 8*yi(2*D+(1:D),:) + yi(3*D+(1:D),:))./h;
    end

    y0 = y0';
    dy = dy';
end

end