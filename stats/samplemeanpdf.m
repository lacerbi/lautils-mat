\function z = samplemeanpdf(varargin)
%SAMPLEMEANPDF  Sampling distribution of the mean.
%   Z = QTRAPZ(Y) computes an approximation of the integral of Y via
%   the trapezoidal method (with unit spacing).  To compute the integral
%   for spacing different from one, multiply Z by the spacing increment.
%
%   For vectors, QTRAPZ(Y) is the integral of Y. For matrices, QTRAPZ(Y)
%   is a row vector with the integral over each column. For N-D
%   arrays, QTRAPZ(Y) works across the first non-singleton dimension.
%
%   Z = QTRAPZ(Y,DIM) integrates across dimension DIM of Y. The length of X 
%   must be the same as size(Y,DIM).
%
%   QTRAPZ is potentially much faster than TRAPZ with *large* matrices.
%
%   See also TRAPZ.

% Set up the defaults
narginchk(2,4);
[x,dx,y,n,nclt,dim,TolZero] = parseinputs(varargin{:});
 
% Size of input distribution
if isvector(y); nDims = 1; y = y(:); else nDims = ndims(y); end

% Verify if N is an integer (robust to numerical error)
integern = isnumequal(n,round(n),TolZero);
if integern; n = round(n); end

% Verify that NCLT is Inf or integer
if nclt < Inf
    assert(isnumequal(nclt,round(nclt),TolZero), ...
        'NCLT needs to be Inf or integer-valued.');
    nclt = round(nclt);
end

ngrid = size(y,dim);    % number of grid points

% Compute distribution means
m = bsxfun(@rdivide,qtrapz(bsxfun(@times,y,x),dim),qtrapz(y,dim));

% Number of samples exceed limit for central limit theorem approximation
if n >= nclt
    % Standard deviation of normal approximation is SD(Y)/sqrt(N)
    sd = sqrt(bsxfun(@rdivide,qtrapz(bsxfun(@times,y,x.^2),dim),qtrapz(y,dim)) - m.^2);
    z = bsxfun_normpdf(x,m,sd/sqrt(n))/dx;
    return;

% Return uniform distribution for zero samples
elseif isnumequal(n,0,TolZero)
    z = ones(size(y))/ngrid/dx;
    return;
    
% One sample, return base distribution
elseif isnumequal(n,1,TolZero)
    z = bsxfun(@rdivide,y,qtrapz(y,1))/dx;
    return;
    
% Interpolate between uniform and base distribution for less than 1 sample
elseif n < 1
    z = pdfinterp(ones(size(y)), bsxfun(@rdivide,y,qtrapz(y,dim)), ceil(n)-n);
    z = bsxfun(@rdivide,z,qtrapz(z,1))/dx;
    return;
end

% Compute padding
ntot = 2^(ceil(log2(ceil(n)*ngrid))+1);
pad = ntot - ngrid;

% Shift starting distributions to their mean
switch nDims
    case 1
        index = find(x > m, 1);
        yshifted = [y(index:end); zeros(pad,1); y(1:index-1)];
    case 2
        [~,index] = min(abs(bsxfun(@minus,x,m)),[],1);
        yshifted = zeros(ntot,size(y,2));
        for i = 1:size(y,2)
            yshifted(:,i) = [y(index(i):end,i); zeros(pad,1); y(1:index(i)-1,i)];
        end
    otherwise
        error('Y can have dimension 1 or 2.');
end

if integern
    % Integer number of samples
    z = computesamplemean(n,x,yshifted,dim,m);
else
    % Non-integer number of samples, interpolate b/w adjacent integers
    z1 = computesamplemean(floor(n),x,yshifted,dim,m);
    
    % Check whether central limit theorem approximation applies
    if ceil(n) >= nclt
        sd = sqrt(bsxfun(@rdivide,qtrapz(bsxfun(@times,y,x.^2),dim),qtrapz(y,dim)) - m.^2);
        z2 = bsxfun_normpdf(x,m,sd/sqrt(ceil(n)));
    else
        z2 = computesamplemean(ceil(n),x,yshifted,dim,m);
    end
    z = pdfinterp(z1,z2,ceil(n)-n);
end

% Normalize pdf
z = bsxfun(@rdivide,z,qtrapz(z,1))/dx;

end % SAMPLEMEANPDF

%-------------------------------------------------------------------------%
%% Private functions

%PARSEINPUTS Returns default arguments
function [x,dx,y,n,nclt,dim,TolZero] = parseinputs(varargin)
x = []; dx = 1; nclt = Inf; dim = 1; TolZero = 1e-12;
switch nargin
    case 2
        y = varargin{1}; n = varargin{2};
    case 3
        if isscalar(varargin{2})
            y = varargin{1}; n = varargin{2}; nclt = varargin{3};
        else
            x = varargin{1}(:); y = varargin{2}; n = varargin{3};
        end        
    case 4
            x = varargin{1}(:); y = varargin{2}; n = varargin{3}; ...
                nclt = varargin{4};                
end

if ~isempty(x)
    if isscalar(x)
        dx = x;
    else
        dx = x(2)-x(1);
        if any(abs(diff(diff(x))) > TolZero)
            error('X needs to be an equally spaced grid.');
        end
    end
end
x = (1:size(y,dim))';
end
        
%-------------------------------------------------------------------------%
%COMPUTESAMPLEMEAN Compute sample mean for provided pdf.
function z = computesamplemean(n,x,yshifted,dim,m)
    ntot = size(yshifted,1);
    ngrid = size(x,1);

    % Base interpolation range
    xq = (-ntot/2+1:ntot/2)';
    
    half1_i = (1:ntot/2)';
    half2_i = (ntot/2+1:ntot)';

    yshifted = bsxfun(@rdivide,yshifted,max(yshifted,[],dim));
    
    zz = fft(yshifted,[],dim);
    zz = intpower(zz,n); % INTPOWER is faster than MATLAB's in-built power
    zz = real(ifft(zz,[],dim));
    zz(isnan(zz)) = 0;    
    
    % Expanded interpolation range for z
    zq = n*((1:ntot)-ntot/2)';
    
    z = interp1(xq,[zz(half2_i,:); zz(half1_i,:)],zq);          
    z = z(ntot/2-ngrid/2+1:ntot/2+ngrid/2,:);
    z(isnan(z)) = 0;
    
    m2 = bsxfun(@rdivide,qtrapz(bsxfun(@times,z,x),dim),qtrapz(z,dim));
    xnew = bsxfun(@minus,x,m2-m);
    z = lininterp1(xnew, z, x);
    z(isnan(z)) = 0;
    z = max(z,0);
end

%-------------------------------------------------------------------------%
%PDFINTERP Continuous interpolation of pdfs via geometric mean.
function z = pdfinterp(y1,y2,w)
    z = real(exp(w*log(y1) +(1-w)*log(y2)));
    z = max(z,0);
end

%-------------------------------------------------------------------------%
%NUMEQUAL True if A and B are within numerical zero.
function tf = isnumequal(a,b,TolZero)
    tf = abs(a - b) < TolZero;
end

%-------------------------------------------------------------------------%
%LININTERP1 1-D data interpolation for arrays.
function Vout = lininterp1(X,V,Xq)

if isvector(X); X = X(:); end
if isvector(V); V = V(:); end

% Gets the x spacing
idx = 1./(X(2,:)-X(1,:));           % one over to perform divide only once
Xq = bsxfun(@minus,Xq,X(1,:));      % subtract minimum of x

Xqi = floor(bsxfun(@times,Xq,idx))+1;           % indices of nearest-lower-neighbors
flag = Xqi<1 | Xqi>size(X,1)-1 | isnan(Xqi);    % finds indices out of bounds
V = [V; NaN(1, size(V,2))];
Xqi(flag) = size(V,1)-1;

delta = Xqi - bsxfun(@times,Xq,idx);

linind1 = bsxfun(@plus, Xqi, size(V,1)*(0:size(V,2)-1));
linind2 = bsxfun(@plus, Xqi + 1, size(V,1)*(0:size(V,2)-1));

Vout = bsxfun(@times,delta,V(linind1)) + bsxfun(@times,1-delta,V(linind2));
Vout(flag) = NaN;

end