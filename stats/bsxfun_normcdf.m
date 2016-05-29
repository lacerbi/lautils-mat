function p = bsxfun_normcdf(x,mu,sigma)
%BSXFUN_NORMCDF Vectorized normal cumulative distribution function (cdf).
%   P = BSXFUN_NORMCDF(X,MU,SIGMA) returns the cdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. Dimensions of X, MU, and SIGMA must either match, or 
%   be equal to one. Computation of the pdf is performed with singleton
%   expansion enabled via BSXFUN. The size of Y is the size of the input 
%   arguments (expanded to non-singleton dimensions).
%
%   All elements of SIGMA are assumed to be non-negative (no checks).
%
%   See also BSXFUN, BSXFUN_NORMPDF, NORMCDF.

%   Author: Luigi Acerbi
%   Release date: 15/07/2015

if nargin<3
    error('bmp:bsxfun_normcdf:TooFewInputs','Input argument X, MU or SIGMA are undefined.');
end

try
    if isscalar(mu)
        z = bsxfun(@rdivide, x-mu, sigma);
    elseif isscalar(sigma)
        z = bsxfun(@minus, x, mu)./sigma;
    else
        z = bsxfun(@rdivide, bsxfun(@minus, x, mu), sigma);
    end
    
    % Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))), 
    % to produce accurate near-zero results for large negative x.
    p = 0.5 * erfc(-z ./ sqrt(2));

catch
    error('bmp:bsxfun_normpdf:InputSizeMismatch',...
          'Non-singleton dimensions must match in size.');
end
