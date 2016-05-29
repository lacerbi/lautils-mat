function y = normlogpdf(x,mu,sigma)
%NORMLOGPDF Normal log probability density function (log pdf).
%   Y = NORMLOGPDF(X,MU,SIGMA) returns the log pdf of the normal 
%   distribution with mean MU and standard deviation SIGMA, evaluated at 
%   the values in X. The size of Y is the common size of the input 
%   arguments.  A scalar input functions as a constant matrix of the same 
%   size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   See also NORMPDF.

if nargin<1
    error('normlogpdf:TooFewInputs','Input argument X is undefined.');
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    y = -0.5 * ((x - mu)./sigma).^2 - log(sigma) - 0.5*log(2*pi);
catch
    error('normlogpdf:InputSizeMismatch',...
          'Non-scalar arguments must match in size.');
end
