function T = trap(funfcn,a,b,n)
%TRAP   Numerically evaluate integral, trapezoidal rule.
%   T = TRAP(FUN,A,B) approximates the integral of scalar-valued
%   function FUN from A to B using the trapezoidal rule. The interval
%   A-to-B is divided by default in 100 steps. 
%   The function Y=FUN(X) should accept a vector argument X and return a 
%   vector result Y, the integrand evaluated at each element of X.  
%
%   T = TRAP(FUN,A,B,N) divides the interval in N steps. 
%
%   Use array operators .*, ./ and .^ in the definition of FUN
%   so that it can be evaluated with a vector argument.
%
%   Notes:
%   TRAP is fast but it might yield an enormous error.  
%   Might implement later an adaptive method.
%
%   See also QUAD, TRAPZ, FUNCTION_HANDLE.

f = fcnchk(funfcn);
if nargin < 4 || isempty(n), n = 100; end;
if ~isscalar(a) || ~isscalar(b)
   error('MATLAB:trap:scalarLimits',...
         'The limits of integration must be scalars.');     
end

T = 0;
Z = linspace(a, b, n);
dz = (b - a)/(n-1);
if (size(Z, 2) < 2); return; end
T = dz*qtrapz(funfcn(Z));

end