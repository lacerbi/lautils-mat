function x = bayes1(x1,x2,f,lossfun,n)
%BAYES1 1-d Gaussian mixture (gmm) probability density function (pdf).
%   Y = GAUSSMIXPDF(X,W,MU,SIGMA) returns the pdf of the 1-dimensional
%   Gaussian mixture model (gmm) with mixing weights W, means MU and 
%   standard deviations SIGMA, evaluated at the values in X.
%   The size of Y is the size of X. W, MU and SIGMA can either be vectors 
%   of the same size or scalars. A scalar input parameter functions as a 
%   constant vector of the same size as the other input parameters. 
%   
%   W and SIGMA must be positive-valued, otherwise a vector of NaNs is 
%   returned.
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   [Y,G]= GAUSSMIXPDF(...) returns the first derivative G of the gmm pdf, 
%   evaluated at the values in X.
%
%   [Y,G,H]= GAUSSMIXPDF(...) returns the second derivative H of the gmm 
%   pdf, evaluated at the values in X.
%
%   The gmm (unnormalized) pdf has shape: 
%      gmm(x) = sum(W.*normpdf(x, MU, SIGMA))
%   The pdf is normalized only if the weight vector W sums to one.
%
%   See also GMM1CDF, GMM1MAX, NORMPDF.

%   Copyright (c) by Luigi Acerbi, August 2013

if nargin < 4; lossfun = 'map'; end
if nargin < 5; n = 1000; end

% x = linspace(a, b, n);

if ~ischar(lossfun); error('a'); end
  
switch lower(lossfun)

    case 'map'
        % Brute force max-finding
        xmesh = linspace(x1, x2, n);
        [fx index] = max(f(xmesh));
        dx = diff(xmesh(1:2));        
        [xa, fa] = fminbnd(@(x) (-f(x)), max(xmesh(index)-dx,x1), min(xmesh(index)+dx,x2));
        
        % Precise max-finding (although it may fail!)
        [xb, fb] = fminbnd(@(x) (-f(x)),x1,x2);
        
        % Use the best result
        if fb < fa; x = xb; else x = xa; end
    
    case {'mean', 'mse'}
        xmesh = linspace(x1, x2, n);
        fx = f(xmesh);
        x = qtrapz(fx.*xmesh)/qtrapz(fx);
        % x = quad(@(x) x.*f(x), x1, x2)/quad(f, x1, x2);
            
    case 'median'
        
    otherwise
        
    
end



end