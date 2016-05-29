function alpha = qrandvm_disc(theta,kappa,n,drop)
%QRANDDVM Quick discretised Von Mises distributed pseudorandom numbers.
%
%   ALPHA = QRANDDVM(THETA,KAPPA,N) returns a column vector of N pseudorandom 
%   angles drawn from the discretized Von Mises distribution with circular 
%   mean THETA and concentration parameter KAPPA. THETA and KAPPA can be 
%   scalars or matrices with distinct parameter values.
%   N can be a scalar or a vector. If N is a vector (e.g., [2 10 8]), the 
%   function creates a matrix output with the respective dimensionality. 
%   If N is not specified, the default value is 1 if both THETA and KAPPA 
%   are scalars, otherwise the size of THETA and/or KAPPA.
%   The returned angles are expressed in radians in the range [-pi, pi].
%   The angles are discretized in increments of 1 deg (pi/180 radians).
%   QRANDDVM uses a lookup table that is computed once and stored in 
%   memory (as a persistent variable) to speed up subsequent function
%   calls. The lookup table occupies about 70 MB of memory.
%
%   ALPHA = QRANDDVM(THETA,KAPPA,N,1) clears the lookup table after the 
%   function call (not recommended).
%
%   For maximal speed, take advantage of the vectorization of QRANDDVM 
%   by generating all required Von Mises random numbers in a single call.
%
%   See also QRANDVM, QTRAPZ.

% Author:   Luigi Acerbi
% Date:     April 24, 2016
% Email:    <luigi.acerbi@gmail.com>

if nargin < 4 || isempty(drop); drop = 0; end

persistent lookuptable;

kappamin = 1e-3;
kappamax = 1e3;

if isempty(lookuptable)
    lookuptable.Nkappa = 5e4;
    lookuptable.Nalpha = 181;
    lookuptable.logkappa = linspace(log(kappamin),log(kappamax),lookuptable.Nkappa)';
    lookuptable.kappa = exp(lookuptable.logkappa);
    lookuptable.alpha = linspace(0,pi,lookuptable.Nalpha);
    
    % Compute Von Mises PDF on a finer grid first
    Nfine = 10;
    alpha_fine = linspace(0,pi,lookuptable.Nalpha*Nfine);
    logpdf_fine = bsxfun(@times, lookuptable.kappa, cos(alpha_fine));
    logpdf_fine = bsxfun(@minus, logpdf_fine, max(logpdf_fine,[], 2));
    pdf_fine = exp(logpdf_fine);
    pdf = zeros(lookuptable.Nkappa, lookuptable.Nalpha);    
    % Use QTRAPZ for quick numerical integration if installed
    if exist('qtrapz','file'); trapzfun = @qtrapz; else trapzfun = @trapz; end
    
    % Compute PDF on coarser grid
    for i = 1:lookuptable.Nalpha
        pdf(:,i) = trapzfun(pdf_fine(:,(1:Nfine) + Nfine*(i-1)),2);
    end
    % Compute discrete CDF and normalise
    lookuptable.cdf = cumsum(pdf,2);
    lookuptable.cdf = bsxfun(@rdivide,lookuptable.cdf,lookuptable.cdf(:,end));
    
    clear Nfine alpha_fine logpdf_fine pdf_fine pdf;
end

% Default parameters
if nargin < 1 || isempty(theta) 
    theta = 0;
end
if nargin < 2 || isempty(kappa)
    kappa = 1;
end

n1 = size(theta);
n2 = size(kappa);

if nargin < 3 || isempty(n)
    if isscalar(theta) && isscalar(kappa)
        n = 1;
    elseif isscalar(theta)
        n = size(kappa);
    else
        n = zeros(1,max(numel(n1),numel(n2)));
        n(1:numel(n1)) = n1;
        n(1:numel(n2)) = max(n(1:numel(n2)),n2);
    end
end

dmax = max([numel(n1),numel(n2),numel(n)]);

n1 = [n1,ones(1,dmax-numel(n1))];
n2 = [n2,ones(1,dmax-numel(n2))];
n = [n,ones(1,dmax-numel(n))];
ntot = prod(n);

if ~isscalar(theta) && ~all(n1 == n | n1 == 1)
    error('THETA must be a scalar or have the same size as KAPPA and the size specified by N.');
end

if ~isscalar(kappa) && ~all(n2 == n | n2 == 1)
    error('KAPPA must be a scalar or have the same size as THETA and the size specified by N.');
end

kappa = repmat(kappa, n./n2);   % Expand KAPPA to all dimensions
kappa = kappa(:);

logkappa = log(kappa);

% Compute position in the lookup table
dlogkappa = lookuptable.logkappa(2)-lookuptable.logkappa(1);
pos = 1 + round( (logkappa - lookuptable.logkappa(1))/dlogkappa );

smallkappa = pos < 1;
largekappa = pos > lookuptable.Nkappa;
medkappa = ~smallkappa & ~largekappa;

alpha = zeros(ntot,1);

% Small kappa: draw from a uniform distribution, then discretise
if any(smallkappa)
    alpha(smallkappa) = round((lookuptable.Nalpha-1)*rand([sum(smallkappa,1),1]));
end

% Large kappa: draw from a discretized Gaussian, then discretise
if any(largekappa)
    sigma = 1./sqrt(kappa(largekappa));
    alpha(largekappa) = round((lookuptable.Nalpha-1)/pi*sigma.*abs(randn([sum(largekappa,1),1])));
end

% Intermediate kappa: draw from lookup table
if any(medkappa)
    posmed = pos(medkappa);
    if 0
        r = mnrnd_private(lookuptable.cdf(posmed,:));
    else
        cdf = lookuptable.cdf(posmed,:);
        m = size(cdf,1);
        samp_k = bsxfun(@gt, cdf, rand(m,1));
        idx = bsxfun(@times, samp_k, 1:size(samp_k,2));
        [~,r] = min(idx,[],2);
    end
    
    alpha(medkappa) = lookuptable.alpha(r);
end

% Randomize direction
flip = rand(size(alpha,1),1) < 0.5;
alpha(flip) = -alpha(flip);

% Reshape ALPHA to the desired matrix shape
alpha = reshape(alpha,n);

alpha = mod(bsxfun(@plus, alpha, theta)+pi,2*pi)-pi;
alpha = round((lookuptable.Nalpha-1)/pi*alpha)/(lookuptable.Nalpha-1)*pi;

end

%--------------------------------------------------------------------------

%MNRND_PRIVATE Random sample from multinomial CDF.
function r = mnrnd_private(cdf)
    m = size(cdf,1);
    samp_k = bsxfun(@gt, cdf, rand(m,1));
    idx = bsxfun(@times, samp_k, 1:size(samp_k,2));
    [~,r] = min(idx,[],2);
    % idx = diff([zeros(m,1), samp_k],[],2);
    % [~,r] = max(idx,[],2);
end
