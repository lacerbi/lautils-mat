function kappa = fisher2kappa(J,drop)
%FISHER2KAPPA Fisher information to Von Mises concentration parameter.
%
%   KAPPA = FISHER2KAPPA(J) converts the Fisher information value J into
%   Von Mises's concentration parameter KAPPA. J can be a scalar or a
%   matrix of Fisher information values. KAPPA has the same size as J.
%   FISHER2KAPPA uses a lookup table that is computed once and stored in 
%   memory (as a persistent variable) to speed up subsequent function
%   calls. The lookup table occupies about 32 MB of memory.
%
%   KAPPA = FISHER2KAPPA(J,1) clears the lookup table after the function 
%   call (not recommended).
%
%      The relationship between J and kappa is: 
%         J = kappa .* I_1(kappa) ./ I_0(kappa)
%      where I_0 (resp. I_1) is the modified Bessel function of the first 
%      kind of order 0 (resp. 1).
%
%   Reference: 
%   Keshvari, S., van den Berg, R., & Ma, W. J. (2012). Probabilistic 
%   computation in human perception under variability in encoding precision. 
%   PLoS One, 7(6), e40216.
%
%   See also BESSELI.

%   Author: Luigi Acerbi
%   Email: luigi.acerbi@gmail.com
%   Date: Apr 26, 2016

if nargin < 2 || isempty(drop); drop = 0; end

persistent lookuptable;

% Compute lookup table with mapping from J to kappa
if isempty(lookuptable)
    % Change from big and small lookup table
    lookuptable.bigthreshold = 1;
    
    % Lookup table for large values of J (almost linear)
    kapparangeLarge = linspace(lookuptable.bigthreshold,1e4,1e5);
    JrangeLarge = kapparangeLarge.*besseli(1,kapparangeLarge,1)./besseli(0,kapparangeLarge,1);
    lookuptable.JrangeLarge = linspace(lookuptable.bigthreshold,JrangeLarge(end),1e6)';
    lookuptable.kappainvLarge = interp1(JrangeLarge,kapparangeLarge,lookuptable.JrangeLarge,'pchip');

    % Lookup table for small values of J
    kapparangeSmall = linspace(0,lookuptable.bigthreshold+2,1e5);
    JrangeSmall = kapparangeSmall.*besseli(1,kapparangeSmall,1)./besseli(0,kapparangeSmall,1);
    lookuptable.JrangeSmall = linspace(0,lookuptable.bigthreshold,1e6)';
    lookuptable.kappainvSmall = interp1(JrangeSmall,kapparangeSmall,lookuptable.JrangeSmall,'pchip');
end

sz = size(J);
J = J(:);

if any(J < 0)
    error('The Fisher information J needs to be non-negative.');
end

% Compute kappa depending on the case
kappa = zeros(size(J));

% Very large J: kappa(J) is approximately linear (plus an offset)
idx = (J >= lookuptable.JrangeLarge(end));
if any(idx)
    kappa(idx) = lookuptable.kappainvLarge(end) + (J(idx) - lookuptable.JrangeLarge(end));
end
    
% J greater than 1 (high precision)   
idx = (J >= lookuptable.bigthreshold & J < lookuptable.JrangeLarge(end));
if any(idx)
    dJ = lookuptable.JrangeLarge(2)-lookuptable.JrangeLarge(1);
    idx_large = 1 + (J(idx) - lookuptable.JrangeLarge(1))/dJ; 
    idx1 = floor(idx_large);  idx2 = ceil(idx_large);
    w = idx_large - idx1;
    kappa(idx) = (1-w).*lookuptable.kappainvLarge(idx1) + w.*lookuptable.kappainvLarge(idx2);    
end
    
% J between 1e-3 and 1 (low precision)
idx = (J >= 1e-3 & J < lookuptable.bigthreshold);
if any(idx)
    dJ = lookuptable.JrangeSmall(2)-lookuptable.JrangeSmall(1);
    idx_small = 1 + (J(idx) - lookuptable.JrangeSmall(1))/dJ; 
    idx1 = floor(idx_small);  idx2 = ceil(idx_small);
    w = idx_small - idx1;
    kappa(idx) = (1-w).*lookuptable.kappainvSmall(idx1) + w.*lookuptable.kappainvSmall(idx2);        
end

% Very small J: use series expansion for small kappa
idx = J < 1e-3;
if any(idx)
    kappa(idx) = 2*sqrt(1-(sqrt(1 - J(idx))));
end

kappa = reshape(kappa,sz);

% Clear persistent variable
if drop
    lookuptable = [];
end

end