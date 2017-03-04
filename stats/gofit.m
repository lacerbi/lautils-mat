function [g,gbar,output] = gofit(x,ll,method)
%GOFIT Absolute goodness of fit for discrete variables.
%   G = GOFIT(X,LL) returns the absolute goodness of fit for a discrete
%   array of counts X and estimated log likelihood LL. X is a matrix whose
%   elements are data counts; each row is a trial type, and each column a
%   data bin. X can also be a cell array of such data matrices.
%   LL is the total estimated log likelihood of the data (ideally, computed
%   via cross-validation methods, such as leave-one-out).
%
%   G = GOFIT(X,LL,METHOD) specifies the method used to estimate the entropy 
%   of the data. METHOD can be:
%       - 'grassberger' for the Grassberger estimator [1] (default)
%       - 'nsb'         for the Nemenman-Shafee-Bialek (NSB) estimator [2]
%       - 'mle'         for the maximum likelihood or 'plug-in' estimator
%
%   [G,GBAR] = GOFIT(...) also returns the absolute goodness of fit of the
%   histogram model GBAR (with frequencies defined from the empirical 
%   frequencies across batches).
%
%   References:
%   [1] Grassberger, P. (2003). Entropy estimates from insufficient samplings. 
%   arXiv preprint physics/0307138.
%   [2] Nemenman, I., Shafee, F., & Bialek, W. (2001). Entropy and inference, 
%   revisited. arXiv preprint physics/0108025.

%   Author: Luigi Acerbi
%   Date: 16/Dec/2016

if ~iscell(x); x = {x}; end
if nargin < 3 || isempty(method); method = 'grassberger'; end

Nblocks = numel(x);
% Ns = 1e4;
% Hbootstrap = zeros(Ns,1);

H = 0;
Hvar = 0;
ll_chance = 0;
Ntot = 0;
for i = 1:Nblocks

    % Number of counts per condition
    Nj = nansum(x{i},2);
    Ntot = Ntot + sum(Nj);

    % Empirical probabilities
    pk = bsxfun(@rdivide, x{i}, Nj);
    
    % Estimate entropy of the data
    [M,K] = size(x{i});
    
    switch lower(method(1:3))
        case 'nsb'                      % Nemenman-Shafee-Bialek estimator
            [H1,Hvar1] = nsbent(x{i});
            H = H + nansum(Nj.*H1,1);
            Hvar = Hvar + nansum(Nj.*Hvar1,1);
        case 'gra'                      % Grassberger estimator            
            H = H + nansum(Nj.*grassent(x{i}),1);
        case 'mle'                      % Maximum likelihood estimator
            H = H + nansum(Nj.*nansum(-pk .* log(pk),2),1);
    end
        
%     for j = 1:M
%         if any(isnan(pk(j,:))); continue; end
%         D = mnrnd_private(pk(j,:),Nj(j),Ns);            
%         Hbootstrap = Hbootstrap + Nj(j).*grassent(D);         
%     end        
    
    % Chance log likelihood
    ll_chance = ll_chance + log(1/K) * sum(Nj,1);
end

% Goodness of fit
[g,g_ci] = my_goodness(H,Hvar,ll,ll_chance,Ntot);

if nargout > 1
    % Goodness of fit of histogram model
    ll_hist = 0;
    for i = 1:Nblocks
        qbar = nansum(x{i},1);
        qbar = qbar ./ sum(qbar);
        ll_hist = ll_hist + sum(sum(bsxfun(@times, x{i}, log(qbar)),1),2);
    end
    [gbar,gbar_ci] = my_goodness(H,Hvar,ll_hist,ll_chance,Ntot);
end

if nargout > 2
    output.H = H;
    output.Hvar = Hvar;
    output.ll = ll;
    output.ll_chance = ll_chance;
    output.g_ci = g_ci;
    output.gbar_ci = gbar_ci;
end

%--------------------------------------------------------------------------
function [g,g_ci] = my_goodness(H,Hvar,ll,ll_chance,Ntrials)
%MY_GOODNESS Compute absolute goodness of fit and credible interval

if Hvar == 0
    g = 1 - (-H -ll ) ./ (-H -ll_chance );
    g_ci = [];
else
    Hrange = sqrt(Hvar)*linspace(-5,5,101) + H;
    dH = Hrange(2)-Hrange(1);
    npdf = 1/sqrt(2*pi*Hvar)*exp(-0.5*(Hrange-H).^2/Hvar);
    g = 1 - qtrapz(npdf .* max(-Hrange -ll,0) ./ max(-Hrange -ll_chance,eps), 2)*dH;
    % g = 1 - (-Hrange -ll) ./ (-Hrange -ll_chance);
    
    % 95% credible intervals (fast but wrong when below chance level)
    g_ci(1) = 1 - max(-H + 1.96*sqrt(Hvar) -ll,0) ./ max(-H  + 1.96*sqrt(Hvar) -ll_chance, eps);
    g_ci(2) = 1 - max(-H - 1.96*sqrt(Hvar) -ll,0) ./ max(-H  - 1.96*sqrt(Hvar) -ll_chance, eps);
end

%--------------------------------------------------------------------------
function y = nansum(x,dim)
%NANSUM Sum, ignoring NaNs.

x(isnan(x)) = 0;
if nargin == 1
    y = sum(x);
else
    y = sum(x,dim);
end

%--------------------------------------------------------------------------
function r = mnrnd_private(p,n,m)
%MNRND_PRIVATE Random sample from the multinomial distribution.

if nargin < 2 || isempty(n); n = 1; end
if nargin < 3 || isempty(m); m = size(p,1); end
cdf = cumsum(p,2);
rr = rand(m,1,n);
samp_k = bsxfun(@gt, cdf, rr) & bsxfun(@le, [zeros(size(cdf,1),1),cdf(:,1:end-1)], rr);
r = sum(samp_k,3);