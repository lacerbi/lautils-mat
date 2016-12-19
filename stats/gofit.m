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

H = 0;
Hvar = 0;
ll_chance = 0;
mnf = 0;
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
    end
    % H = H - nansum(Ncounts.*nansum(pk .* log(pk),2),1);
    
    % Chance log likelihood
    q = 1 / K;
    ll_chance = ll_chance + log(q)*sum(Nj,1);
    
    % Multinomial normalization factor correction
    tmp = 0;
    if K == 2
        for j = 1:M
            if Nj(j) > 0
                kk = 0:Nj(j);
                binoln = gammaln(Nj(j)+1) ...
                    - gammaln(kk+1) - gammaln(Nj(j)-kk+1);
                tmp = tmp + nansum(exp(binoln).*binoln.*pk(j,1).^kk.*pk(j,2).^(Nj(j)-kk),2);
            end
        end
    else
        tmp = sum(gammaln(Nj+1));        
        for j = 1:M
            if Nj(j) > 0
                kk = 0:Nj(j);
                bino = exp(gammaln(Nj(j)+1) ...
                    - gammaln(kk+1) - gammaln(Nj(j)-kk+1));
                for k = 1:K
                    tmp = tmp - nansum(bino.*gammaln(kk+1).*pk(j,k).^kk.*(1-pk(j,k)).^(Nj(j)-kk),2);                    
                end
            end
        end
    end
    mnf = mnf + tmp;
end

mnf = 0;

% Goodness of fit
[g,g_ci] = my_goodness(H,Hvar,ll,ll_chance,mnf,Ntot);

if nargout > 1
    % Goodness of fit of histogram model
    ll_hist = 0;
    for i = 1:Nblocks
        qbar = nansum(x{i},1);
        qbar = qbar ./ sum(qbar);
        ll_hist = ll_hist + sum(sum(bsxfun(@times, x{i}, log(qbar)),1),2);
    end
    [gbar,gbar_ci] = my_goodness(H,Hvar,ll_hist,ll_chance,mnf,Ntot);
end

if nargout > 2
    output.H = H;
    output.Hvar = Hvar;
    output.ll = ll + mnf;
    output.ll_chance = ll_chance + mnf;
    output.mnf = mnf;
    output.g_ci = g_ci;
    output.gbar_ci = gbar_ci;
    % output.model_correct = exp((-H -ll -mnf )/Ntot);
    % output.data_correct = exp((-H)/Ntot);
end

%--------------------------------------------------------------------------
function [g,g_ci] = my_goodness(H,Hvar,ll,ll_chance,mnf,Ntrials)
%MY_GOODNESS Compute absolute goodness of fit and credible interval

if Hvar == 0
    g = 1 - (-H -ll -mnf ) ./ (-H -ll_chance -mnf );
    g_ci = [];
else
    Hrange = sqrt(Hvar)*linspace(-5,5,101) + H;
    dH = Hrange(2)-Hrange(1);
    npdf = 1/sqrt(2*pi*Hvar)*exp(-0.5*(Hrange-H).^2/Hvar);
    g = 1 - qtrapz(npdf .* max(-Hrange -ll -mnf,0) ./ max(-Hrange -ll_chance -mnf,eps), 2)*dH;
    
    % 95% credible intervals (fast but wrong when below chance level)
    g_ci(1) = 1 - max(-H + 1.96*sqrt(Hvar) -ll -mnf,0) ./ max(-H  + 1.96*sqrt(Hvar) -ll_chance -mnf, eps);
    g_ci(2) = 1 - max(-H - 1.96*sqrt(Hvar) -ll -mnf,0) ./ max(-H  - 1.96*sqrt(Hvar) -ll_chance -mnf, eps);
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