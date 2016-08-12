function [mu,sigma,post] = psychofit(x,y,mu_vec,logsigma_vec)
%PSYCHOFIT Maximum likelihood fit of yes/no psychometric function.
%   [MU,SIGMA] = PSYCHOFIT(X,Y) returns the mean MU and standard deviation
%   SIGMA of a sigmoidal curve (normal cdf) fitted to psychometric data
%   with stimulus value X and response Y. Y is expected to be binary.

MINP = sqrt(eps);

ndef = 51;

[x_vec(1,1,:),~,idx] = unique(x);
nx = numel(x_vec);
ny1 = zeros(1,1,nx);
ny0 = zeros(1,1,nx);

for i = 1:nx
    yidx = y(idx == i);
    ny1(i) = sum(yidx == 1);
    ny0(i) = sum(yidx == 0);    
end

if nargin < 3 || isempty(mu_vec)
    delta = x_vec(end)-x_vec(1);
    mu_vec = [x_vec(1)-0.5*delta,x_vec(end)+0.5*delta];
end
if numel(mu_vec) == 2; mu_vec = linspace(mu_vec(1),mu_vec(2),ndef); end

if nargin < 4 || isempty(logsigma_vec)
    logsigma_vec = log(x_vec(end)-x_vec(1)) + [-log(1e3),log(2)];
end
if numel(logsigma_vec) == 2; logsigma_vec = linspace(logsigma_vec(1),logsigma_vec(2),ndef); end

mu_vec = mu_vec(:);
logsigma_vec = logsigma_vec(:)';

sigma_vec = exp(logsigma_vec);


% Likelihood
logf1 = log1p(erf(bsxfun(@rdivide,bsxfun(@minus,x_vec,mu_vec),sqrt(2)*sigma_vec)));
logf0 = log1p(-erf(bsxfun(@rdivide,bsxfun(@minus,x_vec,mu_vec),sqrt(2)*sigma_vec)));

% Log likelihood
ll = sum(bsxfun(@times, logf1, ny1) + bsxfun(@times, logf0, ny0),3);

% Maximum likelihood
[maxll,idx1] = max(ll,[],1);
[~,idx2] = max(maxll);
mu = mu_vec(idx1(idx2));
sigma = sigma_vec(idx2);

options.Display = 'off';
options.DerivativeCheck = 'off';
options.GradObj = 'on';
% options.TolX = 1e-4;

x0 = [mu,log(sigma)];

theta = gpminimize(x0, @nloglike, 2e3);
% tic; theta = fminunc(@nloglike,x0,options); toc
% [theta; theta2]

mu = theta(1);
sigma = exp(theta(2));

if nargin > 2
    % Return Bayesian posterior (with flat prior)

    post = exp(ll - max(ll(:)));
    % nf = qtrapz(qtrapz(post,1),2);
    % post = post./nf;

    if nargin > 1
        %mu = qtrapz(mu_vec.*qtrapz(post,2),1);
        if nargin > 2
            %sigma = qtrapz(sigma_vec.*qtrapz(post,1),2);
        end
    end
end

function [nll,dtheta] = nloglike(theta)
%NLOGLIKE Function to minimize
    mu_new = theta(1);
    sigma_new = exp(theta(2));
    z = (x_vec-mu_new)/(sqrt(2)*sigma_new);
    erz = erf(z);
    logf1 = log1p(erz);
    logf0 = log1p(-erz);
    
    % Negative log likelihood
    nll = -sum(logf1.*ny1 + logf0.*ny0);

    % Derivatives
    tt = exp(-z.^2) .* (-ny1./(1 + erz) + ny0./(1 - erz));
    dtheta(1) = -sqrt(2./pi)/sigma_new .* sum(tt);
    dtheta(2) = -sqrt(4./pi) .* sum(z .* tt);
end

end

