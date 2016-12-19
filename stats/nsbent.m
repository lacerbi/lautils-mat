function [H,Hvar] = nsbent(x,flag)
%NSBENT Nemenman-Shafee-Bialek estimator of entropy for discrete distributions.
%   H = nsbent(X) returns the estimated entropy of a discrete array of 
%   counts X. If X is a matrix, NSBENT computes the entropy for each row.
%
%   Reference:
%   Grassberger, P. (2003). Entropy estimates from insufficient samplings. 
%   arXiv preprint physics/0307138.

%   Author: Luigi Acerbi
%   Date: 17/Dec/2016

if nargin < 2 || isempty(flag); flag = 0; end

if ~ismatrix(x); error('X should be an array or two-dimensional matrix.'); end
if size(x,2)==1; x=x'; end      % Single array

K = size(x,2);
Nj = sum(x,2);
h = unique(x(:));   % List of counts    

a_vec(1,1,:) = linspace(0, 10, 2^10);
da = a_vec(2) - a_vec(1);

% Unnormalized prior over concentration parameter a
upriora_vec = K*psi(1, K*a_vec + 1) - psi(1, a_vec + 1);

% Log likelihood
logpxa = bsxfun(@plus, gammaln(Nj+1), ...
    bsxfun(@plus, gammaln(K*a_vec) - K*gammaln(a_vec), ...
    - gammaln(bsxfun(@plus, Nj, K*a_vec)) ...
    + sum(bsxfun(@minus, gammaln(bsxfun(@plus, x, a_vec)), gammaln(x+1)),2)));

% Log unnormalized posterior
%logupost = bsxfun(@plus, logpxa, log(upriora_vec));
logupost = bsxfun(@plus, sum(logpxa,1), log(upriora_vec));

% Normalize posterior
logupost = bsxfun(@minus, logupost, max(logupost,[],3));
post = exp(logupost);
post(isnan(post)) = 0;
post = bsxfun(@rdivide, post, qtrapz(post,3)*da);

% plot(squeeze(a_vec),squeeze(post),'-k','LineWidth',1); hold on; drawnow;

% Entropy estimate
ntilde = bsxfun(@plus, x, a_vec);
Ntilde_vec = sum(ntilde, 2);

% Precompute polygamma values
psi0N1 = psi(Ntilde_vec + 1);
psi0n1 = psi(ntilde + 1);

% Integrand for the mean
Im1 = bsxfun(@times, post, ...
    psi0N1 - sum(bsxfun(@rdivide, ntilde .* psi0n1, Ntilde_vec),2));
H = qtrapz(Im1, 3)*da;

if nargout > 1
    % Precompute polygamma values
    psi0N2 = psi(Ntilde_vec + 2);
    psi1N2 = psi(1,Ntilde_vec + 2);
    psi0n1b = permute(psi0n1,[1 4 3 2]);

    % Compute variance
    ntilde_b = permute(ntilde, [1 4 3 2]);
    Iik = bsxfun(@minus, ...
        bsxfun(@times, bsxfun(@minus, psi0n1, psi0N2), bsxfun(@minus, psi0n1b, psi0N2)), ...
        psi1N2);
    Ji = bsxfun(@minus, bsxfun(@plus, bsxfun(@minus, psi(ntilde+2), psi0N2).^2, psi(1, ntilde + 2)), psi1N2);

    mask(1,:,1,:) = ones(K) - eye(K); % Zero diagonal

    S = bsxfun(@times, mask, bsxfun(@times, ntilde, bsxfun(@rdivide, ntilde_b, (Ntilde_vec+1).*Ntilde_vec)) .* Iik);
    S2 = bsxfun(@times, bsxfun(@rdivide, (ntilde + 1).*ntilde, (Ntilde_vec+1).* Ntilde_vec), Ji);
    Im2 = bsxfun(@times, post, sum(sum(S,2),4) + sum(S2,2));
    Hvar = qtrapz(Im2, 3)*da - H.^2;
end