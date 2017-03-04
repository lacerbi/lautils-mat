function H = grassent(x)
%GRASSENT Grassberger estimator of entropy for discrete distributions.
%   H = grassent(X) returns the estimated entropy of a discrete array of 
%   counts X. If X is a matrix, GRASSENT computes the entropy for each row.
%
%   H = grassent(X,1) computes the estimated entropy with the Poisson
%   approximation (otherwise uses the correct binomial expansion).
%
%   [H,G] = grassent(...) also returns the Grassberger numbers G(h) used to
%   compute the entropy, defined as:
%
%       G(h) = psi(h) + 1/2 * (-1)^h * ( psi((h+1)/2) - psi(h/2) )
%
%   where psi is the digamma function. Note that G is non-NaN only for the
%   values actually used.
%
%   Reference:
%   Grassberger, P. (2003). Entropy estimates from insufficient samplings.
%   arXiv preprint physics/0307138.

%   Author: Luigi Acerbi
%   Date: 15/Dec/2016

if ~ismatrix(x); error('X should be an array or two-dimensional matrix.'); end
if size(x,2)==1; x=x'; end      % Single array

n = sum(x,2);
h = unique([x(:); n(:)]);   % List of counts, including the G(N_k)

% Grassberger numbers shifted by one to include G(0)
G = NaN(1,1+max(h));
G(h+1) = psi(h) + 0.5*(-1).^h .*(psi(0.5*(h+1)) - psi(0.5*h));
G(1) = 0;

% Entropy estimate
base = zeros(size(n));
base(:) = G(n+1);
H = base - sum(x.*G(x+1),2) ./ n;

if nargout > 1
    G = G(2:end);   % Shift back by one
end