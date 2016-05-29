function nLogL = binolike(p,n,x)
% BINOLIKE Negative log-likelihood for the binomial distribution.
%
%   NLOGL = BINOLIKE(X,N,P) returns the negative of the log-likelihood
%   for the binomial distribution, evaluated at parameters P and N, 
%   given X. 
%
%   See also BINOPDF.

if any(p(:)<0) || any(p(:)>1)
  error('probability should be between 0 and 1');
end

if any(min(x(:))<0)
  error('X should be in the range 0:N')
end

nk  = bsxfun(@minus,bsxfun(@minus,gammaln(n + 1),gammaln(x + 1)),gammaln(bsxfun(@minus,n,x) + 1));
nLogL = -(bsxfun(@plus,nk,bsxfun(@times,x,log(p))) + bsxfun(@times,bsxfun(@minus,n,x),log(1 - p)));