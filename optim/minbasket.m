function idx = minbasket(X,y,k)
%MINBASKET Return candidate minima from basket of points.

[n,d] = size(X);

mu = mean(X,1);
X = bsxfun(@minus,X,mu);
M = X'*X / (n-1);
[U,S] = svd(M);
X = X*U'*diag(1./sqrt(diag(S+eps)));

ymu = mean(y);
ystd = std(y);
y = (y - ymu)/ystd * d; % Same weight as all other dimensions combined

idx = fastkmeans([X, y(:)],k,struct('Preprocessing','none'));


end