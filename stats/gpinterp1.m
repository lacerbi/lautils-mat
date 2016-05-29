function yi = gpinterp1(x, y, xi)
%GPINTERP1 1-D interpolation using Gaussian processes
%   YI = INTERP1(X,Y,XI) interpolates to find YI, the values of the
%   underlying function Y at the points in the array XI. X must be a
%   vector of length N.
%   If Y is a vector, then it must also have length N, and YI is the
%   same size as XI.  If Y is an array of size [N,D1,D2,...,Dk], then
%   the interpolation is performed for each D1-by-D2-by-...-Dk value
%   in Y(i,:,:,...,:).
%   If XI is a vector of length M, then YI has size [M,D1,D2,...,Dk].
%   If XI is an array of size [M1,M2,...,Mj], then YI is of size
%   [M1,M2,...,Mj,D1,D2,...,Dk].

% We use Gaussian processes just to perform interpolation
meanfunc = @meanZero;
likfunc = @likGauss;
covfunc = @covSEiso;

% hyp.mean = -10;

% The correlation length is the average step size, whereas the 
% vertical correlation is equal to the average vertical step.
dx = mean(abs(diff(x)));
dy = mean(abs(diff(y)));
hyp.cov = log([dx; dy]);

% We assume a small nonzero noise to favor conditioning
hyp.lik = log(0.001);

% Calculate the Gaussian process interpolation
yi = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x', y', xi');

end
