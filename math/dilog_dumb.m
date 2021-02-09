function y = dilog_dumb(x,nmax)
%DILOG_DUMB Compute dilogarithm (Spence's function) in [0,1] via series.

if nargin < 2 || isempty(nmax); nmax = 1e6; end

y = zeros(size(x));
k = 1:ceil(nmax);
ik2 = 1./k.^2;

for i = 1:numel(x)
    y(i) = sum(x(i).^k.*ik2);
end

y(x < 0 | x > 1) = NaN;
y(x == 1) = pi^2/6;

end
    

