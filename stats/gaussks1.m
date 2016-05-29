function yi = gaussks1(x,y,xi,h)
%GAUSSKS1 1-D Gaussian kernel smooth function approximation
%   YI = GAUSSKS1(X,Y,XI) interpolates to find YI, the values of the
%   underlying function Y at the points in the array XI, smoothed by a
%   Gaussian kernel with unit standard deviation. X and Y must vectors of 
%   length N.
%
%   YI = GAUSSKS1(X,Y,XI,H) uses a bandwidth parameter H (standard
%   deviation of the Gaussian kernel).

if nargin<1
    error('gaussks1:TooFewInputs','Input argument X is undefined.');
end
if nargin<2
    error('gaussks1:TooFewInputs','Input argument Y is undefined.');
end
if nargin<3
    error('gaussks1:TooFewInputs','Input argument XI is undefined.');
end
if nargin<4
    h = 1;
end

yi = zeros(size(xi,1),size(x,2));
for i = 1:length(xi)
    z = exp(-0.5*((xi(i)-x)/h).^2);
    yi(i) = sum(z.*y)/sum(z);
end
    
end