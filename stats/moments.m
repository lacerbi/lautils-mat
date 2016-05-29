function mms = moments(x, p)
% MOMENTS Central moments of a given distribution.
%
%   MMS = MOMENTS(X, P)
%   returns the mean, standard deviation, skewness and kurtosis of a 
%   discrete probability distribution defined over X with probability P. 
%   X and P must be matrices with the same number of columns. 
%   P can have multiple rows, each one counting as a different distribution 
%  (in which case MMS will be a matrix).

% Check that X and P have the same length
if size(x, 2) ~= size(p, 2)
    error('X and P must have the same number of columns.');
end

% Bring both matrices to the same number of rows
if size(x, 1) == 1; x = ones(size(p, 1), 1)*x; end
if size(p, 1) == 1; p = ones(size(x, 1), 1)*p; end

% Rescale each row -- they must be probability distributions
for i = 1:size(x, 1); p(i, :) = p(i, :)/sum(p(i, :)); end

% Calculate the central moments for each row
mms = zeros(size(x, 1), 4);
for i = 1:size(x, 1)
    mms(i, 1) = sum(p(i, :).*x(i, :));
    mms(i, 2) = sqrt(sum(p(i, :).*(x(i, :).^2)) - mms(i, 1)^2);
    mms(i, 3) = (sum(p(i, :).*(x(i, :).^3)) - 3*mms(i,1)*mms(i,2)^2 - mms(i,1)^3)/mms(i,2)^3;
    mms(i, 4) = (sum(p(i, :).*(x(i, :) - mms(i,1)).^4))/mms(i,2)^4 - 3;
end

end