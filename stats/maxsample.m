function [x,y] = maxsample(nsamples,fun,x0,width,alpha,y0)
%MAXSAMPLE Stochastic exploration of high-valued regions of a function.
%   [X,Y] = MAXSAMPLE(NSAMPLES,FUN,X0) draws NSAMPLES random samples that
%   explore the domain of function FUN looking for maxima and generally
%   high-valued regions. X0 is a row vector or scalar containing the initial
%   value of the random sample sequences. X0 must be within the domain of
%   the target function. NSAMPLES is the number of samples to be generated.
%   FUN is a function handle created using @. FUN takes only one argument
%   as an input and this argument has the same type and size as X0.
%   X is a matrix whose rows are random samples, Y is a column vector whose
%   elements are the values of the target function for each sample.
%   MINSAMPLE tries to avoid resampling from the already explored regions.
%  
%   [X,Y] = MAXSAMPLE(NSAMPLES,FUN,X0,WIDTH) draws random samples for the
%   target function with a typical width WIDTH. W is a scalar or vector. 
%   If it is a scalar, all dimensions are assumed to have the same typical 
%   widths. If it is a vector, each element of the vector is the typical 
%   width scale of the target function in that dimension. The default value
%   of WIDTH is 0.1.
%
%   [X,Y] = MAXSAMPLE(NSAMPLES,FUN,X0,WIDTH,ALPHA) uses weight ALPHA for 
%   the penalty of visiting already explored regions. ALPHA should be of 
%   the magnitude of the variation of FUN along the typical scales WIDTH.
%   The default value of ALPHA is 1.
%
%   [X,Y] = MAXSAMPLE(NSAMPLES,FUN,X0,WIDTH,ALPHA,Y0) continue a previous
%   sampling run with points X0 and function values Y0.

% Copyright Luigi Acerbi 2014

if nargin < 3; help maxsample; end
if nargin < 4; width = 0.1; end
if nargin < 5; alpha = 1; end
if nargin < 6; y0 = []; end

offset = size(y0, 1);
D = size(x0, 2);
x = zeros(offset + nsamples, D);
y = zeros(offset + nsamples, 1);
ww = zeros(offset + nsamples, 1);
rr = rand(nsamples, 1);
nn = randn(nsamples, D);

x0 = x0(:)';
width = width(:)';
if isscalar(width); width = ones(1, D)*width; end

if offset == 0
    x(1, :) = x0;
    y(1) = fun(x0);
    ww(1) = y(1);
    offsetone = 1;
else
    x(1:offset, :) = x0;
    y(1:offset) = y0;
    ww(1) = y(1);
    % Compute weights
    for i = 2:offset; updateweights(i); end
    offsetone = 0;
end

for i = (offset+offsetone+1):(offset+nsamples) 
    % Choose a candidate
    WW = cumsum(exp(ww(1:i-1)));
    ind = find(rr(i-offset)*WW(i-1) < WW, 1);
    ww(ind) = ww(ind) - alpha;

    % Generate new point till you get a good one
    normr = nn(i-offset, :);
    while 1
        x(i, :) = x(ind, :) + normr.*width;
        y(i) = fun(x(i, :));
        if ~isnan(y) & ~isinf(y); break; end
        normr = randn(1, D);
    end
    
    % x(i, :) - x0
    % y(i)
    
    % Update weights
    updateweights(i);
end

return;

    % Update weights for each sample
    function updateweights(n)
        normdist = sum(abs(bsxfun(@rdivide, bsxfun(@minus, x(n, :), x(1:(n-1), :)), width)), 2)/D;        
        dw = alpha*exp(-normdist);
        % dw = alpha*exp(-0.5/D*sum(bsxfun(@rdivide, bsxfun(@minus, x(i, :), x(1:(i-1), :)), 0.25*width).^2, 2));
        ww(n) = y(n) - sum(dw);
        ww(1:n-1) = ww(1:n-1) - dw;
    end
end