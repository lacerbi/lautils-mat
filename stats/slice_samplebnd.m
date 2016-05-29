function [samples fval] = slice_samplebnd(N, burn, logdist, xx, widths, lb, ub, step_out, varargin)
%SLICE_SAMPLEBND simple axis-aligned implementation of slice sampling for vectors
% with constrains.
%
%     samples = slice_sample(N, burn, logdist, xx, widths, step_out, varargin)
%
% Inputs:
%             N  1x1  Number of samples to gather
%          burn  1x1  after burning period of this length
%       logdist  @fn  function logprobstar = logdist(xx, varargin{:})
%            xx  Dx1  initial state (or array with D elements)
%        widths  Dx1  or 1x1, step sizes for slice sampling
%            lb  Dx1  lower bounds
%            ub  Dx1  upper bounds
%      step_out bool  set to true if widths may sometimes be far too small
%      varargin   -   any extra arguments are passed on to logdist
%
% Outputs:
%      samples  DxN   samples stored in columns (regardless of original shape)
%      fval     1xN   function value at the sample
%
% Iain Murray May 2004, tweaks June 2009, a diagnostic added Feb 2010
% See Pseudo-code in David MacKay's text book p375

bb = ub - lb;
xx = asin(2*(xx - lb)./bb-1);

% startup stuff
D = numel(xx);
samples = zeros(D, N);
if numel(widths) == 1
    widths = repmat(widths, D, 1);
end

% output function value only if requested
if nargout > 1; fval = zeros(1, N); end

log_Px = feval(logdist, reconvert(xx), varargin{:});

len_iter = floor(N/20);
% Main loop
for ii = 1:(N+burn)
    if mod(ii - burn, len_iter) == 0
        fprintf('Iteration %d                 \r', ii - burn);
    end
    log_uprime = log(rand) + log_Px;

    % Sweep through axes (simplest thing)
    for dd = 1:D
        x_l = xx;
        x_r = xx;
        xprime = xx;

        % Create a horizontal interval (x_l, x_r) enclosing xx
        rr = rand;
        x_l(dd) = xx(dd) - rr*widths(dd);
        x_r(dd) = xx(dd) + (1-rr)*widths(dd);
        if step_out
            % Typo in early editions of book. Book said compare to u, but it should say u'
            while (feval(logdist, reconvert(x_l), varargin{:}) > log_uprime)
                x_l(dd) = x_l(dd) - widths(dd);
            end
            while (feval(logdist, reconvert(x_r), varargin{:}) > log_uprime)
                x_r(dd) = x_r(dd) + widths(dd);
            end
        end

        % Inner loop:
        % Propose xprimes and shrink interval until good one found
        zz = 0;
        while 1
            zz = zz + 1;
            if zz == 10 || zz >= 20
                fprintf('Iteration %d   Step %d       \r', ii - burn, zz);
            end
            xprime(dd) = rand()*(x_r(dd) - x_l(dd)) + x_l(dd);
            log_Px = feval(logdist, reconvert(xprime), varargin{:});
            if log_Px > log_uprime
                break % this is the only way to leave the while loop
            else
                % Shrink in
                if xprime(dd) > xx(dd)
                    x_r(dd) = xprime(dd);
                elseif xprime(dd) < xx(dd)
                    x_l(dd) = xprime(dd);
                else
                    error('BUG DETECTED: Shrunk to current position and still not acceptable.');
                end
            end
        end
        xx(dd) = xprime(dd);
    end

    % Record samples
    if ii > burn
        samples(:, ii - burn) = reconvert(xx(:));
        if nargout > 1; fval(ii - burn) = log_Px; end
    end
end
fprintf('\n');

return;

    function xold = reconvert(xhat)
        xold = 0.5*(sin(xhat)+1).*bb + lb;
    end

end
