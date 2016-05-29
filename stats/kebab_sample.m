function [samples loglikes] = kebab_sample(N, burn, logdist, xx, widths, step_out, verbose, varargin)
%KEBAB_SAMPLE slice sampling with random rotation every sweep
%
%     [samples loglikes] = slice_sample(N, burn, logdist, xx, widths, step_out, verbose, varargin)
%
% Inputs:
%             N  1x1  Number of samples to gather
%          burn  1x1  after burning period of this length
%       logdist  @fn  function logprobstar = logdist(xx, varargin{:})
%            xx  Dx1  initial state (or array with D elements)
%        widths  Dx1  or 1x1, step sizes for slice sampling
%      step_out bool  set to true if widths may sometimes be far too small
%                     (default 0)
%       verbose bool  set to true prints log to screen (default 1)
%      varargin   -   any extra arguments are passed on to logdist
%
% Outputs:
%      samples  DxN   samples stored in columns (regardless of original shape)
%      loglikes 1xN   log-likelihood of samples (optional)
%
% Iain Murray May 2004, tweaks June 2009, a diagnostic added Feb 2010
% Luigi Acerbi inteface tweaks Jan 2012, loglikes output Feb 2013
% See Pseudo-code in David MacKay's text book p375

% By default do not step out
if ~exist('step_out', 'var'); step_out = []; end
if isempty(step_out); step_out = 0; end

% By default be verbose
if ~exist('verbose', 'var'); verbose = []; end
if isempty(verbose); verbose = 1; end

% startup stuff
xx = xx(:);
D = numel(xx);
samples = zeros(D, N);
if nargout > 1; loglikes = zeros(1, N); end

if numel(widths) == 1
    widths = repmat(widths, D, 1);
end
widths = widths(:);

log_Px = feval(logdist, xx, varargin{:});

len_iter = floor(N/20);
% Main loop
for ii = 1:(N+burn)
    if verbose && mod(ii - burn, len_iter) == 0
        fprintf('Iteration %d                 \r', ii - burn);
    end
    log_uprime = log(rand) + log_Px;
    
    % Choose a random direction
    v1 = randn(1, D);
    v1 = v1/sqrt(v1*v1');
    base=[v1;null(v1)'];

    % Sweep through axes (simplest thing)
    for dd = 1:D
        xprime = xx;

        % Create a horizontal interval (x_l, x_r) enclosing xx
        rr = rand;
        wvec = base(dd, :)'.*widths;
        wsize = 1;
        
        x_l = xx - rr*wvec;
        x_r = xx + (1-rr)*wvec;
        if step_out
            % Typo in early editions of book. Book said compare to u, but it should say u'
            while (feval(logdist, x_l, varargin{:}) > log_uprime)
                x_l = x_l - wvec;
                wsize = wsize + 1;
                rr = rr + 1;
            end
            while (feval(logdist, x_r, varargin{:}) > log_uprime)
                x_r = x_r + wvec;
                wsize = wsize + 1;
            end
        end

        % Inner loop:
        % Propose xprimes and shrink interval until good one found
        zz = 0;
        while 1
            zz = zz + 1;
            if verbose && (zz == 10 || zz >= 20)
                fprintf('Iteration %d   Step %d       \r', ii - burn, zz);
            end
            
            rr2 = rand()*wsize;
            xprime = x_l + rr2*wvec;
            log_Px = feval(logdist, xprime, varargin{:});
            if log_Px > log_uprime
                break % this is the only way to leave the while loop
            else
                % Shrink in
                if rr2 > rr
                    wsize = rr2;                    
                elseif rr2 < rr
                    x_l = xprime;
                    rr = rr - rr2;
                    wsize = wsize - rr2;
                else
                    errorstr = 'BUG DETECTED: Shrunk to current position and still not acceptable. Current position: ';
                    for k = 1:(length(xx)-1); errorstr = [errorstr num2str(xx(k)) ', ']; end
                    error([errorstr num2str(xx(length(xx)))]);        
                end
            end
        end
        xx = xprime;
    end

    % Record samples
    if ii > burn
        samples(:, ii - burn) = xx(:);
        if nargout > 1; loglikes(ii - burn) = log_Px; end
    end
end
if verbose; fprintf('\n'); end

