function [samples,fvals,exitflag,output] = temperslice(logf,x0,N,betas,widths,LB,UB,options,varargin)
%TEMPERSLICE Parallel tempering with local slice sampling moves.
%
%   SAMPLES = TEMPERSLICE(LOGF,X0,N,BETAS) draws N random samples from a 
%   target distribution with the log probability density function LOGF 
%   by combining parallel tempering with the coordinate-wise slice sampling 
%   method for local sampling. X0 is a row vector containing the initial 
%   value of the random sample sequences. X0 must be within the domain of 
%   the target distribution. X0 can also be a matrix of M row vectors, where
%   each row corresponds to the starting point for a separate chain (see 
%   below). N is the number of samples to be returned.
%   LOGF is a function handle created using @. LOGF takes one argument as
%   input that has the same type and number of columns as X0 and returns 
%   the target log density function (minus a constant; that is, the 
%   normalization constant of the pdf needs not be known). BETAS is a vector 
%   of M inverse temperatures used by each parallel sampled chain. Usually 
%   0 <= BETAS(m) <= 1, BETAS(m) > BETAS(m+1) and BETAS(1) = 1 (which 
%   corresponds to the original target density). If BETAS is a scalar and 
%   set to 1, TEMPERSLICE behaves as a simple slice sampling algorithm.
%   SAMPLES is a matrix each row of which corresponds to a sampled point in 
%   the random sequence. By default, TEMPERSLICE only returns the parallel 
%   chain for which BETAS(m) == 1 (if such chain does not exist, it returns 
%   the first chain).
%
%   SAMPLES = TEMPERSLICE(LOGF,X0,N,BETAS,WIDTHS) uses WIDTHS as a scalar or 
%   vector of typical widths for the slice sampling algorithm. If WIDTHS is 
%   a scalar, all dimensions are assumed to have the same typical widths. 
%   If it is a row vector, each element of the vector is the typical width 
%   of the marginal target distribution in that dimension. The default value 
%   of W(i) is (UB(i)-LB(i))/2 if the i-th bounds are finite, or 10 otherwise.
%   If WIDTHS is a scalar or row vector, the WIDTHS for the m-th chain are
%   multiplied by 1/sqrt(BETAS(m)) to account for the relaxation of the
%   density. WIDTHS can also be specified as a matrix with M rows, where
%   each row is the vector of WIDTHS for the m-th chain.
%   By default TEMPERSLICE uses an adaptive widths method during the burn-in 
%   period, so the choice of typical widths is not crucial.
%
%   SAMPLES = TEMPERSLICE(LOGF,X0,N,WIDTHS,LB,UB) defines a set of lower 
%   and upper bounds on the domain of the target density function, which is
%   assumed to be zero outside the range LB <= X <= UB. Use empty matrices 
%   for LB and UB if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded 
%   below; set UB(i) = Inf if X(i) is unbounded above. If LB(i) == UB(i),
%   the variable is assumed to be fixed on that dimension.
%
%   SAMPLES = TEMPERSLICE(LOGF,X0,N,WIDTHS,LB,UB,OPTIONS) samples with 
%   the default sampling parameters replaced by values in the structure 
%   OPTIONS. TEMPERSLICE uses these options:
%
%     OPTIONS.Thin generates random samples with Thin-1 out of Thin values 
%     omitted in the generated sequence (after burn-in). Thin is a positive
%     integer. The default value of Thin is 1.
%
%     OPTIONS.Burnin omits the first Burnin points before starting recording 
%     points for the generated sequence. Burnin is a non-negative integer. 
%     The default value of Burnin is round(0.1*N) (that is, 10% of the 
%     number of recorded samples).
%    
%     OPTIONS.StepOut if set to true, performs the stepping-out action when 
%     the current window does not bracket the probability density; see 
%     Neal (2003) for details. StepOut is a boolean value. The default 
%     value of StepOut is false. 
%
%     OPTIONS.Display defines the level of display. Accepted values for
%     Display are 'iter', 'notify', 'final', and 'off' for no display. The 
%     default value of Display is 'off'. 
%    
%     OPTIONS.LogPrior allows the user to specify a prior over X.
%     OPTIONS.LogPrior is a function handle created using @. LogPrior takes 
%     one argument as input that has the same type and size as X0 and 
%     returns the log prior density function at X. The generated samples
%     will be then drawn from the log density LOGF + OPTIONS.LogPrior.
%
%   SAMPLES = TEMPERSLICE(...,VARARGIN) passes additional arguments
%   VARARGIN to LOGF.
%
%   [SAMPLES,FVALS] = TEMPERSLICE(...) returns the sequence of values 
%   FVALS of the target log pdf at the sampled points. If a prior is
%   specified in OPTIONS.LogPrior, FVALS does NOT include the contribution
%   of the prior.
%
%   [SAMPLES,FVALS,EXITFLAG] = TEMPERSLICE(...) returns an EXITFLAG that
%   describes the exit condition of TEMPERSLICE. Possible values of 
%   EXITFLAG and the corresponding exit conditions are
%
%    0  Target number of recorded samples reached, with no explicit
%       violation of convergence (this does not ensure convergence).
%    1  Detected probable lack of convergence of the sampling procedure.
%    2  Detected lack of convergence of the sampling procedure.
%    3  No explicit violation of convergence detected, but the number of
%       effective (independent) samples in the sampled sequence is much 
%       lower than the number of requested samples N for at least one 
%       dimension.
%    -1 Target number of recorded samples reached, convergence status is
%       unknown (no diagnostics have been run).
%
%   [SAMPLES,FVALS,EXITFLAG,OUTPUT] = TEMPERSLICE(...) returns a structure
%   OUTPUT with the number of evaluations of LOGF in OUTPUT.funcCount, the 
%   value WIDTHS used during sampling in OUTPUT.widths (they can differ from 
%   the initial WIDTHS due to width adaptation during the burn-in stage), 
%   and the sequence of the values of the log prior at the sampled points 
%   in OUTPUT.logpriors (only if OPTIONS.LogPrior is nonempty).
%
%   LOGF can return either a scalar (the value of the log probability 
%   density at X) or a row vector (the value of the log probability density
%   at X for each data point; each column corresponds to a different data
%   point). In the latter case, the total log pdf is obtained by summing
%   the log pdf per each individual data point. If LOGPDF returns a vector,
%   FVALS returned by TEMPERSLICE is a matrix (each row corresponds to a 
%   sampled point, each column to a different data point).
%   Knowing the log pdf of the sampled points per each data point can be 
%   useful to compute estimates of predictive error such as the widely 
%   applicable information criterion (WAIC); see Watanabe (2010).

% Author: Luigi Acerbi
% Email: luigi.acerbi@nyu.edu
% Date: feb/18/2016
%
%   References: 
%   - R. Neal (2003), Slice Sampling, Annals of Statistics, 31(3), p705-67.
%   - D. J. MacKay (2003), Information theory, inference and learning 
%     algorithms, Cambridge university press, p374-7.
%   - S. Watanabe (2010), Asymptotic equivalence of Bayes cross validation
%     and widely applicable information criterion in singular learning 
%     theory, The Journal of Machine Learning Research, 11, p3571-94.

%% Default parameters and options

% By default unbounded sampling
if nargin < 5; widths = []; end
if nargin < 6 || isempty(LB); LB = -Inf; end
if nargin < 7 || isempty(UB); UB = Inf; end
if nargin < 8; options = []; end

% Default options
defopts.Swap = 1;               % Swap frequency
defopts.Thin = 1;               % Thinning frequency
defopts.Burnin = round(0.1*N);  % Burnin
defopts.StepOut = 0;            % Step-out
defopts.Display = 'off';        % Display
defopts.AdaptiveWidths = 1;     % Adaptive widths
defopts.LogPrior = [];          % Log prior over X
defopts.Diagnostics = 1;        % Perform convergence diagnostics at the end
defopts.Debug = 1;              % DEBUG mode. Return all chains.
defopts.Tunnel = [];            % Tunnel samples

storeallchains = 1;             % Store all chains while running

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

%% Startup and initial checks
D = size(x0,2);
assert(isvector(betas), ...
    'BETAS needs to be a vector of inverse temperatures for the parallel tempering.');
betas = betas(:);
M = size(betas,1);
targetidx = find(betas == 1,1);
if isempty(targetidx)
    targetidx = 1;
    warning('No inverse temperature in vector BETAS is equal to 1. The sampler will return the sequence corresponding to the first element of BETAS.');
end
dbetas = diff(betas);
if any(dbetas >= 0)
    warning('The inverse temperatures in BETAS are expected to be monotonically decreasing, that is BETA(m+1) < BETA(m) for 1 <= m < M.');
end

if numel(LB) == 1; LB = repmat(LB, [1,D]); end
if numel(UB) == 1; UB = repmat(UB, [1,D]); end
if size(LB,1) > 1; LB = LB'; end
if size(UB,1) > 1; UB = UB'; end
if size(widths,2) == 1; widths = repmat(widths, [1,D]); end
LB_out = LB - eps(LB);
UB_out = UB + eps(UB);
basewidths = widths;    % User-supplied widths
if isempty(options.LogPrior); doprior = 0; else doprior = 1; end

if options.Burnin == 0 && isempty(widths) && options.AdaptiveWidths
    warning('WIDTHS not specified and adaptation is ON (OPTIONS.AdaptiveWidths == 1), but OPTIONS.Burnin is set to 0. TEMPERSLICE will attempt to use default values for WIDTHS.');
end

if isempty(widths); widths = (UB - LB)/2; end
widths(isinf(widths)) = 10;
if size(x0,1) == 1; x0 = repmat(x0,[M 1]); end
if size(widths,1) == 1
    widths = bsxfun(@rdivide,widths,sqrt(betas));
end
if size(basewidths,1) == 1
    basewidths = bsxfun(@rdivide,basewidths,sqrt(betas));
end
widths = bsxfun(@min,widths,2*(UB - LB));
basewidths = bsxfun(@min,basewidths,2*(UB - LB));

% Check the starting point
assert(size(x0,1) == M, ...
    'The initial matrix X0 needs to be a row vector or a matrix with M rows, where M is the number of parallel chains.');

funccount = 0;

y = zeros(M,1);
logprior = zeros(M,1);
[y(1),temp,logprior(1)] = feval(@logpdfbound, x0(1,:), betas(1), varargin{:});   % First evaluation at X0
nData = numel(temp);    % Number of data points per sample
fval = zeros(M,nData);
fval(1,:) = temp;
for m = 2:M             % Evaluate remaining parallel chains
    [y(m),fval(m,:),logprior(m)] = feval(@logpdfbound, x0(m,:), betas(m), varargin{:});   % First evaluation at X0
end
xx = x0;
if storeallchains
    samples = zeros(N,D,M);
    if nargout > 1; fvals = zeros(N,nData,M); end
    if ~isempty(options.LogPrior) && nargout > 3; logpriors = zeros(N,M); else logpriors = []; end    
else
    samples = zeros(N,D);
    if nargout > 1; fvals = zeros(N, nData); end
    if ~isempty(options.LogPrior) && nargout > 3; logpriors = zeros(N, 1); else logpriors = []; end
end
thin = floor(options.Thin);
burn = floor(options.Burnin);
swap = floor(options.Swap);
assert(all(swap > 0) && isvector(swap), ...
    'The swap frequency OPTIONS.Swap needs to be a positive integer or a vector of positive integers.');
swap = swap(:);
if size(swap,1) == 1; swap = swap*ones(M,1); end

log_Px = y;
widths(LB == UB) = 1;   % WIDTHS is irrelevant when LB == UB, set to 1

% Sanity checks
assert(size(LB,1) == 1 && size(UB,1) == 1 && numel(LB) == D && numel(UB) == D, ...
    'LB and UB need to be empty matrices, scalars or row vectors with the same number of columns as X0.');
assert(all(isfinite(betas(:))) && isreal(betas), ...
    'The vector of inverse temperatures BETAS need to be all real numbers.');
assert(all(UB >= LB), ...
    'All upper bounds UB need to be equal or greater than lower bounds LB.');
assert(all(all(isfinite(bsxfun(@times,betas,LB)) & isfinite(bsxfun(@times,betas,UB)))), ...
    'If any BETAS is zero the corresponding bounds LB and UB need to be finite.');
assert(all(widths(:) > 0 & isfinite(widths(:))) && isreal(widths), ...
    'The matrix of typical widths WIDTHS need to be all positive real numbers.');
assert(all(all(bsxfun(@ge, x0, LB) & bsxfun(@le, x0, UB))), ...
    'The initial starting point X0 is outside the bounds.');
assert(all(isfinite(y)) && isreal(y), ...
    'The initial starting point X0 needs to evaluate to a real number (not Inf or NaN).');
assert(thin > 0 && isscalar(thin), ...
    'The thinning factor OPTIONS.Thin needs to be a positive integer.');
assert(burn >= 0 && isscalar(burn), ...
    'The burn-in samples OPTIONS.Burnin needs to be a non-negative integer.');

effN = N + (N-1)*(thin-1); % Effective samples

nswaps = zeros(1,M-1);  % Number of swaps between i-th and (i+1)-th chains

if any(strcmpi(options.Display,{'final','iter','notify'}))
    fprintf(' Iteration     f-count       log p(x)                   Action\n');
    displayFormat = ' %7.0f     %8.0f    %12.6g    %26s\n';
end

xx_sum = zeros(M,D);
xx_sqsum = zeros(M,D);

if ~isempty(options.Tunnel)
    dotunnel = 1;
    tscale = std(options.Tunnel);
    tsamples = bsxfun(@rdivide,options.Tunnel,tscale);
    tsamples2(1,:,:) = tsamples';
    dd = squeeze(bsxfunandsum(@minus,@power,tsamples,tsamples2,2,2));
    dd(dd == 0) = Inf;
    amdist = mean(min(dd)); % Average minimum distance
    tsqradius = amdist/4;
    tcount = 0;
    tPx = exp(options.TunnelF - max(options.TunnelF));
    tPx = tPx/sum(tPx);

    clear amdist dd tsamples2;
else
    dotunnel = 0; 
end


% Main loop
for ii = 1:(effN+burn)

    if any(strcmpi(options.Display,{'iter','notify'})) && ii == burn+1
        action = 'start recording';
        fprintf(displayFormat,ii-burn,funccount,log_Px(targetidx),action);
    end
    
    %% Local update via slice sampling
    
    for mm = 1:M
        
        log_uprime = log(rand) + log_Px(mm);

        % Axes sweep
        for dd = 1:D
            if LB(dd) == UB(dd); continue; end      % Fixed dimension, skip

            x_l = xx(mm,:);
            x_r = xx(mm,:);
            xprime = xx(mm,:);

            % Create a horizontal interval (x_l, x_r) enclosing xx
            rr = rand;
            x_l(dd) = xx(mm,dd) - rr*widths(mm,dd);
            x_r(dd) = xx(mm,dd) + (1-rr)*widths(mm,dd);

            % Adjust interval to outside bounds for bounded problems
            if isfinite(LB(dd)) || isfinite(UB(dd))        
                if x_l(dd) < LB_out(dd)
                    delta = LB_out(dd) - x_l(dd);
                    x_l(dd) = x_l(dd) + delta;
                    x_r(dd) = x_r(dd) + delta;
                end
                if x_r(dd) > UB_out(dd)
                    delta = x_r(dd) - UB_out(dd);
                    x_l(dd) = x_l(dd) - delta;
                    x_r(dd) = x_r(dd) - delta;
                end
                x_l(dd) = max(x_l(dd),LB_out(dd));
                x_r(dd) = min(x_r(dd),UB_out(dd));
            end

            % Step-out procedure
            if options.StepOut
                steps = 0;
                stepsize = widths(mm,dd);
                while (feval(@logpdfbound, x_l, betas(mm), varargin{:}) > log_uprime)
                    x_l(dd) = x_l(dd) - stepsize;
                    steps = steps + 1;
                end
                while (feval(@logpdfbound, x_r, betas(mm), varargin{:}) > log_uprime)
                    x_r(dd) = x_r(dd) + stepsize;
                    steps = steps + 1;
                end
                if any(strcmpi(options.Display,{'iter','notify'})) && steps >= 10
                    action = ['step-out dim ' num2str(dd) ' (' num2str(steps) ' steps)'];
                    fprintf(displayFormat,ii-burn,funccount,log_Px,action);
                end            
            end

            % Shrink procedure (inner loop)
            % Propose xprimes and shrink interval until good one found
            shrink = 0;
            while 1
                shrink = shrink + 1;

                xprime(dd) = rand()*(x_r(dd) - x_l(dd)) + x_l(dd);
                [log_Px(mm),fval(mm,:),logprior(mm)] = feval(@logpdfbound, xprime, betas(mm), varargin{:});

                if log_Px(mm) > log_uprime
                    break % this is the only way to leave the while loop
                else
                    % Shrink in
                    if xprime(dd) > xx(mm,dd)
                        x_r(dd) = xprime(dd);
                    elseif xprime(dd) < xx(mm,dd)
                        x_l(dd) = xprime(dd);
                    else
                        errorstr = ['Shrunk to current position and proposal still not acceptable. ' ...
                            'Current position: ' num2str(xx(mm,:),' %g') '. ' ...
                            'Log f: (new value) ' num2str(log_Px(mm)), ', (target value) ' num2str(log_uprime) '.'];
                        error(errorstr);                    
                    end
                end
            end

            % Width adaptation (only during burn-in, might break detailed balance)
            if ii <= burn && options.AdaptiveWidths
                delta = UB(dd) - LB(dd);
                if shrink > 3
                    if isfinite(delta)
                        widths(mm,dd) = max(widths(mm,dd)/1.1,eps(delta));
                    else
                        widths(mm,dd) = max(widths(mm,dd)/1.1,eps);                    
                    end
                elseif shrink < 2
                    widths(mm,dd) = min(widths(mm,dd)*1.2, delta);
                end
            end

            if any(strcmpi(options.Display,{'iter','notify'})) && shrink >= 10
                action = ['shrink dim ' num2str(dd) ' (' num2str(shrink) ' steps)'];
                fprintf(displayFormat,ii-burn,funccount,log_Px(mm),action);
            end

            xx(mm,dd) = xprime(dd);
        end
    end
        
    %% Global update via parallel swapping
    for mm = M-1:-1:1
        swaptime = mod(ii, max(swap(mm),swap(mm+1))) == 0;    
        if swaptime
            A = exp((betas(mm)-betas(mm+1))*(sum(fval(mm+1,:))-sum(fval(mm,:))));
            if rand() < A
                nswaps(mm) = nswaps(mm) + 1;
                temp = xx(mm,:);
                xx(mm,:) = xx(mm+1,:);
                xx(mm+1,:) = temp;
                temp = fval(mm,:);
                fval(mm,:) = fval(mm+1,:);
                fval(mm+1,:) = temp;
                temp = logprior(mm);
                logprior(mm)= logprior(mm+1);
                logprior(mm+1) = temp;
                % Update current function value
                log_Px([mm;mm+1]) = betas([mm;mm+1]).*sum(fval([mm,mm+1],:),2) ...
                    + logprior([mm;mm+1]);
            end
        end
    end
    
    %% Global update via tunnelling
    if dotunnel
       for mm = 1:M
           dist = sum(bsxfun(@minus, tsamples, xx(mm,:)./tscale).^2,2);
           f = dist < tsqradius;
           idx = find(f);
           nx = numel(idx);
           if isempty(idx); continue; end
           v = xx(mm,:)./tscale - tsamples(idx(randi(nx)),:);
           % Pick a random tunnel target
           % tidx = find(rand() < cumsum(tPx),1,'last');
           idx_prime = randi(size(tsamples,1));
           xprime = (tsamples(idx_prime,:) + v).*tscale;
           [log_uprime,fval_prime,logprior_prime] = feval(@logpdfbound, xprime, betas(mm), varargin{:});
           dist = sum(bsxfun(@minus, tsamples, xprime./tscale).^2,2);
           nprime = sum(dist < tsqradius);
           A = exp(log_uprime - log_Px(mm))*nx/nprime;
           % [A, xprime]
           if rand() < A
                xx(mm,:) = xprime;
                log_Px(mm) = log_uprime;
                fval(mm,:) = fval_prime;
                logprior(mm) = logprior_prime;           
           end               
           tcount = tcount + f;
       end
    end
    
    %% Record samples and miscellaneous bookkeeping
    
    % Record samples?
    record = ii > burn && mod(ii - burn - 1, thin) == 0;
    if record
        ismpl = 1 + (ii - burn - 1)/thin;
        if storeallchains           % Store all chains while running
            samples(ismpl,:,:) = xx';
            if nargout > 1; fvals(ismpl,:,:) = fval'; end
            if nargout > 3 && doprior; logpriors(ismpl,:) = logprior; end            
        else        
            samples(ismpl,:) = xx(targetidx,:);
            if nargout > 1; fvals(ismpl,:) = fval(targetidx,:); end
            if nargout > 3 && doprior; logpriors(ismpl) = logprior(targetidx); end
        end
    end
    
    % Store summary statistics during burn-in
    if ii <= burn
        xx_sum = xx_sum + xx;
        xx_sqsum = xx_sqsum + xx.^2;
        
        % End of burn-in, update WIDTHS if using adaptive widths method
        if ii == burn && options.AdaptiveWidths
            newwidths = 10*sqrt(xx_sqsum/burn - (xx_sum/burn).^2);
            newwidths = bsxfun(@min, newwidths, UB_out - LB_out);
            if isempty(basewidths)
                widths = newwidths;
            else
                % Max between new widths and geometric mean with user-supplied 
                % widths (i.e. bias towards keeping larger widths)
                widths = bsxfun(@max, newwidths, sqrt(newwidths.*basewidths));
            end
        end
    end
        
    if strcmpi(options.Display,'iter')
        if ii <= burn; action = 'burn';
        elseif ~record; action = 'thin';
        else action = 'record';
        end
        fprintf(displayFormat,ii-burn,funccount,log_Px(targetidx),action);
    end
    
end

if any(strcmpi(options.Display,{'iter','final','notify'}))
    if thin > 1
        thinmsg = ['\n   and keeping 1 sample every ' num2str(thin)];
    else
        thinmsg = '';
    end    
    fprintf(['\nSampling terminated:\n * %d samples obtained after a burn-in period of %d samples' thinmsg ',\n   for a total of %d function evaluations.'], ...
        N, burn, funccount);
end

% Diagnostics
if options.Diagnostics && ...
        (nargout > 2 || any(strcmpi(options.Display,{'iter','final','notify'})))
    for m = 1:M; [exitflag(m),R(m,:),Neff(m,:),tau(m,:)] = diagnose(samples,m,options); end
    diagstr = [];
    if any(exitflag == 1) || any(exitflag == 2)
        diagstr = [diagstr '\n * Try sampling for longer, by increasing N or OPTIONS.Thin, or try increasing mixing with more parallel chains.'];
    elseif any(exitflag == 3)
        diagstr = [diagstr '\n * Try increasing OPTIONS.Thin to obtain more uncorrelated samples.'];
    elseif all(exitflag == 0)
        diagstr = '\n * No violations of convergence have been detected (this does NOT guarantee convergence).';
    end
    if any(strcmpi(options.Display,{'iter','final','notify'})) && ~isempty(diagstr)
        fprintf(diagstr);
    end
    
    exitflag = exitflag(targetidx);
    
end

% Return only target chain
if ~options.Debug
    samples = samples(:,:,targetidx);
    fvals = fvals(:,:,targetidx);
    logpriors = logpriors(:,targetidx);
end

if any(strcmpi(options.Display,{'iter','final','notify'}))
    fprintf('\n\n');
end

if nargout > 3
    output.widths = widths;
    output.logpriors = logpriors;
    output.funccount = funccount;
    output.nswaps = nswaps;
    if options.Diagnostics
        output.R = R;
        output.Neff = Neff;
        output.tau = tau;
    end
    if dotunnel
        output.tcount = tcount;
    end
end

%--------------------------------------------------------------------------
function [y,fval,logprior] = logpdfbound(x,beta,varargin)
%LOGPDFBOUND Evaluate log pdf with bounds and prior.

    fval = [];
    logprior = [];
    
    if any(x < LB | x > UB)
        y = -Inf;
    else
        
        if doprior
            logprior = feval(options.LogPrior, x);
            if isnan(logprior)
                y = -Inf;
                warning('Prior density function returned NaN. Trying to continue.');
                return;
            elseif ~isfinite(logprior)
                y = -Inf;
                return;
            end
        else
            logprior = 0;
        end
        
        fval = logf(x,varargin{:});
        funccount = funccount + 1;

        if any(isnan(fval))
            y = -Inf;
            warning('Target density function returned NaN. Trying to continue.');
        else
            y = beta*sum(fval) + logprior;        
        end
    end
    
end

end

function [exitflag,R,Neff,tau] = diagnose(samples,m,options)

    N = size(samples,1);
    exitflag = -1;    

    try        
        warning_orig = warning;
        warning('off','all');
        [R,Neff,~,~,~,tau] = psrf(samples(1:floor(N/2),:,m), samples(floor(N/2)+(1:floor(N/2)),:,m));
        warning(warning_orig);
        
        diagstr = [];
        if any(R > 1.5)
            diagstr = ['\n * Chain ' num2str(m) ': Detected lack of convergence! (max R = ' num2str(max(R),'%.2f') ' >> 1, mean R = ' num2str(mean(R),'%.2f') ').'];
            exitflag = 1;
        elseif any(R > 1.1)
            diagstr = ['\n * Chain ' num2str(m) ': Detected probable lack of convergence (max R = ' num2str(max(R),'%.2f') ' > 1, mean R = ' num2str(mean(R),'%.2f') ').'];
            exitflag = 2;
        end
        if any(Neff < N/10)
            diagstr = [diagstr '\n * Chain ' num2str(m) ': Low number of effective samples (min Neff = ' num2str(min(Neff), '%.1f') ...
                ', mean Neff = ' num2str(mean(Neff),'%.1f') ', requested N = ' num2str(N,'%d') ').'];
            if exitflag == 0; exitflag = 3; end
        end
        if isempty(diagstr) && exitflag == -1
            exitflag = 0;
        end
        
        if any(strcmpi(options.Display,{'iter','final','notify'})) && ~isempty(diagstr)
            fprintf(diagstr);
        end
        
    catch
        warning('Error while computing convergence diagnostics with PSRF.');
        R = NaN;
        Neff = NaN;
    end

end
