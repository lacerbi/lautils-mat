function [x,fval,exitflag,output,cache] = fminmulti(fun,LB,UB,nStarts,options,vararginarray)
%FMINMULTI Function minimization via recursive multi-start method.
%   X = FMINMULTI(FUN,LB,UB,NSTARTS) finds a global minimum X of the 
%   function FUN within lower bounds LB and upper bounds UB. FUN accepts a 
%   row vector X as input and returns a scalar function value F evaluated 
%   at X. LB and UB are row vectors of the same length as the input 
%   accepted by FUN. NSTARTS is an array with the number of starting points 
%   per each minimization cycle. At each cycle i, FMINMULTI starts 
%   NSTARTS(i) local optimizations from the NSTARTS(i) best points from the 
%   previous minimization cycle (or randomly chosen points within the 
%   optimization bounds).
%
%   X = FMINMULTI(FUN,LB,UB,NSTARTS,OPTIONS) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS. 
%   The accepted fields for OPTIONS are:
%
%FMINMULTI PARAMETERS for MATLAB
%Method (*)      - Method for local optimization
%                   [ 'bps' | {'fminsearch'} | 'fmincon' | 'feval' | 'patternsearch' ]
%MaxFunEvals (*) - Maximum number of function evaluations allowed
%                   [ positive integer | {200*numberofvariables} ]
%MaxIter (*)     - Maximum number of iterations allowed
%                   [ positive integer | {200*numberofvariables} ]
%TolFun (*)      - Termination tolerance on the function value
%                   [ positive scalar | {1e-4} ]
%TolX (*)        - Termination tolerance on X
%                   [ positive scalar | {1e-6} ]
%OptimOptions (*)- Local optimizer options struct (fields in 'OptimOptions' 
%                  ovverride global settings in OPTIONS)
%                   [ options struct | {[]} ]
%Display (*)     - Level of display
%                   [ {'off'} | 'iter' | 'notify' | 'final' ]
%OutputFcn      - Output function called by FMINMULTI at each iteration 
%                   [ function handle | {[]} ] 
%InitRange      - Initial/reasonable range of values for optimization;
%                 first row lower bounds, second row upper bounds.
%                 'InitRange' should be specified as the region that 
%                 is likely to contain the global minimum (determined via 
%                 educated guesses and prior knowledge). This region is 
%                 usually smaller than the region bracketed by LB and UB.
%                   [ matrix | {[LB; UB]} ]
%InitialPoints  - Matrix array of starting points (each as a row vector)
%                   [ matrix | {[]} ]
%MidpointStart  - Add midpoint to starting set
%                   [ {'on'} | 'off' ]
%XScale         - Typical length scale 
%                 (used to judge closeness between points)
%                   [ vector | {diff(InitialRange,[],1)/sqrt(12)} ]
%FvalScale      - Difference between function values considered significant
%                 (used to judge closeness between points)
%                   [ positive scalar | {1} ]
%RescaleVars    - Linear/nonlinear rescaling of variables to improve search
%                 (variables are automatically rescaled back to normal
%                 coordinates - however be aware that there might be 
%                 roundoff errors which may affect equality comparisons)
%                   [ 'on' | {'off'} ]
%Cache          - Keep a CACHE of the points evaluated so far to decide 
%                 future starting points. This options could be expensive 
%                 in terms of memory and speed.
%                   [ {'on'}, 'off' ]
%CacheSize      - Limit the CACHE size. A typical choice of 1e5 is 
%                 suggested but it depends on computer speed and memory. 
%                 Used if 'Cache' option is 'on'.
%                   [ positive integer | {1e5} ]
%
%   (*) These options accept a matrix array or cell array as value, in
%   which case the i-th element of the array is used in the i-th
%   optimization cycle (the last element of the array is used if the
%   number of cycles is greater than the array length).
%
%   Example
%      A typical usage of FMINMULTI:
%         camel6 = @(x) (4-2.1*x(1)^2+(x(1)^4)/3) * x(1)^2 + x(1)*x(2) + (-4+4*x(2)^2) * x(2)^2;
%         localopt.Algorithm = 'sqp';
%         options.Method = {'feval','fmincon','patternsearch'};
%         options.OptOptions = {[],localopt,[]};
%         [x,fval] = fminmulti(camel6,[-3,-2],[3,2],[100,10,1],options)
%      optimizes the 'six-hump camel' function by first evaluating it
%      at 100 random starting points; then it runs FMINCON with the
%      'sqp' algorithm starting from the 10 most promising (lowest-value), 
%      non-overlapping points; finally it takes the best point so far
%      and finalizes the optimization with a run of PATTERNSEARCH. (This
%      is overkilling for the camel6 function.)
%
%   X = FMINMULTI(FUN,LB,UB,NSTARTS,OPTIONS,VARARGINARRAY) passes a cell
%   array of additional arguments to the objective function. VARARGINARRAY
%   is a cell array of cell arrays. VARARGINARRAY{1}{:} is passed as an 
%   additional list of arguments to FUN in the 1st cycle,
%   VARARGINARRAY{2}{:} is passed to FUN on the 2nd cycle, and so on.
%   VARARGINARRAY{end}{:} is passed to FUN for any cycle greater than the
%   length of VARARGINARRAY.
%
%   [X,FVAL] = FMINMULTI(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINMULTI(...) returns an EXITFLAG that describes 
%   the exit condition of FMINMULTI. Currently EXITFLAG is always 0. This
%   output is used for compatibility with other minimization functions.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINMULTI(...) returns a structure 
%   OUTPUT with information such as total number of iterations, and final 
%   objective function value.
%
%   [X,FVAL,EXITFLAG,OUTPUT,CACHE] = FMINMULTI(...) returns the content
%   of the CACHE structure (the 'Cache' option must be 'on').
%
%   FMINMULTI uses FMINSEARCHBND to perform bounded optimization with
%   the Nelder-Mead algorithm since the standard FMINSEARCH does not
%   support optimization bounds.
%
%   See also FMINCON, FMINSEARCH, FMINSEARCHBND, PATTERNSEARCH.

% Author: Luigi Acerbi
% e-mail: luigi.acerbi@gmail.com
% Release: 0.9
% Release date: 17/Aug/2015

if nargin < 5; options = []; end
if nargin < 6; vararginarray = []; end

LB = LB(:)'; UB = UB(:)';   % Lower and upper bounds as row vectors
nvars = length(LB);         % Number of parameters

if ~isfield(options,'InitRange') || isempty(options.InitRange)
    warning('''InitRange'' field of OPTIONS structure not specified. Using full LB and UB range.');
end

defopts.Method = 'fminsearch';          % Optimization method
defopts.MaxFunEvals = 200*nvars;        % Max function evaluations
defopts.MaxIter = 200*nvars;            % Max iterations
defopts.TolFun = 1e-4;                  % TolFun
defopts.TolX = 1e-6;                    % TolX
defopts.OptimOptions = [];              % Local optimizer options
defopts.Display = 'off';                % Display level
defopts.OutputFcn = [];                 % Output function
defopts.InitRange = [LB; UB];           % Initial range
defopts.InitialPoints = [];             % Starting minimization
defopts.SobolInit = 'on';               % Initialize with Sobol grid
defopts.SobolSeed = [];                 % Chosen seed for Sobol grid
defopts.MidpointStart = 'on';           % Include midpoint
defopts.XScale = [];                    % Characteristic length scale in X
defopts.FvalScale = 1;                  % Characteristic function amplitude
defopts.RescaleVars = 'off';            % Rescale variables
defopts.Cache = 'on';                   % Cache intermediate fcn evaluations
defopts.CacheSize = 1e5;                % Cache size
defopts.BPSUseCacheEpochs = 2;          % Use cached fcn evals for BPS
defopts.LoadFile    = '';               % Load interrupted sampling from file
defopts.SaveFile    = '';               % Save sampling to file
defopts.SaveTime    = 1e4;              % Save every this number of seconds

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

% Cellify fields
for f = {'Method','Display','OptimOptions'}
    if ~iscell(options.(f{:})); options.(f{:}) = {options.(f{:})}; end
end
if ~isempty(vararginarray) && ~iscell(vararginarray); vararginarray = {vararginarray}; end

% Matrixfy fields
for f = {'MaxFunEvals','MaxIter','TolFun','TolX'}
    if iscell(options.(f{:})); options.(f{:}) = cell2mat(options.(f{:})); end
end

% Convert string, numerical or logical variables to 'on' and 'off'
for f = {'MidpointStart','RescaleVars','Cache'}
    options.(f{:}) = onoff(options.(f{:}));
end

if ~isempty(options.LoadFile) && exist(options.LoadFile,'file')
    %% Load interrupted execution from file
    
    optionsnew = options;           % Store current options    
    load(options.LoadFile,'-mat');  % Load run from recovery file
    % Start epoch and inner loop from recovered iteration
    iEpoch0 = iEpoch;
    iStart0 = iStart;
    % if trace > 1; 
        fprintf('Loading ongoing optimization from file ''%s''.\n', options.LoadFile); 
    % end

    % Copy some new OPTIONS to the old options structure
    copyfields = {'LoadFile','SaveFile','SaveTime','Display'};    
    for f = 1:numel(copyfields)    
        options.(copyfields{f}) = optionsnew.(copyfields{f});
    end
    clear copyfields f optionsnew;
else

    % Reasonable lower and upper bounds
    RLB = options.InitRange(1,:); RUB = options.InitRange(2,:);
    if isscalar(RLB); RLB = RLB*ones(1,nvars); end
    if isscalar(RUB); RUB = RUB*ones(1,nvars); end
    RLB = max(RLB,LB); RUB = min(RUB,UB);

    % Characteristic length scale from reasonable range
    if isempty(options.XScale); options.XScale = (RUB - RLB)/sqrt(12); end
    
    x0 = options.InitialPoints;
    
    % Sobol initialization
    if size(x0,1) < nStarts(1) && strcmpi(options.SobolInit,'on')
        seed = options.SobolSeed;
        if isempty(seed); seed = randi(1e4); end
        r = i4_sobol_generate(nvars,nStarts(1)-size(x0,1),seed)';
        x0 = [x0; bsxfun(@plus, RLB, bsxfun(@times, r, RUB-RLB))];
    end
    
    maxEpochs = length(nStarts);
    output.x = []; output.startx = []; output.nruns = 0;
    x = NaN(1,nvars);
    fvalmin = Inf;

    if ~isempty(options.OutputFcn)
        optimValues.funccount = 0;
        optimValues.x = x;
        optimValues.fval = fvalmin;
        optimValues.epoch = 0;
        optimValues.iteration = 0;
        optimValues.procedure = [];
        stop = options.OutputFcn(x0,optimValues,'init');
    end

    % Caching intermediate points
    if strcmpi(options.Cache,'on')
        cache.x = NaN(options.CacheSize,nvars);
        cache.fval = NaN(options.CacheSize,1);
        cache.index = 1;
    else
        cache.x = [];
        cache.fval = [];
        cache.index = 0;
    end

    % Rescale variables if requested
    if strcmp(options.RescaleVars,'on')
        octstruct = optimct(nvars,LB,UB,RLB,RUB);
        LB = octstruct.lb(:)';
        UB = octstruct.ub(:)';
        RLB = octstruct.plb(:)';
        RUB = octstruct.pub(:)';
        XScale = 2*octstruct.gamma(:)'.*(options.XScale./diff(options.InitRange,[],1));
        options.InitRange = [octstruct.plb(:)'; octstruct.pub(:)'];
        if ~isempty(x0)
            for i = 1:size(x0,1); x0(i,:) = optimct(x0(i,:),octstruct); end
        end
    else
        octstruct = [];
    end
    
    iEpoch0 = 1;
    iStart0 = 0;
end

lastsave = tic;    % Keep track of time

% Loop over optimization epochs
for iEpoch = iEpoch0:maxEpochs
    
    if iStart0 == 0
        % Initial points
        if isempty(x0) && strcmpi(options.MidpointStart,'on') % Reasonable starting point
            x0(1,:) = 0.5*(RLB+RUB);
        end

        if size(x0,1) < nStarts(iEpoch)
            index = size(x0,1)+1:nStarts(iEpoch);
            % x0(index,:) = bsxfun(@plus,bsxfun(@times,rand(length(index),nvars),RUB-RLB),RLB);
            temp = cache.x(~isnan(cache.fval),:);
            X = [output.x;output.startx;temp;x0(1:size(x0,1),:)];  
            x0(index,:) = randx(length(index),RLB,RUB,options.XScale/2,X);
        end
        x0 = x0(1:nStarts(iEpoch),:);

        maxeval = options.MaxFunEvals(min(iEpoch,end));
        maxiter = options.MaxIter(min(iEpoch,end));
        tolfun = options.TolFun(min(iEpoch,end));
        tolx = options.TolX(min(iEpoch,end));
        disp = options.Display{min(iEpoch,end)};
        method = options.Method{min(iEpoch,end)};

        if ~isempty(options.OutputFcn)
            optimValues.funccount = 0;
            optimValues.x = NaN(1,nvars);
            optimValues.fval = fvalmin;
            optimValues.epoch = iEpoch;
            optimValues.iteration = 0;
            optimValues.procedure = method;
            stop = options.OutputFcn(x0,optimValues,'iter');
        end
        
        iStart0 = 1;
    end
    
    bpscache = [];
        
    % Inner optimization loop
    for iStart = iStart0:nStarts(iEpoch)
                        
        % Save current progress to file
        if ~isempty(options.SaveFile) && toc(lastsave) > options.SaveTime
            % if trace > 1;
                fprintf('Saving temp file ''%s''...\n', options.SaveFile); 
            % end
            save(options.SaveFile);
            lastsave = tic;
        end
        
        if iEpoch > 1; display([iEpoch,iStart]); end
        
        xstart = x0(iStart,:);
        
        if ~isempty(options.OutputFcn)
            optimValues.funccount = 0;
            optimValues.x = optimct(x,octstruct,1);
            optimValues.fval = fvalmin;
            optimValues.epoch = iEpoch;
            optimValues.iteration = iStart;
            optimValues.procedure = method;
            stop = options.OutputFcn(optimct(xstart,octstruct,1),optimValues,'iter');
        end        

        if cache.index
            % Internal objective function that caches fcn evaluations
            funfcn = @objectivefunc;
        else
            % Pass additional arguments to function            
            if isempty(vararginarray)
                funfcn = fun;
            else
                funfcn = @(x_) fun(x_, vararginarray{min(iEpoch,end)}{:});
            end
        end
        
        switch method
            case 'feval'
                fval = funfcn(xstart);
                x = xstart; exitflag = 0; output1 = [];
                funcCount = 1;
            case 'fminsearch'
                optoptions = optimset('Display',disp,'TolX',tolx,'TolFun',tolfun,'MaxFunEvals',maxeval,'MaxIter',maxiter);
                optoptions = optimset(optoptions,options.OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = fminsearchbnd(funfcn,xstart,LB,UB,optoptions);
                funcCount = output1.funcCount;
            case 'fmincon'
                assert(exist('patternsearch.m','file') == 2, ...
                    'FMINCON optimization method not installed (requires Optimization Toolbox).');
                optoptions = optimset('Display',disp,'TolX',tolx,'TolFun',tolfun,'MaxFunEvals',maxeval,'MaxIter',maxiter);
                optoptions = optimset(optoptions,options.OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = fmincon(funfcn,xstart,[],[],[],[],LB,UB,[],optoptions);
                funcCount = output1.funcCount;
            case 'patternsearch'
                assert(exist('patternsearch.m','file') == 2, ...
                    'PATTERNSEARCH optimization method not installed (requires Global Optimization Toolbox).');
                optoptions = psoptimset('Display',disp,'TolX',tolx,'TolFun',tolfun,'MaxFunEvals',maxeval,'MaxIter',maxiter);
                optoptions = psoptimset(optoptions,options.OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = patternsearch(funfcn,xstart,[],[],[],[],LB,UB,[],optoptions);
                funcCount = output1.funccount;
            case 'bps'
                assert(exist('bps.m','file') == 2, ...
                    'BPS optimization method not installed.');
                optoptions = options.OptimOptions{min(iEpoch,end)};
                optoptions.Display = disp;
                optoptions.TolX = tolx;
                optoptions.TolFun = tolfun;
                optoptions.MaxFunEvals = maxeval;
                optoptions.MaxIter = maxiter;
                % Feed cached function values to BPS
                if any(options.BPSUseCacheEpochs == iEpoch) && isempty(bpscache)
                    bpscache.X = cache.x(1:cache.index,:);
                    bpscache.Y = cache.fval(1:cache.index,:);
                end
                if any(options.BPSUseCacheEpochs == iEpoch)
                    optoptions.FunValues = bpscache;
                end
                PLB = options.InitRange(1,:);
                PUB = options.InitRange(2,:);
                [x,fval,exitflag,output1] = bps(funfcn,xstart,LB,UB,PLB,PUB,optoptions);
                funcCount = output1.FuncCount;                
        end

        % Store optimization results
        output = storeOptimization(output,...
            xstart,x,fval,funcCount,exitflag,output1);
        
        if ~isempty(options.OutputFcn)
            optimValues.funccount = funcCount;
            optimValues.x = optimct(x,octstruct,1);
            optimValues.fval = fval;
            optimValues.epoch = iEpoch;
            optimValues.iteration = iStart;
            optimValues.procedure = method;
            stop = options.OutputFcn(optimct(x,octstruct,1),optimValues,'iter');
        end
                
        if fval < fvalmin
            fvalmin = fval;
            xmin = x;
        end        
        
    end

    % Recompute x0 for next iteration
    if iEpoch < maxEpochs && size(output.x,1) > 1
        % Build optimization matrix combining optimal points and log likelihood
        xScale = options.XScale(:)'/10*sqrt(nvars);
        fval = output.fval(:);
        xx = [bsxfun(@rdivide,output.x,xScale), fval/options.FvalScale];
        
        % Cluster points close in x-Fval space, keep only minima
        listindex = 1:size(fval,1);        
        for i = 1:size(xx,1)
            if ~isnan(fval(i))
                dd = sum(bsxfun(@minus, xx, xx(i,:)).^2,2);
                clusterindex = find(dd < 1);
                [~,index] = min(fval(clusterindex));
                temp = listindex(clusterindex);
                clusterindex = setdiff(clusterindex,temp(index));
                xx(clusterindex,:) = NaN;
                fval(clusterindex) = NaN;
            end
        end
                
        % Order remaining points and remove discarded ones
        [fval,index] = sort(fval,'ascend');
        x0 = output.x(index,:);
        x0(isnan(fval),:) = [];
    else
        x0 = output.x;        
    end
    
    iStart0 = 0;

end

x = optimct(xmin,octstruct,1);
fval = fvalmin;
exitflag = 0;

if ~isempty(options.OutputFcn)
    optimValues.funccount = funcCount;
    optimValues.x = optimct(x,octstruct,1);
    optimValues.fval = fval;
    optimValues.epoch = iEpoch;
    optimValues.iteration = iStart;
    optimValues.procedure = [];
    stop = options.OutputFcn(optimValues.x,optimValues,'done');
end

% Convert function outputs to original coordinate system
if strcmp(options.RescaleVars,'on')
    if nargout > 3
        for i = 1:size(output.x,1)
            output.x(i,:) = optimct(output.x(i,:),octstruct,1);
        end
        for i = 1:size(output.startx,1)
            output.startx(i,:) = optimct(output.startx(i,:),octstruct,1);
        end
    end

    if nargout > 4
        for i = 1:size(options.Cache.x,1)
            options.Cache.x(i,:) = optimct(options.Cache.x(i,:),octstruct,1);
        end
    end
end
     
    %OBJECTIVEFUNC Objective function with caching
    function fval = objectivefunc(x)
                
        if strcmp(options.RescaleVars,'on'); x = optimct(x,octstruct,1); end
        
        if isempty(vararginarray)
            fval = fun(x);
        else
            fval = fun(x, vararginarray{min(iEpoch,end)}{:});
        end
        
        cache.x(cache.index,:) = x;
        cache.fval(cache.index) = fval;
        cache.index = max(mod(cache.index + 1,options.CacheSize+1), 1);
    end

end

%--------------------------------------------------------------------------
function optimization = storeOptimization(optimization,startx,x,fval,funcalls,exitflag,output)
%STOREOPTIMIZATION Save results from optimization in struct

    if isempty(optimization)
        optimization.nruns = 1;
    else
        optimization.nruns = optimization.nruns + 1;
    end
    index = optimization.nruns;
    optimization.startx(index, :) = startx(:)';
    optimization.x(index, :) = x(:)';
    optimization.fval(index) = fval;
    optimization.funcalls(index) = funcalls;
    optimization.exitflag(index) = exitflag;
    optimization.output{index} = output;
end

%--------------------------------------------------------------------------
function y = onoff(x)
%ONOFF Convert a char, numerical or logical input to 'on'/'off'

flag = 0;   % Error flag

if isscalar(x) && (islogical(x) || isnumeric(x))
    if x; y = 'on'; else y = 'off'; end;
elseif ischar(x)
    switch lower(x)
        case 'on'; y = 'on';
        case 'off'; y = 'off';
        otherwise; flag = 1;
    end
else
    flag = 1;
end

if flag
    error('OPTION struct value should be ''on'' or ''off''.');    
end

end

%--------------------------------------------------------------------------
function varargout = optimct(varargin)
%OPTIMCT Coordinate transform for guided optimization.

NumEps = 1e-10; % Accepted numerical error

% MaxPrecision = 17; % Maximum precision for a double

if nargin >= 5

    nvars = varargin{1};
    lb = varargin{2}(:);
    ub = varargin{3}(:);
    plb = varargin{4}(:);
    pub = varargin{5}(:);
    if nargin < 6; octstruct.logct = []; else octstruct.logct = varargin{6}(:); end

    % Convert scalar inputs to column vectors
    if isscalar(lb); lb = lb*ones(nvars,1); end
    if isscalar(ub); ub = ub*ones(nvars,1); end
    if isscalar(plb); plb = plb*ones(nvars,1); end
    if isscalar(pub); pub = pub*ones(nvars,1); end

    assert(~any(isinf(([plb; pub]))), ...
        'Plausible interval ranges PLB and PUB need to be finite.');

    plb = max(plb,lb);
    pub = min(pub,ub);
    
    if isempty(octstruct.logct)
        % A variable is converted to log scale if all bounds are positive and 
        % the plausible range spans at least one order of magnitude
        octstruct.logct = all([lb, ub, plb, pub] > 0, 2) & (pub./plb >= 10);    
    elseif isscalar(octstruct.logct)
        octstruct.logct = ones(nvars,1);
    end
    
    % Transform to log coordinates
    octstruct.oldbounds.lb = lb;
    octstruct.oldbounds.ub = ub;
    octstruct.oldbounds.plb = plb;
    octstruct.oldbounds.pub = pub;    
    octstruct.lb = lb; octstruct.ub = ub; octstruct.plb = plb; octstruct.pub = pub;
    octstruct.lb(octstruct.logct) = log(octstruct.lb(octstruct.logct));
    octstruct.ub(octstruct.logct) = log(octstruct.ub(octstruct.logct));
    octstruct.plb(octstruct.logct) = log(octstruct.plb(octstruct.logct));
    octstruct.pub(octstruct.logct) = log(octstruct.pub(octstruct.logct));

    octstruct.mu = 0.5*(octstruct.plb + octstruct.pub);
    octstruct.gamma = 0.5*(octstruct.pub - octstruct.plb);
        
    z = ['( (x - ' vec2str(octstruct.mu) ') ./ ' vec2str(octstruct.gamma) ' )'];
    zlog = ['( (log(abs(x) + (x == 0)) - ' vec2str(octstruct.mu) ') ./ ' vec2str(octstruct.gamma) ' )'];
    
    switch sum(octstruct.logct)
        case 0
            octstruct.g = ['@(x) linlog(' z ')'];
            octstruct.ginv = ['@(y) ' vec2str(octstruct.gamma) ' .* linexp(y) + ' vec2str(octstruct.mu) ];        
            
        case nvars
            octstruct.g = ['@(x) linlog(' zlog ')'];
            octstruct.ginv = ['@(y) exp(' vec2str(octstruct.gamma) ' .* linexp(y) + ' vec2str(octstruct.mu) ')'];        
            
        otherwise
            octstruct.g = ['@(x) (1-' vec2str(octstruct.logct) ').*linlog(' z ')' ... 
                '+ ' vec2str(octstruct.logct) '.*linlog(' zlog ')'];
            octstruct.ginv = ['@(y) (1-' vec2str(octstruct.logct) ') .* (' vec2str(octstruct.gamma) ' .* linexp(y) + ' vec2str(octstruct.mu) ') + ' ...
                vec2str(octstruct.logct) ' .* exp(' vec2str(octstruct.gamma) ' .* linexp(y) + ' vec2str(octstruct.mu) ')'];        
    end
    
    % Convert starting values to transformed coordinates
    octstruct.g = str2func(octstruct.g);
    octstruct.ginv = str2func(octstruct.ginv);
    
    if ischar(octstruct.g); g = str2func(octstruct.g); else g = octstruct.g; end
    if ischar(octstruct.ginv); ginv = str2func(octstruct.ginv); else ginv = octstruct.ginv; end
    
    % Check that the transform works correctly in the range
    t(1) = all(abs(ginv(g(lb)) - lb) < NumEps);
    t(2) = all(abs(ginv(g(ub)) - ub) < NumEps);
    t(3) = all(abs(ginv(g(plb)) - plb) < NumEps);
    t(4) = all(abs(ginv(g(pub)) - pub) < NumEps);    
    assert(all(t), 'Cannot invert the transform to obtain the identity at the provided boundaries.');
    
    octstruct.lb = g(lb);
    octstruct.ub = g(ub);
    octstruct.plb = g(plb);
    octstruct.pub = g(pub);
    
    varargout{1} = octstruct;

elseif nargin >= 2
    
    octstruct = varargin{2};    

    if isempty(octstruct)
        varargout{1} = varargin{1}; % Return untransformed input
    else
        if (nargin < 3); direction = 'd'; else direction = varargin{3}(1); end
        
        if direction == 'd' || direction == 'D'
            x = varargin{1};
            if ischar(octstruct.g); g = str2func(octstruct.g); else g = octstruct.g; end
            y = g(x(:));
            y = min(max(y,octstruct.lb),octstruct.ub);    % Force to stay within bounds

            varargout{1} = reshape(y,size(x));
        else
            y = varargin{1};
            if ischar(octstruct.ginv); ginv = str2func(octstruct.ginv); else ginv = octstruct.ginv; end
            x = ginv(y(:));
            x = min(max(x,octstruct.oldbounds.lb),octstruct.oldbounds.ub);    % Force to stay within bounds
            varargout{1} = reshape(x,size(y));
        end
    end
    
end

end

%--------------------------------------------------------------------------
function s = vec2str(v)
% Convert numerical vector to string

MaxPrecision = 17;  % Maximum precision for a double
if size(v,1) > 1; transp = ''''; else transp = []; end
s = '[';
for i = 1:length(v)-1; s = [s, num2str(v(i),MaxPrecision), ',']; end
s = [s, num2str(v(end),MaxPrecision), ']' transp];

% s = ['[' num2str(v(:)',MaxPrecision) ']' transp];

end



