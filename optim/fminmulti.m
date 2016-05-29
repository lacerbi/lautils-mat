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

% Default options
Method = cgetfield(options,'Method',{'fminsearch'});        % Optimization method
MaxFunEvals = cgetfield(options,'MaxFunEvals',200*nvars);   % Max function evaluations
MaxIter = cgetfield(options,'MaxIter',200*nvars);           % Max iterations
TolFun = cgetfield(options,'TolFun',1e-4);                  % TolFun
TolX = cgetfield(options,'TolX',1e-6);                      % TolX
OptimOptions = cgetfield(options,'OptimOptions',[]);        % Local optimizer options
Display = cgetfield(options,'Display','off');               % Display level
OutputFcn = cgetfield(options,'OutputFcn',[]);              % Output function
InitRange = cgetfield(options,'InitRange',[LB;UB]);         % Initial range
x0 = cgetfield(options,'InitialPoints',[]);                 % Starting minimization
MidpointStart = cgetfield(options,'MidpointStart','on');    % Include midpoint
XScale = cgetfield(options,'XScale',[]);                    % Characteristic length scale in X
FvalScale = cgetfield(options,'FvalScale',1);               % Characteristic function amplitude
RescaleVars = cgetfield(options,'RescaleVars','off');       % Rescale variables
Cache = cgetfield(options,'Chache','on');                   % Cache intermediate fcn evaluations
CacheSize = cgetfield(options,'ChacheSize',1e5);            % Cache size

if ~iscell(Method); Method = {Method}; end
if ~iscell(Display); Display = {Display}; end
if ~iscell(OptimOptions); OptimOptions = {OptimOptions}; end
if ~isempty(vararginarray) && ~iscell(vararginarray); vararginarray = {vararginarray}; end

if iscell(MaxFunEvals); MaxFunEvals = cell2mat(MaxFunEvals); end
if iscell(MaxIter); MaxIter = cell2mat(MaxIter); end
if iscell(TolFun); TolFun = cell2mat(TolFun); end
if iscell(TolX); TolX = cell2mat(TolX); end

% Convert string, numerical or logical variables to 'on' and 'off'
MidpointStart = onoff(MidpointStart);
RescaleVars = onoff(RescaleVars);
Cache = onoff(Cache);

% Reasonable lower and upper bounds
RLB = InitRange(1,:); RUB = InitRange(2,:);
if isscalar(RLB); RLB = RLB*ones(1,nvars); end
if isscalar(RUB); RUB = RUB*ones(1,nvars); end
RLB = max(RLB,LB); RUB = min(RUB,UB);

% Characteristic length scale from reasonable range
if isempty(XScale); XScale = (RUB - RLB)/sqrt(12); end

maxEpochs = length(nStarts);
output.x = []; output.startx = []; output.nruns = 0;
x = NaN(1,nvars);
fvalmin = Inf;

if ~isempty(OutputFcn)
    optimValues.funccount = 0;
    optimValues.x = x;
    optimValues.fval = fvalmin;
    optimValues.epoch = 0;
    optimValues.iteration = 0;
    optimValues.procedure = [];
    stop = OutputFcn(x0,optimValues,'init');
end

% Caching intermediate points
if strcmpi(Cache,'on')
    cache.x = NaN(CacheSize,nvars);
    cache.fval = NaN(CacheSize,1);
    cache.index = 1;
else
    cache.x = [];
    cache.fval = [];
    cache.index = 0;
end


if strcmp(RescaleVars,'on')
    octstruct = optimct(nvars,LB,UB,RLB,RUB);
    LB = octstruct.lb(:)';
    UB = octstruct.ub(:)';
    RLB = octstruct.plb(:)';
    RUB = octstruct.pub(:)';
    XScale = 2*octstruct.gamma(:)'.*(XScale./diff(InitRange,[],1));
    InitRange = [octstruct.plb(:)'; octstruct.pub(:)'];
    if ~isempty(x0)
        for i = 1:size(x0,1); x0(i,:) = optimct(x0(i,:),octstruct); end
    end
else
    octstruct = [];
end

% Loop over optimization epochs
for iEpoch = 1:maxEpochs
    % Initial points
    if isempty(x0) && strcmpi(MidpointStart,'on') % Reasonable starting point
        x0(1,:) = 0.5*(RLB+RUB);
    end
    
    if size(x0,1) < nStarts(iEpoch)
        index = size(x0,1)+1:nStarts(iEpoch);
        % x0(index,:) = bsxfun(@plus,bsxfun(@times,rand(length(index),nvars),RUB-RLB),RLB);
        temp = cache.x(~isnan(cache.fval),:);
        X = [output.x;output.startx;temp;x0(1:size(x0,1),:)];  
        x0(index,:) = randx(length(index),RLB,RUB,XScale/2,X);
    end
    x0 = x0(1:nStarts(iEpoch),:);

    maxeval = MaxFunEvals(min(iEpoch,end));
    maxiter = MaxIter(min(iEpoch,end));
    tolfun = TolFun(min(iEpoch,end));
    tolx = TolX(min(iEpoch,end));
    disp = Display{min(iEpoch,end)};
    method = Method{min(iEpoch,end)};
    
    if ~isempty(OutputFcn)
        optimValues.funccount = 0;
        optimValues.x = NaN(1,nvars);
        optimValues.fval = fvalmin;
        optimValues.epoch = iEpoch;
        optimValues.iteration = 0;
        optimValues.procedure = method;
        stop = OutputFcn(x0,optimValues,'iter');
    end
        
    for iStart = 1:nStarts(iEpoch)
        xstart = x0(iStart,:);
        
        if ~isempty(OutputFcn)
            optimValues.funccount = 0;
            optimValues.x = optimct(x,octstruct,1);
            optimValues.fval = fvalmin;
            optimValues.epoch = iEpoch;
            optimValues.iteration = iStart;
            optimValues.procedure = method;
            stop = OutputFcn(optimct(xstart,octstruct,1),optimValues,'iter');
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
                optoptions = optimset(optoptions,OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = fminsearchbnd(funfcn,xstart,LB,UB,optoptions);
                funcCount = output1.funcCount;
            case 'fmincon'
                assert(exist('patternsearch.m','file') == 2, ...
                    'FMINCON optimization method not installed (requires Optimization Toolbox).');
                optoptions = optimset('Display',disp,'TolX',tolx,'TolFun',tolfun,'MaxFunEvals',maxeval,'MaxIter',maxiter);
                optoptions = optimset(optoptions,OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = fmincon(funfcn,xstart,[],[],[],[],LB,UB,[],optoptions);
                funcCount = output1.funcCount;
            case 'patternsearch'
                assert(exist('patternsearch.m','file') == 2, ...
                    'PATTERNSEARCH optimization method not installed (requires Global Optimization Toolbox).');
                optoptions = psoptimset('Display',disp,'TolX',tolx,'TolFun',tolfun,'MaxFunEvals',maxeval,'MaxIter',maxiter);
                optoptions = psoptimset(optoptions,OptimOptions{min(iEpoch,end)});
                [x,fval,exitflag,output1] = patternsearch(funfcn,xstart,[],[],[],[],LB,UB,[],optoptions);
                funcCount = output1.funccount;
            case 'bps'
                assert(exist('bps.m','file') == 2, ...
                    'BPS optimization method not installed.');
                optoptions = OptimOptions{min(iEpoch,end)};
                optoptions.Display = disp;
                optoptions.TolX = tolx;
                optoptions.TolFun = tolfun;
                optoptions.MaxFunEvals = maxeval;
                optoptions.MaxIter = maxiter;
                PLB = options.InitRange(1,:);
                PUB = options.InitRange(2,:);
                [x,fval,exitflag,output1] = bps(funfcn,xstart,LB,UB,PLB,PUB,optoptions);
                funcCount = output1.FuncCount;                
        end

        % Store optimization results
        output = storeOptimization(output,...
            xstart,x,fval,funcCount,exitflag,output1);
        
        if ~isempty(OutputFcn)
            optimValues.funccount = funcCount;
            optimValues.x = optimct(x,octstruct,1);
            optimValues.fval = fval;
            optimValues.epoch = iEpoch;
            optimValues.iteration = iStart;
            optimValues.procedure = method;
            stop = OutputFcn(optimct(xstart,octstruct,1),optimValues,'iter');
        end
                
        if fval < fvalmin
            fvalmin = fval;
            xmin = x;
        end        
        
    end

    % Recompute x0 for next iteration
    if iEpoch < maxEpochs && size(output.x,1) > 1
        % Build optimization matrix combining optimal points and log likelihood
        xScale = XScale(:)'/10*sqrt(nvars);
        fval = output.fval(:);
        xx = [bsxfun(@rdivide,output.x,xScale), fval/FvalScale];
        
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

end

x = optimct(xmin,octstruct,1);
fval = fvalmin;
exitflag = 0;

if ~isempty(OutputFcn)
    optimValues.funccount = funcCount;
    optimValues.x = optimct(x,octstruct,1);
    optimValues.fval = fval;
    optimValues.epoch = iEpoch;
    optimValues.iteration = iStart;
    optimValues.procedure = [];
    stop = OutputFcn(optimValues.x,optimValues,'done');
end

% Convert function outputs to original coordinate system
if strcmp(RescaleVars,'on')
    if nargout > 3
        for i = 1:size(output.x,1)
            output.x(i,:) = optimct(output.x(i,:),octstruct,1);
        end
        for i = 1:size(output.startx,1)
            output.startx(i,:) = optimct(output.startx(i,:),octstruct,1);
        end
    end

    if nargout > 4
        for i = 1:size(Cache.x,1)
            Cache.x(i,:) = optimct(Cache.x(i,:),octstruct,1);
        end
    end
end
     
    %OBJECTIVEFUNC Objective function with caching
    function fval = objectivefunc(x)
                
        if strcmp(RescaleVars,'on'); x = optimct(x,octstruct,1); end
        
        if isempty(vararginarray)
            fval = fun(x);
        else
            fval = fun(x, vararginarray{min(iEpoch,end)}{:});
        end
        
        cache.x(cache.index,:) = x;
        cache.fval(cache.index) = fval;
        cache.index = max(mod(cache.index + 1,CacheSize+1), 1);
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
