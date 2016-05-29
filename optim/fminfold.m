function [xs,fvals,exitflag,output,ftests] = fminfold(fun,x0,optfun,trainmask,options)
%FMINTASK Perform minimization of multiple related tasks.

if nargin < 5; options = []; end

nvars = numel(x0);

% Default options
defopts.TestMask = [];              % Mask of test trials for each fold
defopts.Permute = 0;                % Randomly permute order of folds
defopts.TolFun  = 1e-6;             % Tolerance of function value
defopts.MaxFunEvals = 1000*nvars;   % Max number of fcn evaluations per fold
defopts.DerivativeCheck = 'off';    % Check analytical vs numerical gradient

testflag = nargout > 4;         % Compute function value on test data

for f = fields(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

testmask = options.TestMask;

nFolds = size(trainmask,1);
nData = size(trainmask,2);
assert(isvector(x0) && all(isfinite(x0)) && isreal(x0), ...
    'The starting point X0 should be a real-valued scalar or vector.');
if testflag && isempty(testmask)
    testmask = (trainmask == 0);    % TEST is complement of TRAIN
end
if testflag
    assert(size(testmask,1) == nFolds && size(testmask,2) == nData, ...
        'TESTMASK must be a matrix of the same size as TRAINMASK.');
end
assert(all(sum(abs(trainmask),2) > 0), ...
    'At least one row of TRAINMASK contains all zeroes.');

xs = ones(nFolds,1)*x0(:)';
fvals = realmax*ones(1,nFolds);     % Current best value
exitflag = zeros(1,nFolds);
output.funcCount = zeros(1,nFolds);

fold2run = true(1,nFolds);          % Folds that need to be run

% Randomly permute fold order?
r = [];
if ~isempty(options.Permute) && any(options.Permute)
    if strcmpi(options.Permute,'allbutlast')
        r = [randperm(nFolds-1),nFolds];
    elseif strcmpi(options.Permute,'allbutfirst')
        r = [1,1+randperm(nFolds-1)];
    elseif any(strcmpi(options.Permute,{'no','none','off'}))
        
    else
        r = randperm(nFolds);
    end    
    [~,invr] = sort(r);
    trainmask = trainmask(r,:);
end

% No permutation
if isempty(r)
    r = 1:size(trainmask,1);
    invr = r;
end

% Main loop
iFold = nFolds;
while any(fold2run)
    iFold = mod(iFold,nFolds) + 1;
    if ~fold2run(iFold); continue; end
    
    % Stop optimization if current fold budget has been exhausted
    if output.funcCount(invr(iFold)) >= options.MaxFunEvals
        fold2run(iFold) = 0;
        continue;
    end
   
    % Initialize wrapper function
    fminfold_wrap([],fun,trainmask,iFold);
    
    % Start from best point for this fold
    x0 = xs(iFold,:);
    
    % Derivative check
    switch lower(options.DerivativeCheck)
        case {'on','yes'}
            % dy_num = fgrad(@fminfold_wrap,x0,'five-points');
            dy_num = gradest(@fminfold_wrap,x0);
            [~,dy_ana] = fminfold_wrap(x0);
                        
            err_rel = (dy_num - dy_ana)./dy_num;
            err_abs = dy_num - dy_ana;
            [~,idx_rel] = max(abs(err_rel));
            [~,idx_abs] = max(abs(err_abs));
            
            fprintf('Maximum error between numerical and analytical gradient:\n - Relative error: %10.3g (parameter #%d);\n - Absolute error: %10.3g (parameter #%d).\n', ...
                err_rel(idx_rel), idx_rel, err_abs(idx_abs), idx_abs);
            
        case {'off','no'}
            % Do nothing
            
        otherwise
            error('OPTIONS.DerivativeCheck should be ''on'' or ''off'' for gradient checking.');            
    end
    
    % Run optimizer    
    [x,fval,exitflag(iFold)] = optfun(@fminfold_wrap,x0);
   
    % Get log of run
    foldlog = fminfold_wrap();
    
    output.funcCount(invr(iFold)) = output.funcCount(invr(iFold)) + foldlog.FuncCount;
   
    % Find minimum function value per each fold
    [tempfvals,idx] = min(foldlog.Y,[],1);
    for i = 1:nFolds
        if tempfvals(i) < fvals(i) - options.TolFun
            fvals(i) = tempfvals(i);
            xs(i,:) = foldlog.X(idx(i),:);
            fold2run(i) = true;             % Fold has been updated
        end
    end
    
    % fvals
   
    % This fold has run successfully
    fold2run(iFold) = false;
end

xs = xs(invr,:);
fvals = fvals(invr);
exitflag = exitflag(invr);

% Compute function value on test sets
if testflag
    ftests = zeros(1,nFolds);
    for iFold = 1:nFolds
        if ~any(testmask(iFold,:)); continue; end        
        fminfold_wrap([],fun,testmask,iFold);
        ftests(iFold) = fminfold_wrap(xs(iFold,:));
    end
    output.funcCount = output.funcCount + 1;
end
