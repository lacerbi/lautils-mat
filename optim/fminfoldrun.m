function [xs,fvals,ftests,exitflag,output] = fminfoldrun(optfun,maxfunevals,tolfun,varargin)
%FMINTASK Perform minimization of multiple related tasks.

if nargin < 2; maxfunevals = []; end
if nargin < 3 || isempty(tolfun); tolfun = 1e-3; end

% Get current log
foldlog = fminfoldfun();

if isempty(foldlog)
    error('FMINFOLDFUN not initialized. Ideally, you would run FMINFOLDRUN at the end of a full optimization with FMINFOLDFUN.');
end

nvars = size(foldlog.BestX,2);

% Max number of fcn evaluations per fold
if isempty(maxfunevals); maxfunevals = 2000*nvars; end

testflag = nargout > 4;         % Compute function value on test data

trainmask = foldlog.Mask;
testmask = foldlog.TestMask;

Nfolds = size(trainmask,1);
Ndata = size(trainmask,2);

bestfvals = foldlog.BestFval;       % Current best fold function values
bestfvals(foldlog.NewSearchFlag) = Inf;
exitflag = zeros(1,Nfolds);
output.funcCount = zeros(1,Nfolds);

fold2run = foldlog.NewSearchFlag;   % Folds that need to be run

iter = 1;

% Main loop
while any(fold2run)
    
    fprintf('FMINFOLD iteration %d.\n', iter);
    
    % Go through folds in random order
    for iFold = randperm(Nfolds)
        if ~fold2run(iFold); continue; end

        fprintf('FMINFOLD fold %d.\n', iFold);
        
        % Set current fold on wrapper function
        fminfoldfun([],[],[],[],iFold);

        % Start from best point for this fold
        x0 = foldlog.BestX(iFold,:);

        % Run optimizer
        if nargin > 3
            [x,fval,exitflag(iFold)] = optfun(@(x) fminfoldfun(x,varargin{:}),x0);
        else
            [x,fval,exitflag(iFold)] = optfun(@fminfoldfun,x0);            
        end

        % Get log of run
        foldlog = fminfoldfun();
        
        % Compute improvement for other folds (this one we just ran)
        bestfvals(iFold) = foldlog.BestFval(iFold);
        SearchImprovement = bestfvals - foldlog.BestFval;
        
        fold2run = SearchImprovement > tolfun;
        fold2run(foldlog.FuncCount > maxfunevals) = false;
    end
    
    iter = iter + 1;
end

output.funcCount = foldlog.FuncCount;

xs = foldlog.BestX;
fvals = foldlog.BestFval;
ftests = foldlog.TestFvalAtBestX;
