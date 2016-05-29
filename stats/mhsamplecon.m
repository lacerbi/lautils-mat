function [smpl,smpllogpdf,accept,x0opt,logpdfx0opt] = mhsamplecon(logpdf, nsamples, start, lb, ub, scale, cnst, datatype, opt)
% MHSAMPLECON Generate Markov chain using Metropolis-Hasting algorithm with
% constraints.
%
%   SMPL = MHSAMPLECON(LOGPDF,NSAMPLES,START,LB,UB,SCALE)
%   draws NSAMPLES random samples from a target stationary distribution 
%   LOGPDF using the Metropolis-Hasting algorithm. START is a row vector
%   containing the start value of the Markov Chain. NSAMPLES is an integer
%   specifying the number of samples to be generated. LOGPDF is a function
%   handle created using @. LOGPDF takes one argument as an input and this 
%   argument has the same type and size as START. LOGPDF returns a log 
%   density and it is not necessarily normalized. 
%   MHSAMPLECON uses an internal Gaussian proposal distribution which wraps 
%   around borders; it thus implements the Random Walk Metropolis-Hasting sampling. 
%   LB and UB are respectively lower and upper bounds on the variables.
%   SCALE is a row vector of the same size as START, specifying the size of 
%   the random jump for each dimension (standard deviation of the proposal
%   distribution). 
%   SMPL is a column vector or matrix containing the samples. 
%  
%   SMPL = MHSAMPLECON(LOGPDF,NSAMPLES,START,LB,UB,SCALE,CNST) 
%   checks that additional constraints are respected. CNST is a function
%   that takes as input a row vector of the type and size as START and
%   returns 1 if the constraints are satisfied, 0 if they are not.
%  
%   SMPL = MHSAMPLECON(LOGPDF,NSAMPLES,START,LB,UB,SCALE,CNST,DATATYPE) 
%   specifies the data type for the variables. DATATYPE is a vector of the
%   same length as START, and each entry codes the datatype for the relevant
%   data dimension. By default the data are real numbers (0), 1 represents
%   integer data. DATATYPE can be shorter than START, in which case all
%   the remaining unspecified elements are taken as zeroes.
%  
%   SMPL = MHSAMPLECON(..., OPT) 
%   if OPT is different than 0, execute an optimization of the scale factor
%   before performing the true MCMC. OPT in that case specifies the length
%   of the testing chains.
%
%   [SMPL,SMPLLOGPDF,ACCEPT] = MHSAMPLECON(...) also returns the samples
%   log probability and ACCEPT as the acceptance rate of the proposed 
%   distribution. If NSAMPLES is divisible by ten, ACCEPT is a ten elements 
%   row vector containing the acceptance rate at each stage of the sampling 
%   chain, divided in ten phases of equal length.
%
%   [SMPL,SMPLLOGPDF,ACCEPT,OPT,LOGPDFOPT] = MHSAMPLECON(...) also returns
%   the maximum likelihood sample OPT and its log density LOGPDFOPT.

% Constraint function handle
if ~exist('cnst', 'var'); cnst = []; end

% Data types
if ~exist('datatype', 'var'); datatype = []; end
datatype = [datatype zeros(1, length(start) - length(datatype))];

% Constraint function handle
if ~exist('opt', 'var'); opt = []; end
if isempty(opt); opt = 0; end

% Execute optimization of the scale factor
if opt > 0
    [smpl,smpllogpdf,accept,x0opt,logpdfx0opt] = mhsampleconopt(logpdf, nsamples, start, lb, ub, scale, cnst, datatype, opt);
    return;
end

%x0  is the place holder for the current value
x0 = start;
x0opt = start;
logpdfx0 = logpdf(x0);
logpdfx0opt = logpdfx0;
smpl = zeros(nsamples, length(x0));
smpllogpdf = zeros(nsamples, 1);

% The sampling chain is divided in ten stages
if mod(nsamples, 10) > 0
    display('Warning: NSAMPLES is not divisible by ten. ACCEPT ratio will be a scalar.');
    accept = 0;
    phaselength = nsamples;
else
    accept = zeros(1, 10);
    phaselength = nsamples/10;
    display('Sampling phase 1 out of 10.');
end
phasestart = 1;
phaseindex = 1;

% Metropolis-Hasting Algorithm.
U = log(rand(1,nsamples));
for i = 1:nsamples
    % Sample from proposal distribution
    y = proprnd(x0);    
    % Accept or reject the proposal.
    rho = logpdf(y)-logpdfx0;    
    acc = (U(i)<= min(rho, 0));
    if acc
        x0 = y;
        logpdfx0 = logpdf(x0);
        if logpdfx0 > logpdfx0opt
            logpdfx0opt = logpdfx0;
            x0opt = x0;
        end
        accept(phaseindex) = accept(phaseindex)+1;
    end
    smpl(i, :) = x0;
    smpllogpdf(i) = logpdfx0;
    phasestart = phasestart + 1;
    if phasestart > phaselength
        phaseindex = phaseindex + 1;
        phasestart = 1;
        if phaseindex <= 10; display(['Sampling phase ' num2str(phaseindex) ' out of 10.']); end
    end
end

% Calculate acceptance rate
accept = accept/phaselength;

    % Proposal distribution (Gaussian with wrapping)
    function y = proprnd(x)
        while 1
            y = normrnd(x, scale);
            for j = 1:length(y)
               while 1
                  if y(j) > ub(j)
                      y(j) = 2*ub(j) - y(j); 
                  elseif y(j) < lb(j)
                      y(j) = 2*lb(j) - y(j); 
                  else break;
                  end
               end
            end
            y(datatype == 1) = round(y(datatype == 1));
            % If there are no constraints, or they are satisfied, quit
            if isempty(cnst) || cnst(y); break; end            
        end
    end
end


% MHSAMPLECONOPT Generate Markov chain by finding first optimal scale
% factor.
function [smpl,smpllogpdf,accept,x0opt,logpdfx0opt] = mhsampleconopt(logpdf, nsamples, start, lb, ub, scale, cnst, datatype, opt)
    direction = 0;
    i = 1;
    alpha(1) = 1;
    okrange = [0.23 0.27];
    
    % Convergence parameters (converging from above or below).
    % If during one jump the sign of direction switches, reduces the 
    % jumping factor k.
    k = 0.2;

    while 1
        display(['Testing sampling scale factor: ' num2str(alpha(i))]);

        [smpl,smpllogpdf,accept] = mhsamplecon(logpdf, opt, start, lb, ub, scale*alpha(i), cnst, datatype, 0);            
        accept(i) = mean(accept(6:10));

        if accept(i) > okrange(1) && accept(i) < okrange(2)
            break;
        elseif accept(i) >= okrange(2)
            if direction == -1; k = k*0.75; end
            alpha(i + 1) = alpha(i)*(1 + k);
            direction = 1;
        else
            if direction == 1; k = k*0.75; end
            alpha(i + 1) = alpha(i)/(1+k); 
            direction = -1;
        end
        display(['Current acceptance rate: ' num2str(accept(i))]);
        i = i + 1;
    end
    display(['Final acceptance rate: ' num2str(accept(i))]);

    display('Starting proper sampling...');
    [smpl,smpllogpdf,accept,x0opt,logpdfx0opt] = mhsamplecon(logpdf, nsamples, start, lb, ub, scale*alpha(i), cnst, datatype, 0);            
    display(['Acceptance rate: ' num2str(mean(accept(2:10)))]);

end
