%SLICESAMPLEBND_TEST Test script for SLICESAMPLEBND

%% Example 1: Sample from correlated bivariate normal distribution

display('Sample from a correlated bivariate normal pdf.');
display('(press a key to continue)');
pause;

rho = 0.9;
Mu = [2 2];
Sigma = [1 rho; rho 1];
Lambda = inv(Sigma);

% Correlated bivariate normal distribution (log pdf, without additive constants)
logf = @(x) -0.5*(x-Mu)*Lambda*(x-Mu)';

% Bounds (you must be sure that the probability mass outside these bounds
% is zero or negligible)
LB = [-20 -20];
UB = [20 20];

% Random starting point inside bounds (ideally, pick some starting point in 
% a region with high probability mass)
x0 = rand(1,2).*(UB-LB) + LB;

% We assume we have no idea about the typical length scales of the pdf we
% are sampling from, we just pick something large
widths = (UB-LB)/3;

% We want to get N samples in the end
N = 1e3;

% To stay on the safe side, we burn-in N/2 points (the default is 20% of N)
options.Burnin = N/2;

% We only record one sample every two (we will take 2N samples after 
% burn-in, and keep only half of them)
options.Thin = 2;

% Print major results to screen
options.Display = 'notify';

[samples,fvals,exitflag,output] = slicesamplebnd(logf,x0,N,widths,LB,UB,options);

names = {'x_1','x_2'};

% Plot samples
if exist('cornerplot.m','file')
    % Triangle plot using Will Adler's CORNERPLOT function
    halfwidth = 5*std(samples);
    PLB = mean(samples) - halfwidth;
    PUB = mean(samples) + halfwidth;
    cornerplot(samples,names,[],[PLB;PUB]);    
else
    figure;
    scatter(samples(:,1), samples(:,2));
    box off;
    set(gca,'TickDir','out');
    xlabel(names{1});
    ylabel(names{2});
    display('For better visualization, we recommended Will Adler''s <a href="http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation/content/kde2d.m">cornerplot.m</a> function.');
end
set(gcf,'Color','w');

display('Compare mean and covariance matrix of the target pdf to mean and cov of the samples:')

display(Mu)
display(Sigma)
display(mean(samples))
display(cov(samples))

%% Example 2: Sample from Rosenbrock (banana-shaped) function
display('Now sample from Rosenbrock (banana-shaped) function as log pdf. This will take a while.');
display('(press a key to continue)');
pause;

% More complex case: log pdf ~ negative Rosenbrock function
logf = @(x) -100*(x(2)-x(1).^2).^2-(1-x(1)).^2;

% The log pdf landscape is complex -- we need to sample more
N = 1e4;
options.Burnin = N/2;

[samples,fvals,exitflag,output] = slicesamplebnd(logf,x0,N,widths,LB,UB,options);

% Plot samples
if exist('cornerplot.m','file')
    % Triangle plot using Will Adler's CORNERPLOT function
    halfwidth = 5*std(samples);
    PLB = mean(samples) - halfwidth;
    PUB = mean(samples) + halfwidth;
    cornerplot(samples,names,[],[PLB;PUB]);    
else
    figure;
    scatter(samples(:,1), samples(:,2));
    box off;
    set(gca,'TickDir','out');
    xlabel(names{1});
    ylabel(names{2});
    display('For better visualization, we recommended Will Adler''s <a href="http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation/content/kde2d.m">cornerplot.m</a> function.');
end
set(gcf,'Color','w');

display('Additional information is contained in the OUTPUT structure:');
display(output);