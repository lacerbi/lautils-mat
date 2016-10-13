%TEMPERSLICE_TEST2

% func = @(x) log(p*exp(-0.5*(x-mu1)*L1*(x-mu1)')/sqrtdet1 + (1-p)*exp(-0.5*(x-mu2)*L2*(x-mu2)')/sqrtdet2);
func = @temperslice_test_func;

N = 1e4;
D = 5;
betas = [1 0.33 0.1 0.033 0.01];
% betas = [1];
M = size(betas,2);
widths = 1;
LB = -100*ones(1,D);
UB = 100*ones(1,D);
x0 = bsxfun(@plus, bsxfun(@times, rand(M,D), UB-LB), LB);
options.Display = 'final';
options.Debug = 1;
options.Thin = 1;
options.Swap = 23;
options.Tunnel = []; options.TunnelF = [];
tic
[samples,fvals,exitflag,output] = ...
    temperslice(func,x0,N,betas,widths,LB,UB,options);
toc

fprintf('Press a key to continue...\n');
pause

% Importance resampling
tsamples = samples(:,:,1); tfvals = fvals(:,1,1);
Px = exp(tfvals - betas(1)*(tfvals) - (1-betas(1))*max(tfvals));
Cx = cumsum(Px)/sum(Px);

r = rand(N,1);
[~,idx] = min(abs(bsxfun(@minus,r,Cx')),[],2);
tsamples = tsamples(idx,:);
tfvals = tfvals(idx);

betas = 1;
N = 5e3;
options.Thin = 5;
options.Tunnel = tsamples;
options.TunnelF = tfvals;
[samples,fvals,exitflag,output] = ...
    temperslice(func,x0,N,betas,widths,LB,UB,options);

M = size(betas,2);
figure;
for m = 1:M
    subplot(1,M,m);
    scatter(samples(:,1,m),samples(:,2,m));
    set(gca,'TickDir','out');
end
set(gcf,'Color','w');