%TEMPERSLICE_TEST

mu1 = [-5 -5];
mu2 = [5 5];

s1 = [1 0.3];
L1 = [1/s1(1)^2 0; 0 1/s1(2)^2];
s2 = [0.3 0.3];
L2 = [1/s2(1)^2 0; 0 1/s2(2)^2];
sqrtdet1 = prod(s1);
sqrtdet2 = prod(s2);
p = 0.075;

func = @(x) log(p*exp(-0.5*(x-mu1)*L1*(x-mu1)')/sqrtdet1 + (1-p)*exp(-0.5*(x-mu2)*L2*(x-mu2)')/sqrtdet2);

x0 = mu1;
N = 2e3;
%betas = [1 0.5 0.2 0.05 0.001];
% betas = [1 0.1];
betas = 0.01;
widths = 1;
LB = -100*[1 1];
UB = 100*[1 1];
options.Display = 'final';
options.Debug = 1;
options.Thin = 1;
options.Swap = 23;
options.Tunnel = []; options.TunnelF = [];
[samples,fvals,exitflag,output] = ...
    temperslice(func,x0,N,betas,widths,LB,UB,options);

% Importance resampling
tsamples = samples(:,:,1); tfvals = fvals(:,1,1);
Px = exp(tfvals - betas*(tfvals) - (1-betas)*max(tfvals));
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