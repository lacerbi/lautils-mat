function simulateentropy(K)

lnalpha_range = linspace(-30,30,2^10);
dlnalpha = lnalpha_range(2) - lnalpha_range(1);
alpha_range = exp(lnalpha_range);

meanh = psi(K*alpha_range+1) - psi(alpha_range+1);

prior_nsb = (K*psi(1,K*alpha_range+1) - psi(1,alpha_range+1)) .* alpha_range;
prior_nsb = prior_nsb./(qtrapz(prior_nsb)*dlnalpha);

plot(lnalpha_range, log(prior_nsb),'r','LineWidth', 1); hold on;

varh = (alpha_range + 1)./(K*alpha_range + 1) .* psi(1, alpha_range+1) - psi(1, K*alpha_range+1);

gprime = 1./sqrt(varh);
g = qcumtrapz(gprime.*alpha_range)*dlnalpha;
dg = g(2)-g(1);

% p(alpha) dalpha = p(beta)  = p(H) dH/dbeta dbeta/dalpha 

%plot(lnalpha_range, meanh.*alpha_range, 'k'); hold on;
%plot(lnalpha_range, std(varh).*alpha_range, 'b'); hold on;
%pause

prior_nsb2 = (K*psi(1,K*alpha_range+1) - psi(1,alpha_range+1)) ./ gprime .* alpha_range;
prior_nsb2 = prior_nsb2./(qtrapz(prior_nsb2)*dlnalpha);

plot(lnalpha_range, log(prior_nsb2),'k','LineWidth', 1);

box off;
set(gca,'TickDir','out');
xlabel('$\log \alpha$', 'Interpreter', 'LaTeX');
ylabel('$p(\log \alpha)$', 'Interpreter', 'LaTeX');
set(gcf,'Color','w');

figure(2);
subplot(1,2,1);
sample_entropy(K,alpha_range,prior_nsb);

subplot(1,2,2);
sample_entropy(K,alpha_range,prior_nsb2);

end

function sample_entropy(K,alpha_range,prior)

prior = prior/sum(prior);

[~,idx] = max(mnrnd_private(prior,1,1e5),[],2);
p = gamrnd(alpha_range(idx)'*ones(1,K),1);
p = bsxfun(@rdivide, p, sum(p,2));
H = -sum(p.*log(p),2);
hist(H,50);


end

%--------------------------------------------------------------------------
function r = mnrnd_private(p,n,m)
%MNRND_PRIVATE Random sample from the multinomial distribution.

if nargin < 2 || isempty(n); n = 1; end
if nargin < 3 || isempty(m); m = size(p,1); end
cdf = cumsum(p,2);
rr = rand(m,1,n);
samp_k = bsxfun(@gt, cdf, rr) & bsxfun(@le, [zeros(size(cdf,1),1),cdf(:,1:end-1)], rr);
r = sum(samp_k,3);

end