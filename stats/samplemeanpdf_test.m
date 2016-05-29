function samplemeanpdf_test(N,fun,genfun,titlestring)

if nargin < 1 || isempty(N); N = 5e5; end
if nargin < 2; fun = []; end
if nargin < 4; titlestring = 'Custom pdf'; end

%% Run standard battery test
if isempty(fun)

    % Test normal distribution
    %fun = @(x) normpdf(x,0,1);
    %genfun = @(n,N) randn(n,N);
    %titlestring = 'Normal distribution';
    %test_pdf(N,fun,genfun,titlestring);

    % Test uniform distribution
    fun = @(x) (x >= -3 & x <= 3)./trapz(x >= -3 & x <= 3);
    genfun = @(n,N) (rand(n,N)-0.5)*6;
    titlestring = 'Uniform distribution';
    test_pdf(N,fun,genfun,titlestring);
    
    % Test bimodal distribution
    fun = @(x) 0.5*(normpdf(x,-3,1)+normpdf(x,3,1));
    genfun = @(n,N) 3*(2*(rand(n,N) < 0.5)-1)+randn(n,N);
    titlestring = 'Bimodal distribution';
    test_pdf(N,fun,genfun,titlestring);

    % Test skewed distribution (exponential)
    fun = @(x) 1/2*exp(-(x-7.5)/2).*(x > -7.5);
    genfun = @(n,N) -(log(1-rand(n,N)))*2-7.5;
    titlestring = 'Exponential distribution';
    test_pdf(N,fun,genfun,titlestring);
    
    % Test vectorization (fixed n, multiple pdfs)
    x = linspace(-10,10,1024)';
    mu = linspace(-7,7,1000);
    y = bsxfun(@rdivide, bsxfun(@ge, x, -2 + mu) & bsxfun(@le, x, 2 + mu), qtrapz(bsxfun(@ge, x, -2 + mu) & bsxfun(@le, x, 2 + mu),1));
    titlestring = 'Multiple uniform distributions';
    % y = 0.5*bsxfun_normpdf(x,mu-3,1) + 0.5*bsxfun_normpdf(x,mu+3,1);    
    test_vectorization(x,y,titlestring);
    
else
%% Test user-provided pdf

    test_pdf(N,fun,genfun,titlestring);

end

end % SAMPLEMEANPDF_TEST

%-------------------------------------------------------------------------%
%% Private functions

%TEST_PDF Test a provided pdf
function test_pdf(N,fun,genfun,titlestring)

% Initialize variables
x = linspace(-10, 10, 1024)';
dx = x(2)-x(1);
edges = linspace(min(x),max(x),200);
binx = edges(1:end-1) + 0.5*diff(edges);
z = zeros(size(x,1),9); h = zeros(length(binx),9); clt = zeros(size(x,1),9);
nmax = 12;

y = fun(x);
nrange = [1 2 3 4, 5 6 9 12, 20 30 50 100];

display([titlestring ': Compute sampling pdf of the mean with samplemeanpdf...']);
tic;
for n = 1:nmax; z(:,n) = samplemeanpdf(y,nrange(n)); end
toc;
z = bsxfun(@rdivide, z, qtrapz(z,1)*dx);                    % Normalize

% Monte Carlo approximation
display([titlestring ': Compute sampling pdf of the mean via Monte Carlo sampling...']);
tic;
for n = 1:nmax
    m = mean(genfun(nrange(n),N),1);
    h(:,n) = histcounts(m, edges);
end
toc;
h = bsxfun(@rdivide, h, qtrapz(h,1)*(binx(2)-binx(1)));     % Normalize

% CLT normal approximation
for n = 1:nmax; clt(:,n) = samplemeanpdf(y,nrange(n),0); end
clt = bsxfun(@rdivide, clt, qtrapz(clt,1)*dx);     % Normalize

% Plot figure
figure;
set(gcf,'Color','w');
for n = 1:nmax
    subplot(3,4,n);
    plot(x,z(:,n),'k','LineWidth',2);
    hold on;
    plot(binx,h(:,n),'r:','LineWidth',2);
    plot(x,clt(:,n),'b-.','LineWidth',1);
    if n == 1; title(titlestring); end
    if n == 4
        hl = legend('samplemeanpdf', 'Monte Carlo','Normal approximation');
        set(hl,'Location','NorthEast');
    end        
    if any(n == [9 10 11]); xlabel('x'); end
    if n == 5; ylabel('sample-mean pdf'); end
    axis([min(x) max(x) 0 1]);
    set(gca,'TickDir','out');
    box off;
    ht = text(6.5,0.25,['n = ' num2str(nrange(n))]);
end

end

%-------------------------------------------------------------------------%
%TEST_VECTORIZATION Test vectorized pdf and non-integer N
function test_vectorization(x,y,titlestring)

nplot = 4;
nrange = linspace(1,3.75,12);
nmax = length(nrange);
display('Test vectorization and non-integer n.');
figure;
set(gcf,'Color','w');
for n = 1:nmax
    tic;
    display([titlestring ': compute ' num2str(size(y,2)) ' distributions for n = ' num2str(nrange(n),'%.3g') ' (plot only ' num2str(nplot) ').']);
    z = samplemeanpdf(x,y,nrange(n));
    toc;
    tol = 1e-6;
    z(z < tol & circshift(z,1) < tol & circshift(z,-1) < tol) = NaN;
    subplot(3,4,n);
    index = round(linspace(1,size(z,2),nplot));    
    plot(x,z(:,index),'k','LineWidth',2);
    set(gca,'TickDir','out');
    box off;    
    if n == 1; title(titlestring); end
    if n >= 9; xlabel('x'); end
    if any(n == [1 5 9]) ; ylabel('sample-mean pdf'); end
    axis([min(x) max(x) 0 1]);
    ht = text(6.5,0.85,['n = ' num2str(nrange(n),'%.3g')]);
end

end