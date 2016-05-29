function theta = qrandvm(mu,kappa,n,kappamin,kappamax)
%QRANDVM Quick Von Mises distributed pseudorandom numbers.
%
%   THETA = QRANDVM(MU,KAPPA,N) returns a column vector of N pseudorandom 
%   angles drawn from the Von Mises distribution with circular mean MU
%   and concentration parameter KAPPA. MU and KAPPA can be scalars or
%   matrices with distinct parameter values.
%   N can be a scalar or a vector. If N is a vector (e.g., [2 10 8]), the 
%   function creates a matrix output with the respective dimensionality. 
%   If N is not specified, the default value is 1 if both MU and KAPPA 
%   are scalars, otherwise the size of MU and/or KAPPA.
%   The returned angles are expressed in radians in the range [-pi, pi].
%
%   THETA = QRANDVM(MU,KAPPA,N,KAPPAMIN,KAPPAMAX) uses the uniform
%   approximation for KAPPA < KAPPAMIN and the normal approximation for
%   KAPPA > KAPPAMAX. Default values are KAPPAMIN=1e-6 and KAPPAMAX=1e3.
%
%   For maximal speed, take advantage of the full vectorization of QRANDVM 
%   by generating all required Von Mises random numbers in a single call.
%
%   Improvements over CIRC_VMRND from the Circular Statistics Toolbox:
%     * Fully vectorized (1-2 orders of magnitude faster with large N).
%       As a comparison, try: (run it several times)
%       tic;circ_vmrnd(0,1,1e4);t1=toc;tic;qrandvm(0,1,1e4);t2=toc;t1/t2
%     * Supports matrix values of MU and KAPPA.
%     * N can be a matrix with any number of dimensions.
%     * THETA is always in the range [-pi,pi], even for small KAPPA.
%     * Uses normal approximation (faster) for large values of KAPPA.
%
%   See also CIRC_VMRND.

%   References:
%     Statistical analysis of circular data, Fisher, sec. 3.3.6, p. 49.
%   
% Author:   Luigi Acerbi
% Date:     May 4, 2016
% Email:    <luigi.acerbi@gmail.com>
%
% Notes: Inspired by 'circ_vmrnd.m' from the Circular Statistics Toolbox
% for Matlab by Philipp Berens and Marc J. Velasco, 2009.

% Default parameters
if nargin < 1 || isempty(mu) 
    mu = 0;
end
if nargin < 2 || isempty(kappa)
    kappa = 1;
end
if nargin < 3 || isempty(n)
    if isscalar(mu) && isscalar(kappa)
        n = 1;
    elseif isscalar(mu)
        n = size(kappa);
    else
        n = size(mu);
    end
end
if nargin < 4 || isempty(kappamin)
    kappamin = 1e-6;
end
if nargin < 5 || isempty(kappamax)
    kappamax = 1e3;
end

if isscalar(n)
    n1 = find(size(mu) > 1,1);
    n2 = find(size(kappa) > 1,1);
    d = min([n1(:);n2(:)]);
    if isempty(d); d = 1; end
    dmax = max(ndims(mu),ndims(kappa));
    m = ones(1,dmax);
    m(d) = n;
else
    m = n;
    n = prod(n);
end

if ~isscalar(mu) && ~all(size(mu) == m | size(mu) == 1)
    error('MU must be a scalar or have the same size as KAPPA and the size specified by N.');
end

if ~isscalar(kappa) && ~all(size(kappa) == m | size(kappa) == 1)
    error('KAPPA must be a scalar or have the same size as MU and the size specified by N.');
end

mu = mu(:);
kappa = kappa(:);

% Small and large KAPPA cases
smallkappa = kappa < kappamin;
largekappa = kappa > kappamax;

% If KAPPA is scalar and small/large, use approximations and quit
if isscalar(kappa)
    if smallkappa
        if numel(m) == 1
            theta = 2*pi*rand(n,1)-pi;
        else
            theta = 2*pi*rand(m)-pi;
        end
        return;
    elseif largekappa
        sigma = 1./sqrt(kappa);
        theta = bsxfun(@plus, mu, sigma*randn(m));
        theta = mod(theta + pi, 2*pi) - pi;
        return;
    end        
end

% Other cases
medkappa = ~smallkappa & ~largekappa;

a = 1 + sqrt((1+4*kappa.^2));
b = (a - sqrt(2*a))./(2*kappa);
r = (1 + b.^2)./(2*b);

if n > 10   % Large N, use fully vectorized code
    
    % Initialize loop
    todo = true(1,n);
    todo(smallkappa | largekappa) = false;   % Treat small KAPPA separately
    f = zeros(n,1);

    % Loop until there are no more angles to draw
    while any(todo)
        u = rand(sum(todo),2);

        z = cos(pi*u(:,1));
        if numel(kappa) > 1
            f(todo) = (1+r(todo).*z)./(r(todo)+z);
            c = kappa(todo).*(r(todo)-f(todo));
        else
            f(todo) = (1+r.*z)./(r+z);
            c = kappa*(r-f(todo));
        end
        todo(todo) = (u(:,2) >= (c .* (2-c))) & (log(c)-log(u(:,2)) + 1 - c < 0);    
    end

    % If there is any small or large KAPPA, treat them now
    if any(smallkappa | largekappa)
        theta = zeros(n,1);
        
        % Uniform approximation for small kappa
        theta(smallkappa) = 2*pi*rand(sum(smallkappa),1)-pi;
        
        % Normal approximation for large kappa
        sigma = 1./sqrt(kappa(largekappa));
        if ~isscalar(mu); thetalarge = mu(largekappa); else thetalarge = mu; end
        temp = bsxfun(@plus, thetalarge, sigma.*randn(sum(largekappa),1));
        theta(largekappa) = mod(temp + pi, 2*pi) - pi;
        
        % Other cases
        if ~isscalar(mu); mu = mu(medkappa); end
        theta(medkappa) = mu + sign(rand(sum(medkappa),1) - 0.5) .* acos(f(medkappa));
        theta(medkappa) = angle(exp(1i*theta(medkappa)))';
    else
        theta = mu + sign(rand(n,1) - 0.5) .* acos(f);
        theta = angle(exp(1i*theta))';
    end
    
else    % Small N, use alternative implementation (faster for small N)
    
    if isscalar(mu); mu = mu*ones(n,1); end
    if isscalar(kappa)
        medkappa = true(n,1);
        kappa = kappa*ones(n,1);
        r = r*ones(n,1);
    end
    
    theta = zeros(n,1);    
    for j = find(medkappa)'
        while true
            u = rand(3,1);
            z = cos(pi*u(1));
            f = (1+r(j)*z)/(r(j)+z);
            c = kappa(j)*(r(j)-f);

            if u(2) < c * (2-c) || ~(log(c)-log(u(2)) + 1 -c < 0)
                break
            end
        end

        theta(j) = mu(j) +  sign(u(3) - 0.5) * acos(f);
        theta(j) = angle(exp(1i*theta(j)));
    end
    
    % If there is any small KAPPA, treat them now
    if any(smallkappa)
        theta(smallkappa) = 2*pi*rand(1,sum(smallkappa))-pi;
    end
    % If there is any large KAPPA, treat them now
    if any(largekappa)
        sigma = 1./sqrt(kappa(largekappa));
        if ~isscalar(mu); mu = mu(largekappa); end
        temp = bsxfun(@plus, mu, sigma.*randn(sum(largekappa),1));
        theta(largekappa) = mod(temp + pi, 2*pi) - pi;        
    end
end

% Reshape THETA to the desired matrix shape
if numel(m) > 1
    theta = reshape(theta,m);
end