function y = temperslice_test_func(x)

LogEps = 30;
persistent data;

p = 0.075;

if 0
if isempty(data) || size(x,2) ~= data.D
    data.D = size(x,2);
    data.mu1 = -5*ones(1,data.D);
    data.mu2 = 5*ones(1,data.D);

    data.s1 = 0.4*ones(1,data.D);
    data.L1 = diag(1./data.s1.^2);
    data.s2 = 0.3*ones(1,data.D);
    data.L2 = diag(1./data.s2.^2);
    data.sqrtdet1 = prod(data.s1);
    data.sqrtdet2 = prod(data.s2);
    
    z1 = (x-data.mu1)*data.L1*(x-data.mu1)';
    z2 = (x-data.mu2)*data.L2*(x-data.mu2)';

    if z1 - z2 > LogEps
        y = log(p) -0.5*z1 - log(data.sqrtdet1);    
    elseif z1 - z2 < -LogEps
        y = log(1-p) -0.5*z2 - log(data.sqrtdet2);
    else
        y = log(p*exp(-0.5*z1)/data.sqrtdet1 + (1-p)*exp(-0.5*z2)/data.sqrtdet2);
    end
    
end


else
    D = size(x,2);
    mu1 = -5*ones(1,D);
    mu2 = 5*ones(1,D);

    s1 = 0.4*ones(1,D);
    L1 = diag(1./s1.^2);
    s2 = 0.3*ones(1,D);
    L2 = diag(1./s2.^2);
    sqrtdet1 = prod(s1);
    sqrtdet2 = prod(s2);
    p = 0.075;
    
    z1 = log(p)   -0.5*(x-mu1)*L1*(x-mu1)' - log(sqrtdet1);
    z2 = log(1-p) -0.5*(x-mu2)*L2*(x-mu2)' - log(sqrtdet2);

    if z1 - z2 > LogEps
        y = z1;    
    elseif z1 - z2 < -LogEps
        y = z2;
    else
        y = log(exp(z1) + exp(z2));
    end
    
    
end


end