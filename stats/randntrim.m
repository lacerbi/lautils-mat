%RANDNTRIM Normally distributed trimmed pseudorandom numbers.
%   R = RANDNTRIM(SIGMA, N) returns an N-by-N matrix containing pseudorandom 
%   values drawn from the standard normal distribution trimmed at SIGMA.  
%   RANDNTRIM(SIGMA, M,N) or RANDN(SIGMA, [M,N]) returns
%   an M-by-N matrix. RANDN(M,N,P,...) or RANDN([M,N,P,...]) returns an
%   M-by-N-by-P-by-... array. RANDN(SIGMA) returns a scalar.  
%   RANDNTRIM(SIGMA, SIZE(A)) returns an array the same size as A.
%
function r = randntrim(sigma, n, n2)

switch (nargin)
    case 1
        r = randn();
        while (abs(r) > sigma); r = randn(); end
        
    case 2        
        r = randn(n);
        t = abs(r) > sigma;

        while any(t(:))
            rtemp = randn(n);
            r(t) = rtemp(t);
            t = abs(r) > sigma;
        end
        
    case 3
        r = randn(n, n2);
        t = abs(r) > sigma;

        while any(t(:))
            rtemp = randn(n, n2);
            r(t) = rtemp(t);
            t = abs(r) > sigma;
        end
        
end

end


    
    
    
    