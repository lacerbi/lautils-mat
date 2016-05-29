% TRIMNORMRND Random arrays from the normal distribution with trimmed tails.
%    R = NORMRND(MU,SIGMA,TRIM) returns an array of random numbers chosen from a
%    normal distribution with mean MU and standard deviation SIGMA.  The size
%    of R is the common size of MU and SIGMA if both are arrays.  If either
%    parameter is a scalar, the size of R is the size of the other
%    parameter. TRIM is the number of sigmas after which the tails are
%    trimmed. The data are not rescaled, but accumulate at the tails' end.
function r = trimnormrnd(mu, sigma, trim)
r = normrnd(mu, sigma);
r = max(r, mu - trim.*sigma);
r = min(r, mu + trim.*sigma);
end