function out = sqmean(X, DIM) 
% SQMEAN  Return root of squared mean of data.
%
% SQMEAN(X) returns the root of the squared mean of X
%
% SQMEAN(X, DIM) root squared mean over dimension DIM.

if nargin==1; DIM = []; end

if isempty(DIM)
    out = sqrt(nanmean(X.^2));
else
    out = sqrt(nanmean(X.^2, DIM));
end