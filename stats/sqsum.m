function out = sqsum(X, DIM) 
% SQSUM  Return root of squared sum of data.
%
% SQSUM(X) returns the root of the squared sum of X
%
% SQSUM(X, DIM) root squared sum over dimension DIM.

if nargin==1; DIM = []; end

out = sqrt(sum(X.^2, DIM));