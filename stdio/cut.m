% CUT Cut a slice out of a matrix.
% 
% Parameters:
% X             (1,nx) Matrix
% SLICE         [min1, max1, ...]
% c             component being sliced (optional, normally 1)
%
% Returns:
% S             The sliced matrix
%
function S = cut(X, SLICE, c)
    if ~exist('c','var'); c = []; end
    if isempty(c); c = 1; end    
    nslices = floor(length(SLICE)/2);
    S = [];
    for i = 1:size(X, 1)
        take = 0;
        for j = 1:nslices            
            if (X(i, c) >= SLICE(j*2-1) && X(i, c) < SLICE(j*2)); take = 1; end
            if take; S = [S; X(i, :)]; end
        end
        
        
    end
       
end