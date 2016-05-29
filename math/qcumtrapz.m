function z = qcumtrapz(y,dim)
%QCUMTRAPZ Quick cumulative trapezoidal numerical integration.
%   Z = QCUMTRAPZ(Y) computes an approximation of the cumulative
%   integral of Y via the trapezoidal method (with unit spacing).  To
%   compute the integral for spacing different from one, multiply Z by
%   the spacing increment.
%
%   For vectors, QCUMTRAPZ(Y) is a vector containing the cumulative
%   integral of Y. For matrices, QCUMTRAPZ(Y) is a matrix the same size as
%   X with the cumulative integral over each column. For N-D arrays,
%   QCUMTRAPZ(Y) works along the first non-singleton dimension.
%
%   Z = QCUMTRAPZ(Y,DIM) integrates along dimension DIM of Y. The length of 
%   X must be the same as size(Y,DIM).
%
%   QCUMTRAPZ is up to 3-4 times faster than CUMTRAPZ for large arrays.
%
%   See also CUMTRAPZ, QTRAPZ.

% By default integrate along the first non-singleton dimension
if nargin < 2; dim = find(size(y)>1,1); end    
    
% Compute dimensions of input matrix    
if isvector(y); n = 1; else n = ndims(y); end

y1 = [];
switch n
    case {1,2}      % 1-D or 2-D array
        switch dim
            case 1; y1 = y(1,:);
            case 2; y1 = y(:,1);
        end
        
    case 3      % 3-D array
        switch dim
            case 1; y1 = y(1,:,:);
            case 2; y1 = y(:,1,:);
            case 3; y1 = y(:,:,1);
        end           
                        
    case 4      % 4-D array
        switch dim
            case 1; y1 = y(1,:,:,:);
            case 2; y1 = y(:,1,:,:);
            case 3; y1 = y(:,:,1,:);
            case 4; y1 = y(:,:,:,1);
        end

    case 5      % 5-D array
        switch dim
            case 1; y1 = y(1,:,:,:,:);
            case 2; y1 = y(:,1,:,:,:);
            case 3; y1 = y(:,:,1,:,:);
            case 4; y1 = y(:,:,:,1,:);
            case 5; y1 = y(:,:,:,:,1);
        end
        
    otherwise
        z = cumtrapz(y,dim);
        return;
end

if isempty(y1)
    error('qcumtrapz:dimMismatch', 'DIM must specify one of the dimensions of Y.');
end

z = cumsum(y,dim) - 0.5*bsxfun(@plus, y, y1);

