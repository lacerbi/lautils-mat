function a = intpower(a,b)
%INTPOWER  Array integer power.
%   C = INTPOWER(A,B) raises each element of A to the corresponding power
%   in B, where B is an integer scalar or matrix array.
%
%   INTPOWER can be considerably faster than POWER for a large matrix A and
%   for exponent B different than 2.
%
%   Example:
%         a = randn(1000,1000); b = 5; 
%         tic; a.^b; toc; tic; intpower(a,b); toc
%      compares the execution time of POWER vs INTPOWER on a large matrix.
%
%   See also POWER.

if numel(a) <= 100 || all(b == 2); a = a.^b; return; end

b = round(b);   % Round B to closest integer

if isscalar(b)
    switch b
        case 1; 
        case 2; a = a.*a;
        case 3; a = a.*a.*a;
        case 4; a = a.*a.*a.*a;
        case 5; a = a.*a.*a.*a.*a;
        case 6; a = a.*a.*a.*a.*a.*a;
        case 7; a = a.*a.*a.*a.*a.*a.*a;
        case 8; a = a.*a.*a.*a.*a.*a.*a.*a;
        case 9; a = a.*a.*a.*a.*a.*a.*a.*a.*a;
        case 10; a = (a.*a.*a.*a.*a).^2;
        case 11; a = (a.*a.*a.*a.*a).^2.*a;
        otherwise
            a = a.^b;
    end
end

