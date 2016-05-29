function z = bsxfunandsum(varargin)
%BSXFUNANDSUM  Element-by-element matrix operation followed by sum.
%   Z = BSXFUNANDSUM(FUN,A,B) applies the element-by-element binary 
%   operation specified by the function handle FUN to arrays A and B, with 
%   singleton expansion enabled, and then returns the sum of the 
%   intermediate array along the first non-singleton dimension.
%
%   BSXFUNANSUM prevents out-of-memory errors and severe slow-downs caused 
%   by virtual memory access by splitting the computation in multiple 
%   intermediate steps with smaller arrays.
%
%   Z = BSXFUNANDSUM(FUN[1],FUN[2],...,FUN[k],A[1],A[2],...A[k],A[k+1]) 
%   applies the element-by-element binary operation specified by the 
%   function handle FUN[1] to arrays A[1] and A[2], yielding intermediate
%   array B[1]. Function FUN[2] is then applied to arrays B[1] and A[3], 
%   yielding array B[2]. FUN[3] is applied to arrays B[2] and A[4}, and so 
%   on. The final intermediate array B[k-1] is summed along the first 
%   non-singleton dimension and returned.
%
%   Z = BSXFUNANDSUM(...,DIM) sums the final intermediate array along the 
%   dimension DIM.
%
%   Z = BSXFUNANDSUM(...,DIM,METHOD) specifies alternate methods for
%   compressing the final array along dimension DIM. The default is sum of 
%   elements. Use an empty array [] to specify the default. Available 
%   methods are:
%
%     'prod'     - product of elements
%     'qtrapz'   - fast trapezoidal integration (with unit spacing)
%     'sum'      - sum of elements
%     'trapz'    - trapezoidal integration (with unit spacing), deprecated
%
%   Z = BSXFUNANDSUM(...,DIM,METHOD,TOLMEM) splits the computation in
%   steps that require at most approximately TOLMEM megabytes of memory. 
%   Default TOLMEM is 2048. TOLMEM should be (much) lower than the 
%   available physical memory (obtainable via MEMORY). It is recommended
%   to use TOLMEM=1024 for systems with less than 8GB of RAM.
%
%   Example:
%   A = ones(1000,2000,1); B = ones(1000,1,3000);
%   tic; Z = bsxfunandsum(@times,A,B,3); toc
%   
%   returns (after a while) a 1000 by 2000 array Z whose elements are all 
%   3000. The intermediate array would have required more than 44GB, as 
%   the following command reveals:  tic; Z = sum(bsxfun(@times,A,B),3); toc
%
%   See also BSXFUN, MEMORY, PROD, QTRAPZ, SUM.

% Get number of matrix inputs (equal to # of function handles plus one)
nmat = find(cellfun(@isnumeric,varargin),1);    
for i = 1:nmat-1; funs{i} = varargin{i}; end
for i = 1:nmat; mat{i} = varargin{nmat-1+i}; end

if nargin < nmat*2; dim = []; else dim = varargin{nmat*2}; end
if nargin < nmat*2+1 || isempty(varargin{nmat*2+1}); method = 'sum'; else method = varargin{nmat*2+1}; end
if nargin < nmat*2+2 || isempty(varargin{nmat*2+2}); TolMem = 2048; else TolMem = varargin{nmat*2+2}; end

% Matrix dimensions
nd = cellfun(@ndims, mat);
nm = max(nd);               % Intermediate matrix M ndims

% Compute matrix sizes
sizemats = ones(nmat,nm);
for i = 1:nmat; sizemats(i,1:nd(i)) = size(mat{i}); end    
sizem = max(sizemats,[],1);   % Intermediate matrix M size

% By default integrate along the first non-singleton dimension of M
if isempty(dim); dim = find(sizem>1,1); end

% Size of matrix M in MB (assuming 8 bytes per element)
mbytes = prod(sizem)/131072;

if mbytes < TolMem          % Intermediate matrix size smaller than TolMem
    y = bsxfun(funs{1},mat{1},mat{2});
    for i = 3:nmat; y = bsxfun(funs{i-1},y,mat{i}); end
    op = str2func(method);
    z = op(y,dim);
else                        % Intermediate matrix size larger than TolMem
    % Divide the computation in this # of chunks
    nchunks = min(1 + floor(mbytes/TolMem), sizem(dim));

    % Compute indices
    for i = 1:nmat
        for iDim = 1:nm; index{i}{iDim} = 1:sizemats(i,iDim); end
        lindex{i} = floor(linspace(1,sizemats(i,dim)+1,nchunks+1));
    end

    % Output matrix
    zdims = max(sizemats,[],1);
    zdims(dim) = 1;

    switch lower(method)
        case 'prod'
            z = ones(zdims);
            for n = 1:nchunks
                indexnow = index;
                for i = 1:nmat; indexnow{i}{dim} = lindex{i}(n):max(1,lindex{i}(n+1)-1); end                    
                y = bsxfun(funs{1},mat{1}(indexnow{1}{:}),mat{2}(indexnow{2}{:}));
                for i = 3:nmat; y = bsxfun(funs{i-1},y,mat{i}(indexnow{i}{:})); end                    
                z = z .* prod(y,dim);
            end
        case 'sum'
            z = zeros(zdims);
            for n = 1:nchunks
                indexnow = index;
                for i = 1:nmat; indexnow{i}{dim} = lindex{i}(n):max(1,lindex{i}(n+1)-1); end                    
                y = bsxfun(funs{1},mat{1}(indexnow{1}{:}),mat{2}(indexnow{2}{:}));
                for i = 3:nmat; y = bsxfun(funs{i-1},y,mat{i}(indexnow{i}{:})); end                    
                z = z + sum(y,dim);
            end
        case {'qtrapz','trapz'}
            z = zeros(zdims);
            op = str2func(method);
            for n = 1:nchunks
                indexnow = index;
                for i = 1:nmat; indexnow{i}{dim} = lindex{i}(n):min(sizemats(i,dim),lindex{i}(n+1)); end                    
                y = bsxfun(funs{1},mat{1}(indexnow{1}{:}),mat{2}(indexnow{2}{:}));
                for i = 3:nmat; y = bsxfun(funs{i-1},y,mat{i}(indexnow{i}{:})); end                    
                z = z + op(y,dim);
            end
        otherwise
            error('bsxfunandsum:unknownMethod',...
                'Supported methods are PROD, QTRAPZ, SUM, TRAPZ.');                
    end
end