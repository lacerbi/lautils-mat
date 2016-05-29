function varargout = coordtransf(varargin)
%COORDTRANF Coordinate transform
%
% 
NumEps = 1e-10; % Accepted numerical error

% MaxPrecision = 17; % Maximum precision for a double

%% Apply coordinate transform
if nargin < 4
    
    coordstruct = varargin{2};    

    if isempty(coordstruct)
        varargout{1} = varargin{1}; % Return untransformed input
    else
        if (nargin < 3); direction = 'd'; else direction = varargin{3}(1); end
        
        switch lower(direction(1))
            case 'd'
                x = varargin{1};
                if ischar(coordstruct.g); g = str2func(coordstruct.g); else g = coordstruct.g; end
                y = g(x(:));
                y = min(max(y,coordstruct.lb),coordstruct.ub);    % Force to stay within bounds
                varargout{1} = reshape(y,size(x));
            case 'i'
                y = varargin{1};
                if ischar(coordstruct.ginv); ginv = str2func(coordstruct.ginv); else ginv = coordstruct.ginv; end
                x = ginv(y(:));
                x = min(max(x,coordstruct.oldbounds.lb),coordstruct.oldbounds.ub);    % Force to stay within bounds
                varargout{1} = reshape(x,size(y));
            otherwise
                error(['Unknown direction ' direction ' for coordinate transform.' ...
                    'Use DIRECTION ''d'' for a direct transform and ''i'' for the inverse transform.']);
        end
    end
    
%% Initialize coordinate struct
elseif nargin >= 4
    
    nvars = varargin{1};
    lb = varargin{2}(:);
    ub = varargin{3}(:);
    plb = varargin{4}(:);
    pub = varargin{5}(:);
    if nargin < 6; coordstruct.logtransf = []; else coordstruct.logtransf = varargin{6}(:); end

    % Convert scalar inputs to column vectors
    if isscalar(lb); lb = lb*ones(nvars,1); end
    if isscalar(ub); ub = ub*ones(nvars,1); end
    if isscalar(plb); plb = plb*ones(nvars,1); end
    if isscalar(pub); pub = pub*ones(nvars,1); end

    %assert(all(isnumeric(scale) & scale > 0 & isreal(scale)), ...
    %    'All SCALE parameters need to be positive real numbers.');

    plb = max(plb,lb);
    pub = min(pub,ub);
    
    assert(~any(isinf(([plb; pub]))), ...
        'Plausible interval bounds PLB and PUB need to be finite.');
    
    %if any(scale > (ub-lb))
    %    scale = min(scale,ub-lb);
    %    warning('Detected SCALE parameter(s) larger than the allowed coordinate bounds; reduced to the maximum coordinate range.');
    %end
        
    if isempty(coordstruct.logtransf)
        % A variable is converted to log scale if all bounds are 
        % non-negative and the plausible range spans at least one orders of magnitude
        coordstruct.logtransf = all([lb, ub, plb, pub] >= 0, 2) & (pub./plb >= 10);    
    elseif isscalar(coordstruct.logtransf)
        coordstruct.logtransf = logical(coordstruct.logtransf*ones(nvars,1));
    end
    
    % Transform bounds to log coordinates
    coordstruct.oldbounds.lb = lb;
    coordstruct.oldbounds.ub = ub;
    coordstruct.oldbounds.plb = plb;
    coordstruct.oldbounds.pub = pub;    
    coordstruct.lb = lb; coordstruct.ub = ub;
    coordstruct.plb = plb; coordstruct.pub = pub;
    coordstruct.lb(coordstruct.logtransf) = log(coordstruct.oldbounds.lb(coordstruct.logtransf));
    coordstruct.ub(coordstruct.logtransf) = log(coordstruct.oldbounds.ub(coordstruct.logtransf));
    coordstruct.plb(coordstruct.logtransf) = log(coordstruct.oldbounds.plb(coordstruct.logtransf));
    coordstruct.pub(coordstruct.logtransf) = log(coordstruct.oldbounds.pub(coordstruct.logtransf));

    coordstruct.mu = 0.5*(coordstruct.plb + coordstruct.pub);
    coordstruct.gamma = 0.5*(coordstruct.pub - coordstruct.plb);
    
    z = ['( (x - ' vec2str(coordstruct.mu) ') ./ ' vec2str(coordstruct.gamma) ' )'];
    zlog = ['( (log(abs(x) + (x == 0)) - ' vec2str(coordstruct.mu) ') ./ ' vec2str(coordstruct.gamma) ' )'];
    
    switch sum(coordstruct.logtransf)
        case 0
            coordstruct.g = ['@(x) ' z ];
            coordstruct.ginv = ['@(y) ' vec2str(coordstruct.gamma) ' .* y + ' vec2str(coordstruct.mu) ];        
            %coordstruct.g = ['@(x) linlog(' z ')'];
            %coordstruct.ginv = ['@(y) ' vec2str(coordstruct.gamma) ' .* linexp(y) + ' vec2str(coordstruct.mu) ];        
            
        case nvars
            %coordstruct.g = ['@(x) linlog(' zlog ')'];
            %coordstruct.ginv = ['@(y) exp(' vec2str(coordstruct.gamma) ' .* linexp(y) + ' vec2str(coordstruct.mu) ')'];        
            coordstruct.g = ['@(x) ' zlog ];
            coordstruct.ginv = ['@(y) exp(' vec2str(coordstruct.gamma) ' .* y + ' vec2str(coordstruct.mu) ')'];        
            
        otherwise
            coordstruct.g = ['@(x) (1-' vec2str(coordstruct.logtransf) ').*' z ... 
                ' + ' vec2str(coordstruct.logtransf) '.*' zlog ];
            coordstruct.ginv = ['@(y) (1-' vec2str(coordstruct.logtransf) ') .* (' vec2str(coordstruct.gamma) ' .* y + ' vec2str(coordstruct.mu) ') + ' ...
                vec2str(coordstruct.logtransf) ' .* exp(' vec2str(coordstruct.gamma) ' .* y + ' vec2str(coordstruct.mu) ')'];        
    end
    
    % Convert starting values to transformed coordinates
    coordstruct.g = str2func(coordstruct.g);
    coordstruct.ginv = str2func(coordstruct.ginv);
    
    if ischar(coordstruct.g); g = str2func(coordstruct.g); else g = coordstruct.g; end
    if ischar(coordstruct.ginv); ginv = str2func(coordstruct.ginv); else ginv = coordstruct.ginv; end
    
    % Check that the transform works correctly in the range
    t(1) = all(abs(ginv(g(lb)) - lb) < NumEps);
    t(2) = all(abs(ginv(g(ub)) - ub) < NumEps);
    t(3) = all(abs(ginv(g(plb)) - plb) < NumEps);
    t(4) = all(abs(ginv(g(pub)) - pub) < NumEps);    
    assert(all(t), 'Cannot invert the transform to obtain the identity at the provided boundaries.');
    
    coordstruct.lb = g(lb);
    coordstruct.ub = g(ub);
    coordstruct.plb = g(plb);
    coordstruct.pub = g(pub);
    
    varargout{1} = coordstruct;

    
end

%--------------------------------------------------------------------------
function s = vec2str(v)
% Convert numerical vector to string

MaxPrecision = 17;  % Maximum precision for a double
if size(v,1) > 1; transp = ''''; else transp = []; end
s = '[';
for i = 1:length(v)-1; s = [s, num2str(v(i),MaxPrecision), ',']; end
s = [s, num2str(v(end),MaxPrecision), ']' transp];




