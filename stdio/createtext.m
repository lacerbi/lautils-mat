% CREATETEXT Creates a text matrix, fill with space.
function txt = createtext(varargin)

nl = length(varargin);

maxl = 0;
for i = 1:nl
   if length(varargin{i}) > maxl; maxl = length(varargin{i}); end    
end

% Create a matrix full of 'spaces'
txt = ones(nl, maxl)*32;

for i = 1:nl
   txt(i, 1:length(varargin{i})) = varargin{i};
end

% Convert back to string
txt = char(txt);

end