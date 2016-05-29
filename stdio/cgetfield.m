%CGETFIELD Conditionally get structure field contents.
%   F = CGETFIELD(S,'field') returns the contents of the specified
%   field.  This is equivalent to the syntax F = S.field.
%   S must be a 1-by-1 structure. If the specified field does not exists,
%   returns the empty array.
%  
%   F = CGETFIELD(S,'field',DEFAULT) returns the contents of the 
%   specified field, or DEFAULT if the field is absent or empty.
%
function f = cgetfield(S,field,default)

% Empty default
if nargin < 3; default = []; end

if isfield(S,field) && ~isempty(S.(field))
    f = S.(field);
else
    f = default;
end
