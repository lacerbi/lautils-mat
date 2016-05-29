function options = setoptions(options,field,value,deflt)
%SETOPTIONS Set field of options structure
%
%   OPTIONS = SETOPTIONS(OPTIONS,FIELD,VALUE) set field FIELD of options
%   structure OPTIONS to VALUE, overwriting previous values.
%
%   OPTIONS = SETOPTIONS(OPTIONS,FIELD,VALUE,1) consider the setting as a
%   default setting -- only overwrites previous default settings, but not
%   previously user-set values.

if nargin < 4 || isempty(deflt); deflt = 0; end

% If DEFLT is 1, do not modify field if previously set by user
if ~deflt || ~isfield(options,'userset_') || ~any(strcmp(field, options.userset_))
    options.(field) = value;
    if ~deflt   % User-specified setting
        if ~isfield(options,'userset_'); options.userset_ = []; end
        if ~any(strcmp(field, options.userset_))
            options.userset_{end+1} = field;
        end
    end

    % Put userset_ as last field
    if isfield(options,'userset_')
        list = fieldnames(options);
        index = find(strcmp(list, 'userset_'),1);
        list(index) = [];
        list{end+1} = 'userset_';
        options = orderfields(options,list);
    end
end

