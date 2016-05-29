function tf = anystrcmpi(s1,s2)
%ANYSTRCMPI Compare sets of strings ignoring case.
%
%   TF = ANYSTRCMPI(S1,S2) returns logical 1 (true) if any of the strings 
%   in cell array S1 equals any of the strings in cell array S2, and 
%   returns logical 0 (false) otherwise. ANYSTRCMPI works also with simple
%   string inputs for S1 and S2.

if iscell(s1) && iscell(s2)
    temp = cellfun(@(x) strcmpi(s1,x), s2, 'UniformOutput', 0);
    tf = any(cellfun(@any, temp));
else
    tf = any(strcmpi(s1,s2));
end