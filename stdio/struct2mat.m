function [m, f, l] = struct2mat(s)
%STRUCT2MAT Convert the contents of a struct into a single vector.
%   M = STRUCT2MAT(S) convert a struct S into a single numerical row vector.
%   Only numerical vector fields from the struct S are considered; 
%   multidimensional vectors (matrices) or other fields are ignored. The 
%   vectors in M keep the same order as in the original struct S.
%
%  [M, F, L] = STRUCTMAT(S) as above, but also returns a cell array of
%  strings S for the fields of the struct and a vector of lengths L,
%  specifying for each field the vector length in the original struct S.
%
%	See also MAT2STRUCT

m = []; f = []; l = [];

% Get the original field names
fn = fieldnames(s);

for i = 1:length(fn)
   c = s.(fn{i});
   if isnumeric(c) && isvector(c)
       if size(c, 1) > 1; c = c'; end
       m = [m c];
       l = [l length(c)];
       f{length(f)+1} = fn{i};       
   end
end

end