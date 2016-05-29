function s = mat2struct(m, f, l)
%MAT2STRUCT Convert the contents of a numerical vector into a struct.
%   S = MAT2STRUCT(M, F, L) convert a numerical vector M into a struct
%   with field names taken from the cell array F. The vector L contains the
%   length of each field. Usually M, F, L are produced by STRUCT2MAT.
%
%	See also STRUCT2MAT

% For some reason the struct constructor fails if S is empty
s.e_m_p_t_Y_ = 0;
p = 1;

for i = 1:length(f)
    if p <= length(m)
        s = setfield(s, f{i}, m(p:min(p+l(i)-1, length(m))));
        p = p + l(i);
    else
        s = setfield(s, f{i}, []);
    end        
end

s = rmfield(s, 'e_m_p_t_Y_');

end