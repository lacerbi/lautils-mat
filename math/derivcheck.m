function derivcheck(f,x)
%DERIVCHECK Check analytical vs numerical differentiation for a function

tic
%dy_num = fgrad(f,x,'five-points');
dy_num = gradest(f,x);
toc
tic
[y,dy_ana] = f(x);
toc

fprintf('Relative errors:\n');
(dy_num - dy_ana)./dy_num

fprintf('Absolute errors:\n');
dy_num - dy_ana

end