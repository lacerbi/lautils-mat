% FMINMULTI_TEST

camel6 = @(x) (4-2.1*x(1)^2+(x(1)^4)/3) * x(1)^2 + x(1)*x(2) + (-4+4*x(2)^2) * x(2)^2;
localopt.Algorithm = 'sqp';
options.Display = 'iter';
options.Method = {'feval','bps','patternsearch'};
options.OptOptions = {[],localopt,[]};
options.LoadFile = 'fminmulti_test.mat';
options.SaveFile = 'fminmulti_test.mat';
options.SaveTime = 5;
[x,fval] = fminmulti(camel6,[-3,-2],[3,2],[100,10,1],options)
