function bayes1plot(x1,x2,f,n)

% if nargin < 4; lossfun = 'map'; end
if nargin < 4; n = 1000; end

y = linspace(x1, x2, n); 
xmean = zeros(1, n);
xmap = zeros(1, n);

for i = 1:length(y)
    xmean(i) = bayes1(x1, x2, @(z) f(z, y(i)), 'mean', n);
    xmap(i) = bayes1(x1, x2, @(z) f(z, y(i)), 'map', n);
end

hold on;
plot(y, xmean, 'r', 'LineWidth', 1);
plot(y, xmap, 'b', 'LineWidth', 1);
hold off;

fontsize = 20;
xlabel('Cue position $y$', 'Interpreter', 'LaTeX', 'FontSize', fontsize);
ylabel('Inferred position $\hat{x}$', 'Interpreter', 'LaTeX', 'FontSize', fontsize);

end