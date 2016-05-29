function scatterbars(x, y, xerr, yerr, s, c, m, ax)
%SCATTERBARS Scatter/bubble plot with error bars.
%   SCATTER(X,Y,XERR,YERR,S,C,T,AX) displays colored T at the locations 
%   specified by the vectors X and Y (which must be the same size), with
%   error bars of length XERR and YERR. XERR and YERR can be empty, 
%   scalars, or the same size as X and Y. Empty or zero-valued error bars
%   are not drawn. If they are scalar, the same error bar is applied
%   to all points.
%
%   S determines the area of each marker (in points^2) and the size of the
%   error bars endpoints. S can be a vector the same length a X and Y or a 
%   scalar. If S is a scalar, MATLAB draws all the markers the same size. 
%   If S is empty, the default size is used.
%   
%   C determines the colors of the markers. When C is a vector the
%   same length as X and Y, the values in C are linearly mapped
%   to the colors in the current colormap. When C is a 
%   length(X)-by-3 matrix, it directly specifies the colors of the  
%   markers as RGB values. C can also be a color string. See ColorSpec.
%
%   AX specifies the axes (if left empty, they are automatically set).  
%
%   SCATTERBARS(X,Y,XERR,YERR) draws the markers with specified error bars 
%   in default size and color. 
%
%   See also SCATTER.

% Error bars
if ~exist('xerr', 'var'); xerr = [] ; end
if isempty(xerr); xerr = zeros(1, length(x)); end
if isscalar(xerr); xerr = xerr*ones(1, length(x)); end

if ~exist('yerr', 'var'); yerr = [] ; end
if isempty(yerr); yerr = zeros(1, length(y)); end
if isscalar(yerr); yerr = yerr*ones(1, length(y)); end

% Marker size
if ~exist('s', 'var'); s = [] ; end
if isempty(s); s = 4; end

% Color
if ~exist('c', 'var'); c = [] ; end
if isempty(c); c = [0 0 1]; end

% Graph step size
if ~exist('step', 'var'); step = [] ; end
if isempty(step); step = 1; end

% Marker type
if ~exist('m', 'var'); m = [] ; end
if isempty(m); m = 'o'; end

% Axes size
if ~exist('ax', 'var'); ax = [] ; end

% Size consistency check
if length(xerr) ~= length(x) || length(yerr) ~= length(y)
   error('The length of the data vector does not match the length of the error vector.'); 
end

% Axes boundaries
if isempty(ax)
    xmin = min(x - xerr); xmax = max(max(x + xerr));
    ymin = min(y - yerr); ymax = max(max(y + yerr));
% xmin = floor(xmin/step)*step; xmax = ceil(xmax/step)*step;
% ymin = floor(ymin/step)*step; ymax = ceil(ymax/step)*step;
else
   xmin = ax(1); xmax = ax(2); ymin = ax(3); ymax = ax(4);    
end

%if xmax > xmin && ymax > ymin
    axis([xmin xmax ymin ymax]);
    sizes = [(xmax-xmin) (ymax-ymin)]/200;
%else
%    sizes = [1 1];
%end

% Plot paradigm points for legend
%for j = 1:nConds
%       pss = PSSdata(1, j, 1);
%       xpos = 1;
%       % Plot the central pss point
%       plot(1, x(j), 'o', 'LineWidth', 1, 'Color', colors(j, :), 'MarkerFaceColor', colors(j, :), 'MarkerSize', 8);
%end

for i = 1:length(x)    
    % Plot the x error line
    if xerr(i) > 0
        width = sqrt(s)*sizes(2);
        plot((x(i) + xerr(i))*[1 1], [y(i) - width,  y(i) + width], '-', 'LineWidth', 1, 'Color', c);
        plot((x(i) - xerr(i))*[1 1], [y(i) - width,  y(i) + width], '-', 'LineWidth', 1, 'Color', c);
        plot(x(i) + xerr(i)*[-1, 1], [y(i) y(i)], '-', 'LineWidth', 1, 'Color', c);        
    end
    
    % Plot the y error line
    if yerr(i) > 0
        width = sqrt(s)*sizes(1);
        plot([x(i) - width,  x(i) + width],(y(i) + yerr(i))*[1 1], '-', 'LineWidth', 1, 'Color', c);
        plot([x(i) - width,  x(i) + width],(y(i) - yerr(i))*[1 1], '-', 'LineWidth', 1, 'Color', c);
        plot([x(i) x(i)], y(i) + yerr(i)*[-1, 1], '-', 'LineWidth', 1, 'Color', c);
    end

    % Plot the central point
    plot(x(i), y(i), m, 'LineWidth', 1, 'Color', 'k', 'MarkerFaceColor', c, 'MarkerSize', s);
end

% Draw legend
%leglabel = [];
%for j = 1:nConds
%    leglabel{j} = [num2str(conds(j)) ' ms lag'];
%end

%leghandle = legend(leglabel);
%set(leghandle, 'Interpreter', 'latex');
%set(leghandle, 'FontSize', fontsize);

end