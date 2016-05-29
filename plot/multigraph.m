function hg = multigraph(g, borders, extborders)
%MULTIGRAPH Prepare a boxed, annotated graph.
%


if ~exist('borders', 'var'); borders = []; end
if ~exist('extborders', 'var'); extborders = []; end

% Size constants
if isempty(borders); borders = [0.075, 0.1]; end

if isempty(extborders); extborders = borders; end

xborder = borders(1);
yborder = borders(2);
xextborder = extborders(1);
yextborder = extborders(2);

fontsize = 20;
ticksize = 14;

if ~isnumeric(g) || any(g(:) < 0); error('The graph matrix G must contain integers from 1 to NG.'); end

% Number of subgraphs
ng = max(g(:));
grect = zeros(ng, 4);

% xborder = xborder*size(g, 2);
% yborder = yborder*size(g, 1);

% Check that all subgraphs exist and that they are rectangles
for i = 1:ng
    if all(g ~= i); error(['Subgraph #' num2str(i) ' not in the graph matrix G.']); end; 
    % Find the upper left corner
    for ii = 1:size(g, 1)
        for jj = 1:size(g, 2)
            if g(ii, jj) == i && grect(i, 1) == 0; grect(i, 1:2) = [ii jj]; end
            if g(ii, jj) == i; grect(i, 3:4) = [ii jj]; end
        end
    end    
    if any((grect(i, 3:4) - grect(i, 1:2)) < 0); error('The subgraphs in graph matrix G have inconsistent shapes.'); end 
end

for i = 1:ng
   hg(i) = drawSubgraph(i);
end


% if ~exist('xstr', 'var'); xstr = [] ; end
% if ~exist('ystr', 'var'); ystr = [] ; end
% if ~exist('titlestr', 'var'); titlestr = [] ; end
% if ~exist('fontsize', 'var'); fontsize = [] ; end
% if isempty(fontsize); fontsize = 16; end

% hold on;
% box on;

% if ~isempty(xstr); xlabel(xstr, 'Interpreter', 'latex', 'FontSize', fontsize); end
% if ~isempty(ystr); ylabel(ystr, 'Interpreter', 'latex', 'FontSize', fontsize); end
% if ~isempty(titlestr); title(titlestr, 'Interpreter', 'latex', 'FontSize', fontsize); end

    function h = drawSubgraph(nsub)
        r = grect(nsub, :);
        xslots = r(4) - r(2) + 1;
        yslots = r(3) - r(1) + 1;
        
        nx = size(g,2);
        ny = size(g,1);
        
        xpanelwidth = (1 - 2*xextborder - (nx-1)*xborder)/nx;
        ypanelwidth = (1 - 2*yextborder - (ny-1)*yborder)/ny;     
        
        % xpanelwidth = size(g, 2) + (size(g,2)-1)*xborder + 2*xextborder;
        % ypanelwidth = size(g, 1) + (size(g,1)-1)*yborder + 2*yextborder;
                               
        % xsize = (xslots + (xslots-1)*xborder)/xpanelwidth;
        % ysize = (yslots + (yslots-1)*yborder)/ypanelwidth;
        xpos = (r(2) - 1)*(xpanelwidth + xborder) + xextborder;
        ypos = (ny - r(3))*(ypanelwidth + yborder) + yextborder;
        
        xsize = xslots*xpanelwidth + (xslots-1)*xborder;
        ysize = yslots*ypanelwidth + (yslots-1)*yborder;
                
        h = axes('position', [xpos, ypos, xsize, ysize]);
              
        % xlabel('', 'Interpreter', 'latex', 'FontSize', fontsize);
        % ylabel('', 'Interpreter', 'latex', 'FontSize', fontsize);        
        xlabel('', 'FontSize', fontsize);
        ylabel('', 'FontSize', fontsize);        
        set(h, 'FontSize', ticksize);
    end

end