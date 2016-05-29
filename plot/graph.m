function graph(xstr, ystr, titlestr, fontsize, latex)
%GRAPH Prepare a boxed, annotated graph.
%
% GRAPH(XSTR, YSTR, TITLESTR, FONTSIZE) creates a graph with x-axis caption 
% XSTR, y-axis caption YSTR and title TITLESTR, with font size FONTSIZE.
% If no string or an empty string is provided for a caption, no caption is written.
%
% GRAPH(XSTR, YSTR, TITLESTR) writes the captions in  standard font size.
% 

if ~exist('xstr', 'var'); xstr = [] ; end
if ~exist('ystr', 'var'); ystr = [] ; end
if ~exist('titlestr', 'var'); titlestr = [] ; end
if ~exist('fontsize', 'var'); fontsize = [] ; end
if isempty(fontsize); fontsize = 16; end
if ~exist('latex', 'var'); latex = [] ; end
if isempty(latex); latex = 1; end

hold on;
box on;

switch latex;     
    case 1; interp = 'latex'; 
    case 2; interp = 'tex'; 
    otherwise; interp = 'none';
end

if ~isempty(xstr); xlabel(xstr, 'Interpreter', interp, 'FontSize', fontsize); end
if ~isempty(ystr); ylabel(ystr, 'Interpreter', interp, 'FontSize', fontsize); end
if ~isempty(titlestr); title(titlestr, 'Interpreter', interp, 'FontSize', fontsize); end


end