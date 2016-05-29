function hout=boxplot2(x,varargin)
%BOXPLOT Display boxplots of a data sample.
%   BOXPLOT(X) produces a box and whisker plot with one box for each column
%   of X.  The boxes have lines at the lower quartile, median, and upper
%   quartile values.  The whiskers are lines extending from each end of the
%   boxes to show the extent of the rest of the data.  Outliers are data
%   with values beyond the ends of the whiskers.
%
%   BOXPLOT(X,G) produces a box and whisker plot for the vector X grouped
%   by G.  G is a grouping variable defined as a categorical variable,
%   vector, string matrix, or cell array of strings.  G can also be a cell
%   array of several grouping variables (such as {G1 G2 G3}) to group the
%   values in X by each unique combination of grouping variable values.
%
%   BOXPLOT(...,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%      'notch'       'on' to include notches (default is 'off').
%      'symbol'      Symbol and color to use for all outliers (default is 'r+').
%      'orientation' Box orientation, 'vertical' (default) or 'horizontal'.
%      'whisker'     Maximum whisker length (default 1.5).
%      'labels'      Character array or cell array of strings containing
%                    labels for each column of X, or each group in G.
%      'colors'      A string or a three-column matrix of box colors.  Each
%                    box (outline, median line, and whiskers) is drawn in the
%                    corresponding color.  Default is to draw all boxes with
%                    blue outline, red median, and black whiskers.  Colors are
%                    recycled if necessary.
%      'widths'      A numeric vector or scalar of box widths.  Default is
%                    0.5, or slightly smaller for fewer than three boxes.
%                    Widths are recycled if necessary.  
%      'grouporder'  When G is given, a character array or cell array of
%                    group names, in the order in which the groups in G
%                    should be plotted.  Default is ordering the groups 
%                    the same as given in G.  Ignored when G is not given.
%      'positions'   A numeric vector of box positions (default is 1:n).
%      'plotstyle'   'compact' to render boxes more concisely, which is 
%                    useful when you have many groups.  Boxes are single thick
%                    lines, medians are dots, and outliers default to o
%                    markers.  The outliers are jittered (randomly displaced) 
%                    along the width of the box to minimize exact overlaps.  
%                    The notch and width arguments are ignored. The default 
%                    plotstyle is 'traditional'.
%
%   In a notched box plot the notches represent a robust estimate of the
%   uncertainty about the medians for box-to-box comparison.  Boxes whose
%   notches do not overlap indicate that the medians of the two groups
%   differ at the 5% significance level.  Note that the comparison may be 
%   misleading when the notches extend to the end of the box, as they may 
%   have been truncated to the end of the box; this may occur with small 
%   sample sizes.
%
%   Whiskers extend from the box out to the most extreme data value within 
%   WHIS*IQR, where WHIS is the value of the 'whisker' parameter and IQR 
%   is the interquartile range of the sample.
%
%   BOXPLOT(AX,...) plots into the axes with handle AX.
%
%   H = BOXPLOT(...) returns the handle H to the lines in the box plot.
%   H has one column per box, consisting of the handles for the various
%   parts of the box.  For the 'traditional' plotstyle, the 7 rows of 
%   handles are for the upper whisker, lower whisker, upper adjacent value, 
%   lower adjacent value, box, median, and outliers.  For the 'compact'
%   plotstyle, the 5 rows of handles are for the whisker, box, median
%   outer (ring), median inner (dot), and outliers.
%
%
%   Example:  Box plot of car gas mileage grouped by country
%      load carsmall
%      boxplot(MPG, Origin)
%      boxplot(MPG, Origin, 'sym','r*', 'colors',hsv(7))
%      boxplot(MPG, Origin, 'grouporder', ...
%                   {'France' 'Germany' 'Italy' 'Japan' 'Sweden' 'USA'})
%
%   Example: Plot by median gas mileage
%      [sortedMPG,sortedOrder] = sort(grpstats(MPG,Origin,@median));
%      pos(sortedOrder) = 1:6;
%      boxplot(MPG, Origin, 'position', pos)
%
%   See also ANOVA1, KRUSKALWALLIS, MULTCOMPARE.

%   Older syntax still supported:
%       BOXPLOT(X,NOTCH,SYM,VERT,WHIS)

%   References
%   [1] McGill, R., Tukey, J.W., and Larsen, W.A. (1978) "Variations of
%       Boxplots", The American Statistician, 32:12-16.
%   [2] Velleman, P.F. and Hoaglin, D.C. (1981), Applications, Basics, and
%       Computing of Exploratory Data Analysis, Duxbury Press.
%   [3] Nelson, L.S., (1989) "Evaluating Overlapping Confidence
%       Intervals", Journal of Quality Technology, 21:140-141.

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 2.15.4.24.2.1 $  $Date: 2008/01/30 20:55:29 $

whissw = 0; % 0=don't plot whisker inside the box.

if isscalar(x) && ishandle(x)
    ax = x;
    x = varargin{1};
    varargin(1) = [];
else
    ax = [];
end
error(nargchk(1+length(ax),Inf,nargin));

if isvector(x)
   % Might have one box, or might have a grouping variable. n will be properly
   % set later for the latter.
   x = x(:);
   n = 1; %
else
   % Have a data matrix, use as many boxes as columns.
   n = size(x,2);
end

% Detect if there is a grouping variable by looking at the second input
nargs = nargin - length(ax);
if nargs < 2
   g = [];
else
   g = varargin{1};
   if isempty(g) || isequal(g,1) || isequal(g,0) || (ischar(g) && size(g,1)==1)
      % It's a NOTCH value or a parameter name
      g = [];
   else
      % It's a grouping variable
      if ~isvector(x)
         error('stats:boxplot:VectorRequired',...
               'X must be a vector when there is a grouping variable.');
      end
      varargin(1) = [];
      nargs = nargs - 1;
   end
end

% Set defaults
notch  = 0;
sym    = '';
vert   = 1;
whis   = 1.5;
labels = {};
colors = []; % default is blue box, red median, black whiskers
posns  = []; % default is 1:n
widths = []; % default is 0.5, smaller for three or fewer boxes
grporder = []; % default is 1:n
plotstyle = 0;

% Determine if we have parameter names or the old syntax
if nargs > 1
   if ischar(varargin{1})
      okargs =   {'notch' 'symbol' 'orientation' 'whisker' 'labels' 'colors' 'positions' 'widths' 'grouporder' 'plotstyle'};
      defaults = { notch   sym      vert          whis      labels   colors   posns       widths   grporder     plotstyle};
      [eid,emsg,notch,sym,vert,whis,labels,colors,posns,widths,grporder,plotstyle] = ...
                                     statgetargs(okargs,defaults,varargin{:});
      if ~isempty(eid)
         error(sprintf('stats:boxplot:%s',eid),emsg);
      end
   else
      if (nargs>=2) && ~isempty(varargin{1}), notch = varargin{1}; end
      if (nargs>=3) && ~isempty(varargin{2}), sym   = varargin{2}; end
      if (nargs>=4) && ~isempty(varargin{3}), vert  = varargin{3}; end
      if (nargs>=5) && ~isempty(varargin{4}), whis  = varargin{4}; end
   end
end

% Convert wordy inputs to internal codes

if isequal(notch,'on')
   notch = 1;
elseif isempty(notch) || isequal(notch,'off')
   notch = 0;
elseif ~isscalar(notch) || ~ismember(notch,0:1)
   error('stats:boxplot:InvalidNotch','Invalid value for ''notch'' parameter');
end

if isempty(vert)
   vert = 1;
elseif ischar(vert)
   vert = strmatch(vert,{'horizontal' 'vertical'}) - 1;
end
if isempty(vert) || ~isscalar(vert) || ~ismember(vert,0:1)   
   error('stats:boxplot:InvalidOrientation',...
         'Invalid value for ''orientation'' parameter');
end

if ~isscalar(whis) || ~isnumeric(whis)
   error('stats:boxplot:BadWhisker',...
         'The ''whisker'' parameter value must be a numeric scalar.');
end

switch plotstyle
    case {0,'traditional'}
        plotstyle=0;
    case {1,'compact'}
        plotstyle=1;
    otherwise
        error('stats:boxplot:BadPlotStyle',...
            'Invalid value for ''plotstyle'' parameter.');
end

% Deal with grouping variable before processing more inputs
if ~isempty(g)
   if vert, sep = '\n'; else sep = ','; end
   [g,glabel,gname,multiline] = mgrp2idx(g,size(x,1),sep);
   n = size(gname,1);
   if numel(g) ~= numel(x)
      error('stats:boxplot:InputSizeMismatch',...
            'X and G must have the same length.');
   end
else
    multiline = false;
end

% Reorder the groups if necessary
if ~isempty(g) && ~isempty(grporder)
   if iscellstr(grporder) || ischar(grporder)
      % If we have a grouping vector, grporder may be a list of group names.
      if ischar(grporder), grporder = cellstr(grporder); end
      [dum,grporder] = ismember(glabel,grporder(:));
      % Must be a permutation of the group names
      if ~isequal(sort(grporder),(1:n)')
         error('stats:boxplot:BadOrder', ...
               'The ''grouporder'' parameter value must contain all the unique group names in G.');
      end
   else
      error('stats:boxplot:BadOrder', ...
            'The ''grouporder'' parameter value must be a character array or a cell array of strings.');
   end
   g = grporder(g);
   glabel(grporder) = glabel;
   gname(grporder,:) = gname;
end

% Process the rest of the inputs

if isempty(labels)
   if ~isempty(g)
      labels = glabel;
   end
else
   if ~(iscellstr(labels) && numel(labels)==n) && ...
      ~(ischar(labels) && size(labels,1)==n)
      % Must have one label for each box
      error('stats:boxplot:BadLabels','Incorrect number of box labels.');
   end
   if ischar(labels), labels = cellstr(labels); end
   multiline = false;
end
dfltLabs = (isempty(labels) && isempty(g)); % box labels are just column numbers


if isempty(widths)
   widths = repmat(min(0.15*n,0.5),n,1);
elseif ~isvector(widths) || ~isnumeric(widths) || any(widths<=0)
   error('stats:boxplot:BadWidths', ...
         'The ''widths'' parameter value must be a numeric vector of positive values.');
elseif length(widths) < n
   % Recycle the widths if necessary.
   widths = repmat(widths(:),ceil(n/length(widths)),1);
end
%if too many widths supplied, or recycling overshoots, truncate to correct length
if length(widths)>n
    widths = widths(1:n);
end

if isempty(colors)
   % Empty colors tells boxutil to use defaults.
   whiscolor= 'k';
   boxcolor='b';
   medcolor = 'r';
elseif ischar(colors) && isvector(colors)
    colors = colors(:); % color spec string, make it a column
    if size(colors,1) > 1 && size(colors,1) < n
        % Recycle the colors if necessary.
        colors = repmat(colors,ceil(n/size(colors,1)),1);
    end
    %if too many colors supplied, or if recycling overshoots, truncate to
    %correct length
    if size(colors,1)>n
        colors=colors(1:n);
    end
    whiscolor = colors;
    boxcolor = colors;
    medcolor = colors;
elseif isnumeric(colors) && (ndims(colors)==2) && (size(colors,2)==3)
   % RGB matrix, that's ok
    if size(colors,1) > 1 && size(colors,1) < n
       % Recycle the colors if necessary.
       colors = repmat(colors,ceil(n/size(colors,1)),1);
    end
    %if too many colors supplied, or if recycling overshoots, truncate to
    %correct length
    if size(colors,1)>n
        colors=colors(1:n,:);
    end
    whiscolor = colors;
    boxcolor = colors;
    medcolor = colors;
else
   error('stats:boxplot:BadColors',...
         'The ''colors'' parameter value must be a string or a three-column numeric matrix.');
end

[junklinestyle,markercolor1,markertype1,msg]=colstyle(sym);
if ~isempty(msg)
    error('stats:boxplot:BadSymbol',msg.message);
end

markercolor = markercolor1;
markertype = markertype1;
if strcmp(markertype,'')
    switch plotstyle
        case 0,  
            %If neither color nor symbol specified, default to +.
            %If color but not symbol specified in sym
            %then choose markertype = none (this is deprecated).

            if strcmp(markercolor1,'')
                markertype = '+';
            else
                markertype = 'n';
            end
        case 1,  markertype='o';
        otherwise, error('stats:boxplot:BadPlotStyleInternal',...
                'Invalid value used internally for plotstyle.');
    end
end

if strcmp(markercolor,'')
    switch plotstyle
        case 0,  
            %if neither color nor symbol specified, default to red
            %if symbol but not color specified, choose boxcolor (deprecated)
            if strcmp(markertype1,'')
                markercolor = 'r';
            else
                markercolor = boxcolor;
            end
        case 1,  markercolor=boxcolor;
        otherwise, error('stats:boxplot:BadPlotStyleInternal',...
            'Invalid value used internally for plotstyle.');
    end
end

if isempty(posns)
   posns = 1:n;
elseif ~isvector(posns) || ~isnumeric(posns)
   error('stats:boxplot:BadPositions', ...
         'The ''positions'' parameter value must be a numeric vector.');
elseif length(posns) ~= n
   % Must have one position for each box
   error('stats:boxplot:BadPositions', ...
         'The ''positions'' parameter value must have one element for each box.');
else
   [dum,ord] = sort(posns);
   if isempty(labels) % labels never empty when grouping vector supplied
       % If we have matrix data with no labels, and the positions are not 1:n,
       % we need to force the default column number tick labels 1:n, but in
       % the right positions.
       labels = cellstr(num2str(ord(:)));
   else
       % Permute the labels to match the plot position order.
       labels = labels(ord);
   end
end

%
% Done processing inputs
%

% Put at least the widest box or half narrowest spacing at each margin
if n > 1
    wmax = max(max(widths), 0.5*min(diff(posns)));
else
    wmax = 0.5;
end
xlims = [min(posns)-wmax, max(posns)+wmax];

ymin = nanmin(x(:));
ymax = nanmax(x(:));
if ymax > ymin
   dy = (ymax-ymin)/20;
else
   dy = 0.5;  % no data range, just use a y axis range of 1
end
ylims = [(ymin-dy) (ymax+dy)];



% Scale axis for vertical or horizontal boxes.
ax = newplot(ax);
set(ancestor(ax,'figure'),'CurrentAxes',ax); 
oldstate = get(ax,'NextPlot');
set(ax,'NextPlot','add','Box','on');

if isempty(xlims)
    xlims = [0 1];
end
if isempty(ylims)
    ylims = [0 1];
end
if vert
    axis(ax,[xlims ylims]);
    drawnow;
    nop = getappdata(ax,'NormalizedOuterPosition');
    if isempty(nop)
        % If a boxplot formerly occupied these axes, its outerposition may
        % already have been shrunk to accomodate multi-line group labels.
        % Store a new position if that isn't the case.
        setappdata(ax,'NormalizedOuterPosition',get(ax,'OuterPosition'));
    end
    set(ax,'XTick',sort(posns),'Units','normalized');
    ylabel(ax,'Values');
    if dfltLabs, xlabel(ax, 'Column Number'); end
else
    axis(ax,[ylims xlims]);
    set(ax,'YTick',sort(posns));
    xlabel(ax,'Values');
    if dfltLabs, ylabel(ax,'Column Number'); end
end
if nargout>0
   hout = [];
end

xvisible = NaN(size(x));
notnans = ~isnan(x);
whislo = NaN(n,1);
whishi = NaN(n,1);
boxlo=NaN(n,1);
boxhi=NaN(n,1);
notchlo=NaN(n,1);
notchhi=NaN(n,1);
med=NaN(n,1);
outliervals=cell(1,n);
for i= 1:n
   if ~isempty(g)
      thisgrp = find((g==i) & notnans);
   else
      thisgrp = find(notnans(:,i)) + (i-1)*size(x,1);
   end
   [outliers,whislo(i),whishi(i),outliervals{i},...
       boxlo(i),boxhi(i),notchlo(i),med(i),notchhi(i)] = ...
       boxutil(x(thisgrp),notch,whis,whissw);
   outliers = thisgrp(outliers);
   xvisible(outliers) = x(outliers);
end

if notch
    notchdepth=.5;
else
    notchdepth=0;
end

switch plotstyle
    case 0
    hh=boxrenderer(ax,posns,widths,vert, ...
        'lineAlongResponse',{'locationstart',boxhi,'locationend',whishi,...
            'linestyle','--','linewidth',.5,'linecolor',whiscolor,...
            'tag','Upper Whisker'},...
        'lineAlongResponse',{'locationstart',whislo,'locationend',boxlo,...
            'linestyle','--','linewidth',.5,'linecolor',whiscolor,...
            'tag','Lower Whisker'},...
        'lineAlongFactor',{'location',whishi,'linelength',.5,...
            'linestyle','-','linewidth',.5,'linecolor',whiscolor,...
            'tag','Upper Adjacent Value'},...
        'lineAlongFactor',{'location',whislo,'linelength',.5,...
            'linestyle','-','linewidth',.5,'linecolor',whiscolor,...
            'tag','Lower Adjacent Value'},...
        'lineBoxNotched',{'locationstart',boxlo,'locationend',boxhi,...
            'notchstart',notchlo,'notchmiddle',med,'notchend',notchhi,...
            'notchdepth',notchdepth,'linestyle','-','linewidth',.5,'linecolor',boxcolor,...
            'tag','Box'},...
        'lineAlongFactor',{'location',med,'linelength',1-notchdepth,...
            'linestyle','-','linewidth',.5,'linecolor',medcolor,...
            'tag','Median'},...
        'marker',{'location',outliervals,'jitter',0,'markertype',markertype,...
            'markersize',6,'markercolor',markercolor,...
            'tag','Outliers'} );
    case 1
    hh=boxrenderer(ax,posns,widths,vert, ...
        'lineAlongResponse',{'locationstart',whislo,'locationend',whishi,...
            'linestyle','-','linewidth',.5,'linecolor',boxcolor,...
            'tag','Whisker'},...
        'lineAlongResponse',{'locationstart',boxlo,'locationend',boxhi,...
            'linestyle','-','linewidth',4,'linecolor',boxcolor,...
            'tag','Box'},...
        'marker',{'location',med,'jitter',0,'markertype','o',...
            'markersize',6,'markercolor',boxcolor,'markerfill','b',...
            'tag','MedianOuter'}, ...
        'marker',{'location',med,'jitter',0,'markertype','.',...
            'markersize',6,'markercolor','k',...
            'tag','MedianInner'}, ...
        'marker',{'location',outliervals,'jitter',.5,'markertype',markertype,...
            'markersize',4,'markercolor',markercolor,...
            'tag','Outliers'} ...
            );

end
      
if nargout>0
    hout = hh;
end

if ~isempty(labels)
   if multiline && vert
      % Turn off tick labels and axis label
      set(ax, 'XTickLabel','');
      setappdata(ax,'NLines',size(gname,2));
      xlabel(ax,'');
      ylim = get(ax, 'YLim');
      
      % Place multi-line text approximately where tick labels belong
      ypos = repmat(ylim(1),size(posns));
      text(posns,ypos,labels,'HorizontalAlignment','center', ...
                             'VerticalAlignment','top', 'UserData','xtick');

      % Resize function will position text more accurately
      f = ancestor(ax,'figure');
      set(f, 'ResizeFcn', @resizefcn, ...
               'Interruptible','off', 'PaperPositionMode','auto');
      resizefcn(f);
   elseif vert
      set(ax, 'XTickLabel',labels);
   else
      set(ax, 'YTickLabel',labels);
   end
end
set(ax,'NextPlot',oldstate);

% Store information for gname function
set(ax, 'UserData', {'boxplot' xvisible g vert});

end 

%=============================================================================

function [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = boxutil(x,notch,whis,whissw)
%BOXUTIL Produces a single box plot.

% define the median and the quantiles
pctiles = prctile(x,[25;50;75]);
q1 = pctiles(1,:);
med = pctiles(2,:);
q3 = pctiles(3,:);

% find the extreme values (to determine where whiskers appear)
vhi = q3+whis*(q3-q1);
% upadj = max(x(x<=vhi));
% if (isempty(upadj)), upadj = q3; end
upadj = prctile(x, 97.5);

vlo = q1-whis*(q3-q1);
% loadj = min(x(x>=vlo));
% if (isempty(loadj)), loadj = q1; end
loadj = prctile(x, 2.5);

outlier = x<loadj | x > upadj;
yy = x(outlier);

if whissw == 0
   upadj = max(upadj,q3);
   loadj = min(loadj,q1);
end

if notch
    n1 = med + 1.57*(q3-q1)/sqrt(length(x));
    n2 = med - 1.57*(q3-q1)/sqrt(length(x));
    %prevent notches from extending past edge of box
    if n1>q3, n1 = q3; end
    if n2<q1, n2 = q1; end
else
    n1=med;
    n2=med;
end

end
%=============================================================================

function resizefcn(f,dum) %#ok<INUSD>

% Adjust figure layout to make sure labels remain visible
h = findobj(f, 'UserData','xtick');
if (isempty(h))
   set(f, 'ResizeFcn', '');
   return;
end

% Loop over all axes
allax = findall(f,'Type','Axes');
for j=1:length(allax)
    ax = allax(j);
    nlines = getappdata(ax, 'NLines');
    
    if ~isempty(nlines)
        % Try to retain the original normalized outer position
        nop = getappdata(ax,'NormalizedOuterPosition');
        set(ax,'Units','normalized');
        if isempty(nop)
            nop = get(ax,'OuterPosition');
        end
        set(ax,'OuterPosition',nop);
        
        % Adjust position so the fake X tick labels have room to display
        temp = hgconvertunits(f,[0 0 1 1],'character','normalized',f);
        charheight = temp(4);
        li = get(ax,'LooseInset');
        tickheight = min((nlines+1)*charheight, nop(4)/2);
        p = [nop(1)+li(1)*nop(3),    nop(2)+tickheight, ...
             nop(3)*(1-li(1)-li(3)), nop(4)*(1-li(4))-tickheight];
        p = max(p,0.0001);
        set(ax, 'Position', p);
        
        % The following lines leave the axes in a state such that MATLAB
        % will try to preserve the outerposition rather than the ordinary
        % position
        op = get(ax,'OuterPosition');
        op([1 3]) = nop([1 3]);   % no need to change x positions
        if op(2)<0  % fix if current position yielded a negative value
            op(4) = min(1, op(4) + op(2));
            op(2) = 0;
        end
        set(ax, 'OuterPosition', op);
    end
end


% Position the labels at the proper place
xl = get(ax, 'XLabel');
set(xl, 'Units', 'data');
p = get(xl, 'Position');
ylim = get(ax, 'YLim');
p2 = (p(2)+ylim(1))/2;
for j=1:length(h)
   p = get(h(j), 'Position') ;
   p(2) = p2;
   set(h(j), 'Position', p);
end

end

function hgg=boxrenderer(ax, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin)

if isempty(ax)
    ax=gca;
end
%factorCenterPosition has one element per group
factorCenterPosition = factorCenterPosition(:);
boxDataWidth = boxDataWidth(:);

numgroups = size(factorCenterPosition,1);
if ~ismember(length(boxDataWidth),[1 numgroups])
	error ('stats:boxplot:badBoxDataWidthSize', ...
		'boxDataWidth must be either a scalar or a vector numgroups long');
end

if mod(length(varargin),2) == 1 || ~all(iscellstr(varargin(1:2:end))) ||  ~all(cellfun(@iscell,varargin(2:2:end)))
	error ('stats:boxplot:badSubFuncHandle', ...
		'Parameter names must be strings, and sub-argument lists must be cell arrays');
end

numsubfuncs = length(varargin)/2;

h=NaN(numgroups,numsubfuncs);

for subfuncindex = 1:numsubfuncs
% %%%%%%%%%%%%%%%%
% The next line lists the supported sub-functions
% %%%%%%%%%%%%%%%%%%%%
	if ~ismember(varargin{(subfuncindex-1)*2+1}, ...
		{'marker','lineAlongResponse','lineAlongFactor', ...
		'lineBox','lineBoxNotched'} )
% %%%%%%%%%%%%%%%%%%%%%%%
		error('stats:boxplot:badSubFuncHandle', 'Unrecognized subfunction handle');
	end
	funcname = ['draw_' varargin{(subfuncindex-1)*2+1}];
	h(:,subfuncindex) = feval(funcname,ax, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin{(subfuncindex-1)*2+2});
end

hgg = h';

end

%% draw_marker function
%used for outliers, and for medians drawn as one dot
%
%scalar args apply to all groups
%vector args have one item per group
%cell arrays of vectors have one item per each item in the group
function h=draw_marker(hgg, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin) %#ok<DEFNU>

%varargout = parseSubfuncArgs(argValues, argNames, defaults)
[     location,  jitter,  markertype,  markersize,  markercolor,  markerfill,  tag] = ...
    parseSubfuncArgs(varargin, ...
    {'location','jitter','markertype','markersize','markercolor','markerfill','tag'}, ...
     {[],        0,       'o',         6,          [0 0 1],      'n',         ''  }  );

h=NaN(numgroups,1);
if isempty(location)
    %kick out if no points to plot 
    return;
end

%complain about args that should not have a vector worth of values for each
%group
if iscell(markertype) || iscell(markersize) || iscell(markercolor) || iscell(markerfill)
    error('stats:boxplot:TooManyLineStylesPerGroup',...
        'You can have only one marker style per group');
end

%render as one line per group
for currentGroup = 1:numgroups
    %unpack data for this group, put into l_ variables
    %This unpacking should be put into its own subfunction or something
    %unpack bundle start
    [l_location,l_jitter,l_markertype,l_markersize,l_markerfill,l_boxDataWidth] = ...
        unpackGroup(currentGroup, ...
       location,  jitter,  markertype,  markersize,  markerfill,  boxDataWidth);
    numInGroup = length(l_location);
    if numInGroup == 0
        continue; %don't bother with the rest if this group is empty
    end

    l_markercolor = unpackGroupColor(currentGroup,markercolor);
    if strcmp(l_markertype,'n')
        l_markertype = 'none';
    end

    %compute factor axis position for each point
    %the random number will need to be replaced later with a hash
    %based on some data index, for the repeatability needed for
    %brushing
    if l_jitter==0
        jitter_offset=zeros(numInGroup,1);
    else
        jitter_offset = l_jitter .* (rand(numInGroup,1)-.5);
    end
    factor_value = factorCenterPosition(currentGroup) + ...
        l_boxDataWidth .* jitter_offset;

    if factorsOnXAxis==1
        x=factor_value;
        y=l_location;
    else
        x=l_location;
        y=factor_value;
    end
    switch(l_markerfill)
        case 'n', markerfacecolor = 'none';
        case 'b', markerfacecolor = 'auto';
        case 'f', markerfacecolor = l_markercolor;
        otherwise, 
            error('stats:boxplot:BadMarkerFill',...
                'Marker Fill must be n, b, or f, for none, background, or filled');
 
    end
    h(currentGroup) = line(x,y, ...
        'marker',l_markertype,...
        'markerSize',l_markersize,...
        'markerEdgeColor',l_markercolor,...
        'lineStyle','none',...
        'markerfacecolor',markerfacecolor,...
        'tag',tag,...
        'parent',hgg);

end

 
end

%% draw_lineAlongResponse
%used for whiskers and fixed pixel width boxes
% Width is in pixels
%

function h=draw_lineAlongResponse(hgg, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin) %#ok<DEFNU>

%varargout = parseSubfuncArgs(argValues, argNames, defaults)
[     locationStart,  locationEnd,  lineStyle,  lineWidth,  lineColor, tag] = ...
    parseSubfuncArgs(varargin, ...
    {'locationstart','locationend','linestyle','linewidth','linecolor','tag'}, ...
    {[],              [],          '-',        .5,         [0 0 1],    ''}  );

h=NaN(numgroups,1);
if isempty(locationStart) 
    %kick out if no lines to plot 
    return;
end

%complain about args that should not have a vector worth of values for each
%group
if iscell(lineStyle) || iscell(lineWidth) || iscell(lineColor)
    error('stats:boxplot:TooManyLineStylesPerGroup',...
        ['You can have only one line width and color per group, '...
        ' and one line style for all the groups']);
end

%render as one line per group
for currentGroup = 1:numgroups
    %unpack data for this group, put into l_ variables
    %This unpacking should be put into its own subfunction or something
    %unpack bundle start
    [l_locationStart,l_locationEnd,l_lineWidth,l_boxDataWidth] = ...
        unpackGroup(currentGroup, ...
       locationStart,  locationEnd,  lineWidth,  boxDataWidth); 
    l_boxDataWidth;  %#ok<VUNUS> This line selectively disables the lint 
    %warning from the function above - l_boxDataWidth may come in handy

    numInGroup = length(l_locationStart);
    if numInGroup == 0
        continue; %don't bother with the rest if this group is empty
    end

    l_lineColor = unpackGroupColor(currentGroup,lineColor);

    factor = nan(3*numInGroup-1,1);
    response = nan(3*numInGroup-1,1);
    factor(1:3:end)=factorCenterPosition(currentGroup);
    response(1:3:end) = l_locationStart;
    factor(2:3:end)=factorCenterPosition(currentGroup);
    response(2:3:end) = l_locationEnd;

    if factorsOnXAxis==1
        x=factor;
        y=response;
    else
        x=response;
        y=factor;
    end

    h(currentGroup) = line(x,y, ...
        'Marker','none',...
        'Color',l_lineColor,...
        'LineWidth',l_lineWidth,...
        'LineStyle',lineStyle,...
        'Tag',tag,...
        'Parent',hgg);

end

end


%% draw_lineAlongFactor
%used for medians drawn as lines, and for whisker ends
function h=draw_lineAlongFactor(hgg, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin) %#ok<DEFNU>

%varargout = parseSubfuncArgs(argValues, argNames, defaults)
[     location,  lineLength,  lineStyle,  lineWidth,  lineColor,  tag] = ...
    parseSubfuncArgs(varargin, ...
    {'location','linelength','linestyle','linewidth','linecolor','tag'}, ...
    { [],        1,          '-',        .5,         [0 0 1],    ''}  );
      
h=NaN(numgroups,1);

if isempty(location)
    return;
end

%complain about args that should not have a vector worth of values for each
%group
if iscell(lineStyle) || iscell(lineWidth) || iscell(lineColor)
    error('stats:boxplot:TooManyLineStylesPerGroup',...
        ['You can have only one line width and color per group, '...
        ' and one line style for all the groups']);
end


%render as one line per group
for currentGroup = 1:numgroups
    %unpack data for this group, put into l_ variables
    %This unpacking should be put into its own subfunction or something
    %unpack bundle start
    [l_location,l_lineLength,l_lineWidth,l_boxDataWidth] = ...
        unpackGroup(currentGroup, ...
       location,  lineLength,  lineWidth,  boxDataWidth);
    numInGroup = length(l_location);
    if numInGroup == 0
        continue; %don't bother with the rest if this group is empty
    end

    l_lineColor = unpackGroupColor(currentGroup,lineColor);

    factor = nan(3*numInGroup-1,1);
    response = nan(3*numInGroup-1,1);
    factor(1:3:end)=factorCenterPosition(currentGroup)-l_boxDataWidth.*l_lineLength./2;
    response(1:3:end) = l_location;
    factor(2:3:end)=factorCenterPosition(currentGroup)+l_boxDataWidth.*l_lineLength./2;
    response(2:3:end) = l_location;

    if factorsOnXAxis==1
        x=factor;
        y=response;
    else
        x=response;
        y=factor;
    end

    h(currentGroup) = line(x,y, ...
        'Marker','none',...
        'Color',l_lineColor,...
        'LineWidth',l_lineWidth,...
        'LineStyle',lineStyle,...
        'Tag',tag,...
        'Parent',hgg);

end


end



%% draw_lineBox
%used for data width filled boxes
function h=draw_lineBox(hgg, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin) %#ok<DEFNU>

%varargout = parseSubfuncArgs(argValues, argNames, defaults)
[     locationStart,  locationEnd,  lineStyle,  lineWidth,  lineColor,  tag] = ...
    parseSubfuncArgs(varargin, ...
    {'locationstart','locationend','linestyle','linewidth','linecolor','tag'}, ...
    {[],              [],          '-',         .5,         [0 0 1],   ''}  );
      

h=NaN(numgroups,1);

if isempty(locationStart)
    return;
end

%complain about args that should not have a vector worth of values for each
%group
if iscell(lineStyle) || iscell(lineWidth) || iscell(lineColor)
    error('stats:boxplot:TooManyLineStylesPerGroup',...
        ['You can have only one line width and color per group, '...
        ' and one line style for all the groups']);
end


%render as one line per group
for currentGroup = 1:numgroups
    %unpack data for this group, put into l_ variables
    %This unpacking should be put into its own subfunction or something
    %unpack bundle start
    [l_locationStart,l_locationEnd,l_lineWidth,l_boxDataWidth] = ...
        unpackGroup(currentGroup, ...
       locationStart,  locationEnd,  lineWidth,  boxDataWidth);
    numInGroup = length(l_locationStart);
    if numInGroup == 0
        continue; %don't bother with the rest if this group is empty
    end

    l_lineColor = unpackGroupColor(currentGroup,lineColor);

    factor = nan(6*numInGroup-1,1);
    response = nan(6*numInGroup-1,1);
    %lower left (assuming factorsOnXAxis == 0)
    factor(1:6:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(1:6:end) = l_locationStart;
    %upper left
    factor(2:6:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(2:6:end) = l_locationEnd;
    %upper right
    factor(3:6:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(3:6:end) = l_locationEnd;
    %lower right
    factor(4:6:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(4:6:end) = l_locationStart;
    %lower left, same as starting point
    factor(5:6:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(5:6:end) = l_locationStart;

    if factorsOnXAxis==1
        x=factor;
        y=response;
    else
        x=response;
        y=factor;
    end

    h(currentGroup) = line(x,y, ...
        'Marker','none',...
        'Color',l_lineColor,...
        'LineWidth',l_lineWidth,...
        'LineStyle',lineStyle,...
        'tag',tag,...
        'Parent',hgg);

end

end



%% draw_lineBoxNotched
%used for data width filled boxes
function h=draw_lineBoxNotched(hgg, numgroups, factorCenterPosition, boxDataWidth, factorsOnXAxis, varargin) %#ok<DEFNU>

%varargout = parseSubfuncArgs(argValues, argNames, defaults)
[     locationStart,  locationEnd,  notchStart,  notchMiddle,  notchEnd, ...
        notchDepth,  lineStyle,  lineWidth,  lineColor,  tag] = ...
    parseSubfuncArgs(varargin, ...
    {'locationstart','locationend','notchstart','notchmiddle','notchend',...
        'notchdepth','linestyle','linewidth','linecolor','tag'}, ...
    { [],             [],           [],          [],           [], ...
         .5,         '-',         .5,        [0 0 1],    ''}  );
      

h=NaN(numgroups,1);

if isempty(locationStart)
    return;
end

%complain about args that should not have a vector worth of values for each
%group
if iscell(lineStyle) || iscell(lineWidth) || iscell(lineColor)
    error('stats:boxplot:TooManyLineStylesPerGroup',...
        ['You can have only one line width and color per group, '...
        ' and one line style for all the groups']);
end

%render as one line per group
for currentGroup = 1:numgroups
    %unpack data for this group, put into l_ variables
    %This unpacking should be put into its own subfunction or something
    %unpack bundle start
    [l_locationStart,l_locationEnd,l_notchStart,l_notchMiddle,l_notchEnd,...
        l_notchDepth,l_lineWidth,l_boxDataWidth] = ...
        unpackGroup(currentGroup, ...
       locationStart,  locationEnd,  notchStart,  notchMiddle,  notchEnd, ...
          notchDepth,  lineWidth,  boxDataWidth);
    numInGroup = length(l_locationStart);
    if numInGroup == 0
        continue; %don't bother with the rest if this group is empty
    end

    l_lineColor = unpackGroupColor(currentGroup,lineColor);


    factor = nan(12*numInGroup-1,1);
    response = nan(12*numInGroup-1,1);
    %left middle of box, at the median (assuming factorsOnXAxis ==0)
    factor(1:12:end)=factorCenterPosition(currentGroup) ...
        -(1-l_notchDepth).*l_boxDataWidth./2;
    response(1:12:end) = l_notchMiddle;
    %left of box, at top of notch
    factor(2:12:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(2:12:end) = l_notchEnd;
    %left of box, at top of box
    factor(3:12:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(3:12:end) = l_locationEnd;
    %right of box, at top of box
    factor(4:12:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(4:12:end) = l_locationEnd;
    %right of box, at top of notch
    factor(5:12:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(5:12:end) = l_notchEnd;
    %right middle of box, at the median
    factor(6:12:end)=factorCenterPosition(currentGroup) ...
        +(1-l_notchDepth).*l_boxDataWidth./2;
    response(6:12:end) = l_notchMiddle;
    %right of box, at bottom of notch
    factor(7:12:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(7:12:end) = l_notchStart;
    %right of box, at bottom of box
    factor(8:12:end)=factorCenterPosition(currentGroup)+l_boxDataWidth./2;
    response(8:12:end) = l_locationStart;
    %left of box, at bottom of box
    factor(9:12:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(9:12:end) = l_locationStart;
    %left of box, at bottom of notch
    factor(10:12:end)=factorCenterPosition(currentGroup)-l_boxDataWidth./2;
    response(10:12:end) = l_notchStart;
    %left of box, at middle of notch
    factor(11:12:end)=factorCenterPosition(currentGroup)...
        -(1-l_notchDepth).*l_boxDataWidth./2;
    response(11:12:end) = l_notchMiddle;



    if factorsOnXAxis==1
        x=factor;
        y=response;
    else
        x=response;
        y=factor;
    end

    h(currentGroup) = line(x,y, ...
        'Marker','none',...
        'Color',l_lineColor,...
        'LineWidth',l_lineWidth,...
        'LineStyle',lineStyle,...
        'Tag',tag,...
        'Parent',hgg);

end

end



%%
%expect color to be a length 3 vector, a nx3 matrix (where n is the number of
%groups), or an n-long cell array of mx3 matrices (where m is the number of points 
%in the given group).
function groupColor = unpackGroupColor(groupindex, colorInputArg)
if iscell(colorInputArg)
    groupColor = colorInputArg{groupindex};
    if size(groupColor,1)~=3
        error('stats:boxplot:ColorMustBeLength3Vectors', ...
            'Color must be specified as RGB triples');
    end
elseif ischar(colorInputArg)
    if isscalar(colorInputArg)
        groupColor= colorInputArg;
    elseif isvector(colorInputArg)
        groupColor=colorInputArg(groupindex);
    else
        error('stats:boxplot:unexpectedDataType', ...
            ['color specified as character must be a single character '...
            'or a character array numgroups long']);
    end
else
    sz=size(colorInputArg);
    if ~ismember(3,sz)
        error('stats:boxplot:ColorMustBeLength3Vectors', ...
            'Color must be specified as RGB triple vectors');
    end
    if ismember(1,sz)
        groupColor = colorInputArg;
    else
        if sz(2)~=3
            error('stats:boxplot:MultipleColorsMatrixShape', ...
                'Multiple color must be specified as a 3 column matrix of RGB triples');
        end
        groupColor = colorInputArg(groupindex,:);
    end
end

end

%%
function varargout = unpackGroup (groupindex,varargin)
numvars = length(varargin);

%handle scalar, vector, or cell array of vectors
%verify that all items in this group are the same length, or scalar
if any(cellfun(@iscell,varargin))
    vectlens = ones(numvars+1,1);
    for i=1:numvars
        if iscell(varargin{i})
            vectlens(i) = length(varargin{i}{groupindex});
        end
    end
    %return all empties if one of the inputs is empty
    if min(vectlens)==0
        varargout{numvars}=[]; 
        return;
    end
    %error check
    if length(unique(vectlens))>2
        error('stats:boxplot:vectLenMismatchWithinGroup', ...
            ['when passing in multiple values per group, ', ...
            'there must be the same number of values in each group']);
    end
    %vectlen = max(vectlens);
end

%unpack data for each group, return as a vector for each group
for i=1:numvars
    currentvar = varargin{i};
    if iscell(currentvar)
        varargout{i}=currentvar{groupindex};
    elseif isscalar(currentvar)
        varargout{i}= currentvar;
    elseif isvector(currentvar)
        varargout{i}=currentvar(groupindex);
    else
        error('stats:boxplot:unexpectedDataType', ...
            ['Subfunction arg must be a scalar, a vector, ' , ...
            'or a cell array of vectors']);
    end
end


end


%%
function varargout = parseSubfuncArgs(argValues, argNames, defaults)
numargs = length(defaults);

if ischar(argValues{1}{1})
    varargout{numargs}=[];
	%parse as parameter-value pairs, permitting optional args
	[eid,emsg,varargout{:}] = ...
		statgetargs(argNames,defaults,argValues{1}{:}); 
	if ~isempty(eid) 
		error(sprintf('stats:boxplot:%s',eid),emsg); 
	end 
else
	%parse as fixed parameters, with no optional args and minimal error checking
    varargout = argValues{1};
end

% Do an error check on the args - vectors and cell inputs must be numgroups
% long; scalars are also ok.  
% A tweak is needed for color and linestyle
% length 0 is considered invalid, and will cause an error.
numels = ones(numargs+1,1);
for i=1:numargs
    numels(i)=length(varargout{i});
end
%tweak for certain color cases...
%avoid spurious warning for color specifications, if 1 or 2 colors are
%specificed
% 1 color is a 3x1 or 1x3 vector, 2 colors are a 2x3 matrix.
%these ought to be counted as 1 and 2 respectively, but are instead counted
%as 3 by the code above.
%Do this after the initial check to avoid string matching for all the args
maybeColorArgs = find(numels==3);
if ~isempty(maybeColorArgs)
    for i=maybeColorArgs'
        %c might be cap or lower case, so don't search on it
        %cells and char shortcut names are handled fine with the general
        %case
        if ~isempty(strfind(argNames{i},'olor')) && ~iscell(varargout{i}) && ~ischar(varargout{i})
            sz = size(varargout{i}); %already know one dim is size 3
            if sz(1)==1 || sz(2)==1
                numels(i)=1;
            end
            if sz(1)==2
                numels(i)=2;
            end
        end
    end
end
%tweak for linestyle string; lineStyle must be a scalar
LineStyleArg = strmatch('linestyle',argNames,'exact');
if ~isempty(LineStyleArg)
    lineStyleTemp = ismember(varargout{LineStyleArg},{'-','--',':','-.','none'});
    if length(lineStyleTemp)>1 || ~lineStyleTemp
        error('stats:boxplot:TooManyLineStylesPerAxis',...
            'You can have only one line style for all the groups');
    else
        numels(LineStyleArg)=1; %make the error check pass
    end
end
%tweak for tag string; tag must be a string
tagArg = strmatch('tag',argNames,'exact');
if ~isempty(tagArg)
    tagval = varargout{tagArg};
    if ischar(tagval)
        numels(tagArg)=1;
    else
        error('stats:boxplot:TagMustBeString',...
            'Tag must be set to a string, or '''' ');
    end
end
    
%do the check
if (length(unique(numels)))>2
    error ('stats:boxplot:unequalNumGroups','Vectors and cell arrays of vectors must be uniform length');
end


end

