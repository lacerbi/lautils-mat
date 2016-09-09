function stdfig(fig)
%STDFIG Standardize figure

if nargin < 1 || isempty(fig); fig = gcf; end

screensize = get(groot,'ScreenSize');

c = fig.Children;
for i = 1:numel(c)
   if ~isa(c(i),'matlab.graphics.axis.Axes'); continue; end 
   set(c(i),'Box','off','TickDir','out');    
end

set(fig,'Color','w');
set(fig,'Position',screensize);
