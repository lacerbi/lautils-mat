function equalticklength(h,ticklength)
%EQUALTICKLENGTH Make all tick lengths equal in all figure panels.


if nargin < 2 || isempty(ticklength); ticklength = 0.0035; end

ng = numel(h);

for g = 1:ng
    rect = get(h(g),'Position');
    height = rect(4);
    width = rect(3);
    bottom = rect(2);
    left = rect(1);
    
    %if square(g)
    %    axislen = min([height,width]);        
    %    ticklen = ticklength./axislen^2;        
    %    axislen = height;
    %else
    
    hfig = get(h(g),'Parent');
    rectfig = get(hfig, 'Position');
    width_px = rect(3) * rectfig(3);
    height_px = rect(4) * rectfig(4);
        
    %set(0,'units','pixels');
    %set(h(g),'Units','normalized');
            
        axislen = max([height_px,width_px]);
        ticklen = ticklength/axislen*max(rectfig(3:4));
    %end
    
    set(h(g),'TickLength',ticklen*[1 2.5]);    
end