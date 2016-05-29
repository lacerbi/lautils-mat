% LINFUN Evaluate in X a function specified by an ordered vector of 
% x values FX and associated y values FY. 
% The y values are linearly interpolated. If X is outside the range,
% returns NaN.
%
function Y = linfun(X, FX, FY)
    xmin = FX(1);
    xmax = FX(length(FX));
    if size(FX, 1) > 1; FX = [FX; (xmax + 1)]; else FX = [FX (xmax + 1)]; end
    if size(FY, 1) > 1; FY = [FY; 0]; else FY = [FY 0]; end
    
    Y = zeros(size(X, 1), size(X, 2));
    
    for i = 1:length(X)
        if X(i) < xmin || X(i) > xmax; Y(i) = NaN; continue; end       
        xindex = find(FX > X(i), 1) - 1;
        alpha = (X(i) - FX(xindex))/(FX(xindex+1) - FX(xindex));
        Y(i) = FY(xindex) + alpha*(FY(xindex+1)-FY(xindex));        
    end
end