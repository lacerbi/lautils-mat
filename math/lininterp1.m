function Vout = lininterp1(X,V,Xq)

if isvector(X); X = X(:); end
if isvector(V); V = V(:); end

% Gets the x spacing
dx = 1./(X(2,:)-X(1,:));  % one over to perform divide only once
Xq = bsxfun(@minus,Xq,X(1,:));      % subtract minimum of x

Xqi = bsxfun(@times,Xq,dx)+1;          % indices of nearest-lower-neighbors
flag = Xqi<1 | Xqi>size(X,1) | isnan(Xqi);
                                % finds indices out of bounds
Xqi = floor(Xqi);
                                
V = [V; NaN(1, size(V,2))];
Xqi(flag) = size(V,1)-1;

delta = Xqi - bsxfun(@times,Xq,dx);

linind1 = bsxfun(@plus, Xqi, size(V,1)*(0:size(V,2)-1));
linind2 = bsxfun(@plus, Xqi + 1, size(V,1)*(0:size(V,2)-1));

out1 = reshape(V(linind1),size(Xq));
out2 = reshape(V(linind2),size(Xq));

out1(isnan(out1)) = 0;
out2(isnan(out2)) = 0;

Vout = bsxfun(@times,delta,out1) + bsxfun(@times,1-delta,out2);
Vout(flag) = NaN;

end