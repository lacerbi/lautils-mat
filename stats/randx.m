function xnew = randx(n,lb,ub,scale,X,nsmpl)
%RANDX Uniformly distributed non-overlapping random point
%   XNEW = RANDX(N,LB,UB) draws N random vectors uniformly distributed in 
%   the region with lower bound LB and upper bound UB. Each new vector
%   tries to minimize overlaps with previously drawn points. Points are 
%   considered as hard spheres with radius 1 in each dimension.
%
%   XNEW = RANDX(N,LB,UB,X,SCALE) assumes vectors are hard ellipsoids with 
%   half axis length SCALE.
%
%   XNEW = RANDX(N,LB,UB,SCALE,X) also avoids clashes with points in the
%   matrix array X. 
%
%   XNEW = RANDX(N,LB,UB,SCALE,X,NSMPL) generates NSMPL candidate random
%   vectors per new point (default NSMPL = 50).

%   Author: Luigi Acerbi <luigi.acerbi@gmail.com>
%   Last version: 17/08/2015

lb = lb(:)'; ub = ub(:)';
nvars = length(lb);

if nargin < 4 || isempty(scale); scale = ones(1,nvars); end
if nargin < 5; X = []; end
if nargin < 6 || isempty(nsmpl); nsmpl = 50; end

xnew = zeros(n,nvars);

% Generate candidate points
for i = 1:n
    if size(X, 1) > 0
        newpoints = bsxfun(@plus,bsxfun(@times,rand(1,nvars,nsmpl), ub-lb),lb);

        % Compute normalized distance between candidates and all existing points
        if size(scale,1) > 1; scale = scale(randi(size(scale,1),[1 1 nsmpl]),:); end    
        dd = squeeze(bsxfunandsum(@minus,@rdivide,@power,newpoints,X,scale,2,2));

        % Compute admissible points (hard sphere with unit radius)
        admis = find(all(dd > 2, 1),1);

        if ~isempty(admis)
            xnew(i,:) = squeeze(newpoints(1,:,admis));
        else
            % No admissible point -- take the one furthest from the closest point
            [mindd,index] = max(min(dd,[],1),[],2);
            xnew(i,:) = squeeze(newpoints(1,:,index));
        end
    else
        xnew(i,:) = rand(1,nvars).*(ub-lb) + lb;
    end
    
    X = [X; xnew(i,:)];
end