function [c, xyunique ] = calcCentroidFromTopView(coords, farPoint, L)
% find the centroid of the voxel set "coords" in two steps:
% 0. include only the points whose distnace to the farPoint is more than L
%    L is typically wingLength/3
% 1. consider only the unique (x,y) pairs, ignoing z coordinate. Find their
%    center of mass on the xy plane (xc, yc)
% 2. consider the voxels (x,y,z) for which (x-xc)^2+(y-yc)^2<=delta2.
%    i.e. take only the distance in the xy plane, ignoring their z.
%    Find the center of mass along z of this set only.
%
% retrun the center of wing as well as the unique coordinates in the xy
% plane. this will be used in the function xx to calculate the projection
% of the chord in the xy plane as seen from the top view.

c = zeros(1,3) ;
delta2 = 4 ;

% step 0
% ------
Nvox = size(coords,1) ;
dst = myNorm ( coords - repmat(farPoint,Nvox,1) ) ;
ind0 = dst > L ;

% step 1
% ------
xyunique = unique(coords(ind0,1:2),'rows') ;
cxy      = mean(xyunique) ;
c(1:2)   = cxy ;

% step 2
% ------
coords2 = coords(ind0,:) ;

ind = (coords2(:,1)-c(1)).^2 + (coords2(:,2)-c(2)).^2 <= delta2 ;

c(3) = mean(coords2(ind,3)) ;

end