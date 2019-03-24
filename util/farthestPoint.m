function [p, dst, list] = farthestPoint(coords, p0, L)
% finds the point in coords which is farthest from p0
% returns the point coordinate p and its distance from p0 in dst.
% if there are several points with the same distance, return only one.
% also finds the indices of voxels in coords that whose distance from p is
% smaller or equal to L

% find p
Nvox = size(coords,1) ;
mat1 = double(coords) - double(repmat(p0, Nvox,1)) ;
dst2vec  = sum (mat1 .* mat1, 2) ;
[dst, ind] = max(dst2vec) ;
p = coords(ind,:) ;

% find list
mat1 = double(coords) - double(repmat(p, Nvox, 1)) ;
dst2vec  = sum (mat1 .* mat1, 2) ;
list = (dst2vec <= L^2) ;

end