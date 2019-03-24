function V = im2voxel(A)
%converts a 3D array of image with values representing intensity to a 2D
%array of voxel coordinates
[r,c,v] = ind2sub(size(A),find(A > 0)); %finds index of coloured pixels
V = [r,c,v];
end