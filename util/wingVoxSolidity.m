%--------------------------------------------------------------------------
% function to get solidity (defined as "proportion of the voxels in the
% convex hull that are also in the region") for a given set of wing voxels.
% also return the convex hull, because why not
%
% **NB: the solidity bit of regionprops3 seems to be messed up in some way.
% I'm not sure why, but I often get solidity values greater than one, so we
% need a work-around
%--------------------------------------------------------------------------
function [solidity, convexHull, N_chunks] = wingVoxSolidity(hullVox)
% -----------------------
%% params and boundary definitions
minVoxNum = 400 ;
pad_size = 2 ;

xmin = min(hullVox(:,1)) - pad_size;
xmax = max(hullVox(:,1)) + pad_size;
ymin = min(hullVox(:,2)) - pad_size;
ymax = max(hullVox(:,2)) + pad_size;
zmin = min(hullVox(:,3)) - pad_size;
zmax = max(hullVox(:,3)) + pad_size;

dx = abs(xmax-xmin)+1;
dy = abs(ymax-ymin)+1;
dz = abs(zmax-zmin)+1;

sizevec = [dx dy dz] ;

% ----------------------------
%% BUILD 3D VOLUME FROM HULL
VOL = false(sizevec) ;

for j=1:size(hullVox,1)
    %VOL(hullVox(j,1)+abs(xmin)+1, hullVox(j,2)+abs(ymin)+1,hullVox(j,3)+abs(zmin)+1) = true ; % xxx
    VOL(hullVox(j,1)-xmin+1, hullVox(j,2)-ymin+1, hullVox(j,3)-zmin+1) = true ;
end

% remove small voxel globs
CC = bwconncomp(VOL) ; 
blobSizes = cellfun(@(y) length(y), CC.PixelIdxList) ; 
good_vol_ind = find(blobSizes > minVoxNum) ;
N_chunks = length(good_vol_ind) ; % number of voxel blobs above size thresh

% transform back into coordinates
hullVox_cleaned = [] ; 
for q = 1:N_chunks
    ind_curr = good_vol_ind(q) ; 
    [idx1, idx2, idx3] = ind2sub(size(VOL), CC.PixelIdxList{ind_curr}) ;
    tmp = int16([idx1, idx2, idx3]) + ...
        int16(repmat([xmin-1, ymin-1, zmin-1],size(idx1,1),1));
    hullVox_cleaned = [hullVox_cleaned ; tmp] ;
end
% -----------------------------------
%% find convex hull using alphaShape
shp = alphaShape(double(hullVox_cleaned), Inf) ; 

% get volume from alpha shape and use to calculate solidity 
ch_volume = shp.volume ; 
vox_volume = size(hullVox_cleaned,1) ; 
solidity = vox_volume / ch_volume ; 

% -----------------------------------
%% get coordinates for convex hull
% [qx, qy, qz] = ind2sub(size(VOL), (1:(dx*dy*dz))') ;
% convexHull = int16([qx, qy, qz]) + ...
%         int16(repmat([xmin-1, ymin-1, zmin-1],length(qx),1));
[qx, qy, qz] = meshgrid(xmin:xmax, ymin:ymax, zmin:zmax) ; 
in_shape_idx = inShape(shp, double(qx), double(qy), double(qz)) ; 
convexHull = int16([qx(in_shape_idx), qy(in_shape_idx), qz(in_shape_idx)]) ; 

if (0)
    figure ; hold on
    plot3(hullVox(:,1), hullVox(:,2), hullVox(:,3), 'b.')
    plot3(hullVox_cleaned(:,1), hullVox_cleaned(:,2), hullVox_cleaned(:,3), 'go')
    %plot3(convexHull(:,1), convexHull(:,2), convexHull(:,3), 'rx')
    axis equal
    
    figure ; hold on
     plot3(hullVox_cleaned(:,1), hullVox_cleaned(:,2), hullVox_cleaned(:,3), 'b.')
     plot3(qx(in_shape_idx), qy(in_shape_idx), qz(in_shape_idx), 'rx')
    axis equal
    
    figure ; hold on 
    plot3(hullVox_cleaned(:,1), hullVox_cleaned(:,2), hullVox_cleaned(:,3), 'b.')
    plot(shp)
    axis equal
end

end