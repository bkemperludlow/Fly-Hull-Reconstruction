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
%% clean up wing voxels
[hullVox_cleaned, N_chunks, axis_vec] = removeHullGlobs(hullVox) ; 

% -----------------------------------
%% find convex hull using alphaShape
shp = alphaShape(double(hullVox_cleaned), Inf) ; 

% if the shape is empty, there are very few voxels, so we shouldn't try to
% cluster
if (shp.volume == 0)
    solidity = 1 ; 
    convexHull = hullVox ; 
    return
end

% get volume from alpha shape and use to calculate solidity 
ch_volume = shp.volume ; 
vox_volume = size(hullVox_cleaned,1) ; 
solidity = vox_volume / ch_volume ; 

% -----------------------------------
%% get coordinates for convex hull
% [qx, qy, qz] = ind2sub(size(VOL), (1:(dx*dy*dz))') ;
% convexHull = int16([qx, qy, qz]) + ...
%         int16(repmat([xmin-1, ymin-1, zmin-1],length(qx),1));
[qx, qy, qz] = meshgrid(axis_vec(1):axis_vec(2), axis_vec(3):axis_vec(4),...
    axis_vec(5):axis_vec(6)) ; 
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