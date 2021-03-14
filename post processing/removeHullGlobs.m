%--------------------------------------------------------------------------
% function to remove small blobs from voxel reconstruction of wing. This
% was originally in 'wingVoxSolidity.m' but it might have more general use
%--------------------------------------------------------------------------
function [hullVox_cleaned, N_chunks, axis_vec] = ...
    removeHullGlobs(hullVox, minVoxNum, padSize, watershedFlag)
% -------------------------
%% params/inputs
if ~exist('minVoxNum','var') || isempty(minVoxNum)
    minVoxNum = 200 ; % size cutoff for connected components. NB: use -1 to just take largest
end
if ~exist('padSize','var') || isempty(padSize)
    padSize = 2 ; % number of voxels to use to pad 3D image space
end
if ~exist('watershedFlag','var') || isempty(watershedFlag)
    watershedFlag = false ; % use watershed transform to remove clumps
end

imMinLvl = 1.5 ; 
% ------------------------------------------
%% generate 3D logical array for voxel image
xmin = min(hullVox(:,1)) - padSize;
xmax = max(hullVox(:,1)) + padSize;
ymin = min(hullVox(:,2)) - padSize;
ymax = max(hullVox(:,2)) + padSize;
zmin = min(hullVox(:,3)) - padSize;
zmax = max(hullVox(:,3)) + padSize;

dx = abs(xmax-xmin)+1;
dy = abs(ymax-ymin)+1;
dz = abs(zmax-zmin)+1;

sizevec = [dx dy dz] ;
% store [xmin, xmax, ymin, ymax, zmin, zmax] in vector to pass along for
% (potential) later reconstruction
axis_vec = [xmin, xmax, ymin, ymax, zmin, zmax] ; 
% ----------------------------
%% BUILD 3D VOLUME FROM HULL
VOL = false(sizevec) ;

for j=1:size(hullVox,1)
    %VOL(hullVox(j,1)+abs(xmin)+1, hullVox(j,2)+abs(ymin)+1,hullVox(j,3)+abs(zmin)+1) = true ; % xxx
    VOL(hullVox(j,1)-xmin+1, hullVox(j,2)-ymin+1, hullVox(j,3)-zmin+1) = true ;
end

% use watershed segmentation to break up volume further?
if watershedFlag
    D = -bwdist(~VOL) ;
    mask = imextendedmin(D, imMinLvl) ;
    D2 = imimposemin(D, mask) ;
    L = watershed(D2) ;
    VOL(L == 0) = 0 ;
end

% remove small voxel globs
CC = bwconncomp(VOL) ; 
blobSizes = cellfun(@(y) length(y), CC.PixelIdxList) ; 
if (minVoxNum < 0)
    [~, good_vol_ind] = max(blobSizes) ;
else
    good_vol_ind = find(blobSizes > minVoxNum) ;
end
N_chunks = length(good_vol_ind) ; % number of voxel blobs above size thresh

% ...but if there are no blobs that are large enough, just return original
% image
if (N_chunks < 1)
    hullVox_cleaned = hullVox ; 
    return
end

% transform back into coordinates
hullVox_cleaned = [] ; 
for q = 1:N_chunks
    ind_curr = good_vol_ind(q) ; 
    [idx1, idx2, idx3] = ind2sub(size(VOL), CC.PixelIdxList{ind_curr}) ;
    tmp = int16([idx1, idx2, idx3]) + ...
        int16(repmat([xmin-1, ymin-1, zmin-1],size(idx1,1),1));
    hullVox_cleaned = [hullVox_cleaned ; tmp] ;
end

end