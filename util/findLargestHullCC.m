function [largestCC, otherCC, CC, idx] = findLargestHullCC (hull, conn, sizevec)
% finds the largest connected component in the given 3D hull (x,y,z)
% and returns it in largestCC.
% otherCC contains all the other voxels that were not included in largestCC
% CC contains the connected components analysis returned by bwconncomp
% idx is the indices of the connected components in CC sorted by descending
% size
%
% sizevec is the size of the 3D volume. Dafault value is [512 512 512] 
% conn is the connectivity used for finding connected components in 3D


% -----------
% DEFINITIONS
% -----------
xmin = min(hull(:,1));
xmax = max(hull(:,1));
ymin = min(hull(:,2));
ymax = max(hull(:,2));
zmin = min(hull(:,3));
zmax = max(hull(:,3));

dx = abs(xmax-xmin)+1;
dy = abs(ymax-ymin)+1;
dz = abs(zmax-zmin)+1;


if (~exist('sizevec','var'))
    sizevec = [dx dy dz] ; %[512 512 512] ;
end

if (~exist('conn','var'))
    conn = 6 ;
end

VOL     = false(sizevec) ;
otherCC = [] ;

% -------------------------
% BUILD 3D VOLUME FROM HULL
% -------------------------
for j=1:size(hull,1) 
    %VOL(hull(j,1)+abs(xmin)+1, hull(j,2)+abs(ymin)+1,hull(j,3)+abs(zmin)+1) = true ; % xxx
    VOL(hull(j,1)-xmin+1, hull(j,2)-ymin+1, hull(j,3)-zmin+1) = true ;
end

% ----------------------------
% ANALYZE CONNECTED COMPONENTS
% ----------------------------
CC = bwconncomp(VOL,conn);

Ncc = length(CC.PixelIdxList) ;
if (Ncc>1)
    svec = zeros(Ncc,1) ; % vector containing the size of all connected components
    for j=1:Ncc
        svec(j) = length(CC.PixelIdxList{j}) ;
    end
    [~, idx] = sort(svec,'descend') ;
    
    [idx1, idx2, idx3] = ind2sub(size(VOL), CC.PixelIdxList{idx(1)}) ;
    largestCC = double([idx1 idx2 idx3 ]) ;
    %largestCC = largestCC-repmat(abs(double([xmin+1, ymin+1, zmin+1])),size(largestCC,1),1) ; % xxx
    largestCC = largestCC + repmat(double([xmin-1, ymin-1, zmin-1]),size(largestCC,1),1) ; 
    for j=2:Ncc
        [idx1, idx2, idx3] = ind2sub(size(VOL), CC.PixelIdxList{idx(j)}) ;
        tmp = int16([idx1, idx2, idx3]) + int16(repmat([xmin-1, ymin-1, zmin-1],size(idx1,1),1)); %tmp = int16([idx1, idx2, idx3]) ;
        otherCC = [ otherCC ; tmp ] ;  %#ok<AGROW>
    end
    
else
    % if there is only one connected component
    largestCC = hull ;
    idx = 1 ;
end



return
