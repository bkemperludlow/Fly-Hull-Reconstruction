%--------------------------------------------------------------------------
% Function to determine relevant wing voxels to consider by only looking
% at ones that are a certain distance away from the body's center of mass.
% works on the raw output of hullReconstruction
%--------------------------------------------------------------------------
function [wingLargestCC, farPoint, list, farPointDist] = ...
    refineWingVoxels(wingCoords, bodyCM, lengthThresh)

% get largest connected component within voxel coordinates
wingLargestCC_init = findLargestHullCC (wingCoords);

% from the largest CC, find the farpoint and then select the voxels
% that are close to the farpoint.
[farPoint, farPointDist, list] = farthestPoint(wingLargestCC_init, ...
    bodyCM,lengthThresh) ;


% find the largest connected component within the voxels in list
% list are the indices of the voxles in wing1LargestCC_init whose distance
% from "farPoint" is smaller than LL.
%
% largest component is calculated AGAIN, to account for a case where cutting
% the LL distance left us with two connected components
wingLargestCC = findLargestHullCC (wingLargestCC_init(list,:));

clear wingLargestCC_init

% make sure the far point is included in the connected components
% if not, take another one
if (~ismember(farPoint, wingLargestCC,'rows'))
    tmp = setdiff( wingCoords(list,:), wingLargestCC, 'rows') ;
    wingLargestCC = findLargestHullCC(tmp) ;
    clear tmp
end

% check again
if (~ismember(farPoint, wingLargestCC,'rows'))
    disp('Farthest point not part of wing voxels...');
    %keyboard ;
end


end