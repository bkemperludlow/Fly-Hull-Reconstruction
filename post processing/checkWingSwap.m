% -------------------------------------------------------------------------
% function that takes in fly analysis structure and checks frames to see if
% the left/right wing labels need to be switched
% -------------------------------------------------------------------------
function [data, swapFlag] = checkWingSwap(data,largePertFlag, swapThresh, ...
    bodyFrameFlag) 
% --------------------
%% inputs and params
if ~exist('largePertFlag','var') || isempty(largePertFlag)
   largePertFlag = false ; % do we need to calculate angles specially
end
if ~exist('swapThresh','var') || isempty(swapThresh)
   swapThresh = 10 ; % units=voxel distance. used to see if swapping changes distances enough 
end
if ~exist('bodyFrameFlag','var') || isempty(bodyFrameFlag)
   bodyFrameFlag = false ;  
end

interpFlag = true ; 
defineConstantsScript % get variable names for different body angles
tol = 1e-12 ; % tolerance for check on how close rotated body axis is to xhat

% ---------------------------------------
%% convert to body frame of reference?
if bodyFrameFlag
    data = labToBodyFrame(data, largePertFlag) ;
end
% ---------------------------------------
%% read in pertinent data from structure
Nframes = data.Nimages ; 
rightWingCM = data.rightWingCM ;
leftWingCM = data.leftWingCM ;
AHat = data.AHat ; 
bodyCM = data.bodyCM ; 

% calculate new angle data, just in case
[rhoTimes, rollVectors] = estimateRollVector(data, largePertFlag, interpFlag) ;
data.rhoTimes = rhoTimes ;
data.rollVectors = rollVectors ;

[anglesLabFrame, ~, ~, ~, ~, ~, ~, ~, ~,rotM_YP, rotM_roll] = ...
    calcAnglesRaw_Sam(data, false, largePertFlag) ;
%bodyRoll = (pi/180)*anglesLabFrame(:,RHO) ; 

% -------------------------------------
%% interpolate wing centers of mass
[rightWingCM_interp, good_idx_R, ~] = interpolateWingCM(data,'right') ;
[leftWingCM_interp, good_idx_L, ~] = interpolateWingCM(data,'left') ;

% this interpolation gives us frames where measured and interpolated values
% differ significantly
bad_idx = ~good_idx_R | ~good_idx_L ;
%bad_frames = find(bad_idx) ;

% -------------------------------------------------------------
%% loop through frames and check if L<-->R swap is appropriate
%
% primary check is based on cross product of two vectors: 1) the vector
% connecting the body center of mass to the wing center of mass in the body
% frame of reference (i.e. unyawed, unpitched, unrolled) and 2) the body
% axis vector in the rotated frame (should be equivalent to the unit vector
% on +x axis if angles are correct). 
%
% secondary check based on distance between "raw" and interpolated wing CMs

swapFlag = false(Nframes,1) ; % initialize flags for various conditions
ignoreFramesFlag = false(Nframes,1) ; 
for ind = 1:Nframes

    % get rotation matrix to rotate body axis to [1; 0; 0]
    unyawpitch = squeeze(rotM_YP(:,:,ind)) ; 
%     unroll = eulerRotationMatrix(0, 0, bodyRoll(ind)) ;
    unroll = squeeze(rotM_roll(:,:,ind)) ; 
    rotM = unroll*unyawpitch ; 
    
    % test that body angles are correctly calculated by rotating body axis
    % and seeing if it lies along positive x axis
    AHat_rot = rotM*AHat(ind,:)' ; 
    if abs(AHat_rot(1) - 1) > tol
        disp('error rotating to body frame. check euler angles')
        keyboard
    end
    
    % rotate wing vectors using rotM and normalize
    rightCMVecRot = rotM*(rightWingCM(ind,:) - bodyCM(ind, :))';
    rightCMVecRotHat = rightCMVecRot./norm(rightCMVecRot) ;
    
    leftCMVecRot = rotM*(leftWingCM(ind,:) - bodyCM(ind, :))';
    leftCMVecRotHat = leftCMVecRot./norm(leftCMVecRot) ;
    
    % check cross product with positive x axis
    rightCross = cross([1 ; 0; 0], rightCMVecRotHat) ;
    leftCross = cross([1 ; 0; 0], leftCMVecRotHat) ;
    
    % this is what we would expect if the wings are labeled correctly (i.e.
    % AHat x R should have negative z component and AHat X L should have
    % positive z component
    rightCrossCheck = (rightCross(3) <= 0) ;
    leftCrossCheck = (leftCross(3) >= 0) ;
    
    % now try to account for different possibilities...
    if (~rightCrossCheck && ~leftCrossCheck) || ...
            (~rightCrossCheck && any(isnan(leftCross))) || ...
            (~leftCrossCheck && any(isnan(rightCross)))
        % in this case, we can be pretty certain we should swap
        swapFlag(ind) = true ;
    elseif (rightCrossCheck && leftCrossCheck) || ...
            (rightCrossCheck && any(isnan(leftCross))) || ...
            (leftCrossCheck && any(isnan(rightCross)))
        % in this case, we can be pretty certain we don't need to swap
        swapFlag(ind) = false;
    elseif (~rightCrossCheck && leftCrossCheck) || ...
            (~leftCrossCheck && rightCrossCheck)
        % this case is trickier. could be due to bad roll estimation or
        % maybe the wings voxel reconstruction is bad. so we'll test the
        % interpolated vs raw values for the wing CMs
        dist_mat = pdist2([rightWingCM(ind,:) ; leftWingCM(ind,:)], ...
        	[rightWingCM_interp(ind,:) ; leftWingCM_interp(ind,:)]) ; 
        
        % for this distance matrix, distances along the diagonal correspond
        % to the (R_raw vs R_interp) and (L_raw vs L_interp) conditions.
        % the off diagonal elements are the swaps. If the sum of off
        % diagonal elements is much smaller, we can be reasonably confident
        % that we should swap
        diagSum = sum(diag(dist_mat)) ; 
        offDiagSum = sum(diag(fliplr(dist_mat))) ;
        distCheck = (offDiagSum - diagSum) < -1*swapThresh ; 
        if distCheck 
            swapFlag(ind) = true ;
        elseif ~distCheck && bad_idx(ind)
            swapFlag(ind) = false ;
            ignoreFramesFlag(ind) = true ; 
        end
    else
        disp('mystery condition. check')
        keyboard
        
    end 
end

% -------------------------------------------------------------
%% swap wings that we flagged and add any new frames to ignore
swapFrames = find(swapFlag) ; 
data = swapWingLeftRight(data,swapFrames) ; 

ignoreFrames = find(ignoreFramesFlag) ; 
if isfield(data, 'ignoreFrames') && ~isempty(ignoreFrames)
    data.ignoreFrames = sort(unique([data.ignoreFrames, ignoreFrames'])) ; 
elseif ~isfield(data, 'ignoreFrames') && ~isempty(ignoreFrames)
    data.ignoreFrames = ignoreFrames' ; 
end

% -------------------------------------------------------------
%% convert back to lab frame, if necessary
if bodyFrameFlag
   data = bodyToLabFrame(data) ;  
end
end