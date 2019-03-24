%--------------------------------------------------------------------------
% function to take a second pass at hull reconstruction/analysis in the
% case of large pitch perturbations, which throw off many of the heuristics
% used to assign vectors to wing/body and calculate angles
%
%   TO DO:
%       -Write the part that updates the data structure with new results
%       -Try to account for frames with merged wings that didn't get picked
%           up during original analysis. UPDATE: did some of this, but need
%           to account for case when one wing is a big blob but the other
%           is fine
%       -Try to interpolate the wing centers of mass as I go? that way the
%           nans wouldn't be such an issue
%--------------------------------------------------------------------------
function [data_new, flagStruct] = largePertCorrection(data)

%% first load in relevant data and make sure that body pitch is correct
% right wing
rightWingCM = data.rightWingCM ;
rightSpanHats = data.rightSpanHats ;
rightChordHats = data.rightChordHats ;
rightChordAltHats = data.chord1AltHats ;
rightWingInd = data.rightWingInd ;
rightWingTips = data.rightWingTips ;
%rightDiags = [(data.diag11Right)' , (data.diag12Right)' ] ;

% left wing
leftWingCM = data.leftWingCM ;
leftSpanHats = data.leftSpanHats ;
leftChordHats = data.leftChordHats ;
leftChordAltHats = data.chord2AltHats ;
leftWingInd = data.leftWingInd ;
leftWingTips = data.leftWingTips ;
%leftDiags = [(data.diag21Left)' , (data.diag22Left)' ] ;

defineConstantsScript
% body info
bodyCM = data.bodyCM ;
AHat = data.AHat ;
[bodyPitch, bodyYaw, rotMats_YP] = calcPitchLargePert(data) ;
bodyRoll = data.anglesLabFrame(:,RHO) ;
%bodyYaw = data.anglesLabFrame(:,PHIB) ;
%bodyYawCheck = data.anglesLabFrame(:,PHIB) ;
N_frames = data.Nimages ;
wingLength = data.wingLength ;

% time info
t = (1/data.params.fps) * ...
    (data.params.startTrackingTime : data.params.endTrackingTime) ;
%tms = 1000 * t ;

% default pitch angle to rotate to
deg2rad = (pi/180) ;
rad2deg = (180/pi) ;
thetaB0 = 45 ;
thetaB0rad = deg2rad*thetaB0 ;

% wing cm velocity threshold for swapping chord
velThreshSwap = 5 ; % vox/frame taken from hullAnalysis

% differences in roll angle between rhoTimes should be checked
rollDiffThresh = 0.9 ; %2.5 ; % radians

% interpolate through wingTips so i don't run into nan values (accuracy not
% super important)
%smoothingParams = setSmoothingParams() ;
rightWingTipsFilt = filter_wing_vecs(rightWingTips,3,2,2000) ;
leftWingTipsFilt = filter_wing_vecs(leftWingTips,3,2,2000) ;

%% initialize data containers
N_vox_R = nan(N_frames,1) ;
N_vox_L = nan(N_frames,1) ;

rightCMRot = nan(N_frames, 3) ;
rightSpanRot = nan(N_frames, 3) ;
rightChordRot = nan(N_frames, 3) ;
rightChordAltRot = nan(N_frames, 3) ;
rightVelRot = nan(N_frames, 3) ;

leftCMRot = nan(N_frames, 3) ;
leftSpanRot = nan(N_frames, 3) ;
leftChordRot = nan(N_frames, 3) ;
leftChordAltRot = nan(N_frames, 3) ;
leftVelRot = nan(N_frames, 3) ;

rotMats = nan(3,3,N_frames) ;

% roll is sort of a weird one, since erroneous points can throw off the
% whole spline. so we'll just assume the first one is okay, and then
% calculate our own roll
rhoTimes = 1 ;
rollVectors = zeros(size(data.rollVectors)) ;
rollCurr = deg2rad*bodyRoll(1) ;
rhoDiff = 0 ; 
%---------------------------------------------
% bunch of flags to store info on the process
%---------------------------------------------
% left right swap?
swapFlag = false(N_frames,1) ;
% had to cluster one of the wings?
clusterRightFlag = false(N_frames,1) ;
clusterLeftFlag = false(N_frames,1) ;
% thought we should cluster but it failed?
clusterError1Flag = false(N_frames,1) ;
clusterError2Flag = false(N_frames,1) ;
% wing center of mass velocity too large?
rightVelErrFlag = false(N_frames,1) ;
leftVelErrFlag = false(N_frames,1) ;
% switched chord and alternate chord?
cSwapFlagR = false(N_frames,1) ;
cSwapFlagL = false(N_frames,1) ;
% inverted chord?
cInvFlagR = false(N_frames,1) ;
cInvFlagL = false(N_frames,1) ;
% did we calculate a new rho time?
newRhoFlag = false(N_frames, 1) ;
% struct to store these in later
flagStruct = struct ;
%% get the indices and wing voxel size for every frame
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

for j = 1:N_frames
    row1 = frameStartInd(j) ;
    row2 = frameEndInd(j) ;
    
    %coords = data.res(row1:row2,2:4) ; % xyz positions of each voxel in frame=1st col value
    IDX = data.RESIDX(row1:row2,:) ; % Classifying each voxel as L/R wing or body
    wingRows_R = (IDX(:,rightWingInd)==1) ;
    wingRows_L = (IDX(:,leftWingInd)==1) ;
    N_vox_R(j) = sum(wingRows_R) ;
    N_vox_L(j) = sum(wingRows_L) ;
end

mean_vox_R = nanmean(N_vox_R) ;
std_vox_R = nanstd(N_vox_R) ;
mean_vox_L = nanmean(N_vox_L) ;
std_vox_L = nanstd(N_vox_L) ;

N_vox_R_zScore = (N_vox_R - mean_vox_R)./std_vox_R ;
N_vox_L_zScore = (N_vox_L - mean_vox_L)./std_vox_L ;

%% try to determine a threshold for wing center of mass velocity
% i.e. if the wing cm is moving faster than this, it's likely an error
t0_ind = find(t==0) ;
rightWingCM2Body = rightWingCM - bodyCM ;
leftWingCM2Body = leftWingCM - bodyCM ;

rightWingCMVel = myNorm(diff(rightWingCM2Body(1:t0_ind,:))) ;
leftWingCMVel = myNorm(diff(leftWingCM2Body(1:t0_ind,:))) ;

if (0)
    figure ;
    hold on
    histogram(rightWingCMVel,100)
    histogram(leftWingCMVel,100)
end

velThresh = nanmean([rightWingCMVel ; leftWingCMVel]) + ...
    2*nanstd([rightWingCMVel ; leftWingCMVel]) ;

%% now check through the frames and assess failure cases, including:
% missing voxels for one or both wings
% left and right wings are switched
% chords are swapped around
% ...

% these are the frames where the wing size is either conspicuously small or
% large
clust_R_idx = ((N_vox_R_zScore - N_vox_L_zScore) > 2.5) & ...
    (N_vox_R_zScore > 0) ;
clust_L_idx = ((N_vox_R_zScore - N_vox_L_zScore) < -2.5) & ...
    (N_vox_L_zScore > 0) ;
clust_idx = clust_R_idx | clust_L_idx ;

% these frames lack a span, chord, or cm estimate for one of the wings
R_nan_idx = any(isnan([rightWingCM, rightSpanHats, rightChordHats]),2) ;
L_nan_idx = any(isnan([leftWingCM, leftSpanHats, leftChordHats]),2) ;

% for finding frames where left and right wing are swapped, loop through
% and try to un-yaw then un-pitch
for i = 1:N_frames
    % count off frames
    fprintf('Working on frame %d of %d ... \n',i,N_frames)
    % first check to see if these are frames we can work with
    if (R_nan_idx(i) && L_nan_idx(i)) || (R_nan_idx(i) && ~clust_L_idx(i)) ...
            || (L_nan_idx(i) && ~clust_R_idx(i))
        % should add in extra conditions here to test for merged wings
        continue
    else
        %% find rotation matrix to transform body to consistent frame
        
        % rotation matrix to take the body axis to z axis (and try to unroll)
        M_YP = squeeze(rotMats_YP(:,:,i)) ; % yaw and pitch rotation
        % note that for roll, we're using the most recent estimate of the
        % roll vector
        M_R = eulerRotationMatrix(0,0,rollCurr ) ; % roll rotation
        M1 = M_R*M_YP ; %to strict body axis
        % rotation matrix to pitch down fly by 45 degrees
        M2 = eulerRotationMatrix(0, -thetaB0rad, 0) ; % pitch down by thetaB w.r.t body axis
        M = M2 * M1 ;
        rotMats(:,:,i) = M ;
        
        % rotate long body axis
        AHatRot = M*AHat(i,:)' ;
        % if we don't align along x axis, then we
        
        %% see if it's necessary to cluster one of the wings
        %----------------------------------------------------
        if clust_idx(i)
            % determine the side that likely needs clustering
            if clust_R_idx(i)
                wing_str = 'right' ;
            elseif clust_L_idx(i)
                wing_str = 'left' ;
            else
                disp('this should not happen')
                keyboard
            end
            % perform clustering
            [wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
                clusterWings(data, i, wing_str) ;
            
            % if we didn't get a successful clustering, and we're missing
            % information on one wing, skip this frame
            if badClusterFlag && (R_nan_idx(i) || L_nan_idx(i))
                clusterError1Flag(i) = true ;
                continue
            elseif badClusterFlag && ~(R_nan_idx(i) || L_nan_idx(i))
                clusterError2Flag(i) = true ;
            else
                % if we did successfully cluster, we now need to figure out
                % which is left vs right. this should be done most easily in
                % the rotated frame
                
                centroid1 = centroids(1,:) ;
                centroid2 = centroids(2,:) ;
                rotCentroid1 = M*(centroid1 - bodyCM(i,:))' ;
                rotCentroid2 = M*(centroid2 - bodyCM(i,:))' ;
                
                % take cross product with rotated body axis to check L/R
                CheckVec1 = cross(AHatRot, rotCentroid1) ;
                check1 = sign(CheckVec1(3)) ;
                CheckVec2 = cross(AHatRot, rotCentroid2) ;
                check2 = sign(CheckVec2(3)) ;
                
                if (check1 < 0) && (check2 > 0)
                    % in this case, we have unambiguously that 1=R and 2=L
                    right_idx = 1 ;
                    left_idx = 2 ;
                elseif (check1 > 0) && (check2 < 0)
                    % in this case, we have unambiguously that 1=L and 2=R
                    right_idx = 2 ;
                    left_idx = 1 ;
                else
                    % in this case, we'll check by just finding one that points
                    % along POSITIVE y axis the most, and call that left
                    rotCentroids = [rotCentroid1' ; rotCentroid2'] ;
                    yHatMat = repmat([0, 1, 0],2,1) ;
                    dotWithY = dot(rotCentroids, yHatMat, 2)./myNorm(rotCentroids) ;
                    [~, left_idx] = max(dotWithY) ;
                    [~, right_idx] = min(dotWithY) ;
                end
                
                wingVoxR = wingVox(label_idx == right_idx,:) ;
                rightCM = centroids(right_idx,:) ;
                wingRowsR = wingRows(:,right_idx) ;
                wingVoxL = wingVox(label_idx == left_idx,:) ;
                leftCM = centroids(left_idx,:) ;
                wingRowsL = wingRows(:,left_idx) ;
                
                % now calculate wing vectors and store them in data struct so
                % we can proceed with other checks
                rightRefVecs = [bodyCM(i,:); bodyCM(i-1,:); ...
                    rightWingTipsFilt(i-1,:)] ;
                [spanHatR, chordHatR, chordAltHatR, ~] = ...
                    estimate_wing_vecs(wingVoxR, rightRefVecs, wingLength, ...
                    [], [], rightCM) ;
                
                leftRefVecs = [bodyCM(i,:); bodyCM(i-1,:); ...
                    leftWingTipsFilt(i-1,:)] ;
                [spanHatL, chordHatL, chordAltHatL, ~] = ...
                    estimate_wing_vecs(wingVoxL, leftRefVecs, wingLength, ...
                    [], [], leftCM) ;
                
                % add to storage arrays
                rightWingCM(i,:) = rightCM ;
                rightSpanHats(i,:) = spanHatR ;
                rightChordHats(i,:) = chordHatR ;
                rightChordAltHats(i,:) = chordAltHatR ;
                
                leftWingCM(i,:) = leftCM ;
                leftSpanHats(i,:) = spanHatL ;
                leftChordHats(i,:) = chordHatL ;
                leftChordAltHats(i,:) = chordAltHatL ;
                
                %...including fucking voxels, blurg
                row1 = frameStartInd(i) ;
                row2 = frameEndInd(i) ;
                data.RESIDX(row1:row2,2:3) = [wingRowsR, wingRowsL] ;
                
                % make sure to store whether or not we clustered
                if strcmp(wing_str,'right')
                    clusterRightFlag(i) = ~badClusterFlag ;
                elseif strcmp(wing_str,'left')
                    clusterLeftFlag(i) = ~badClusterFlag ;
                end
                
            end
        end
        
        %% now rotate wing and body vectors:
        %-------------------------------------
        % right wing
        [rotCR, rotSpanR, rotChordR, rotChordAltR] = ...
            rotateWingVecs(rightWingCM(i,:), rightSpanHats(i,:), ...
            rightChordHats(i,:), rightChordAltHats(i,:), bodyCM(i,:), M) ;
        % left wing
        [rotCL, rotSpanL, rotChordL, rotChordAltL] = ...
            rotateWingVecs(leftWingCM(i,:), leftSpanHats(i,:), ...
            leftChordHats(i,:), leftChordAltHats(i,:), bodyCM(i,:), M) ;
        
        %----------------------------------------------------------------------
        % plot check to make sure that rotations are working correctly
        if (0)
            % hat vectors
            figure ;
            hold on
            plot3([0, AHatRot(1)], [0, AHatRot(2)], [0, AHatRot(3)],...
                'k-','linewidth',3)
            %plot3([0, AHat(i,1)], [0, AHat(i,2)], [0, AHat(i,3)],'k--','linewidth',2)
            % spans
            plot3(AHatRot(1)/2 + [0, rotSpanR(1)], AHatRot(2)/2 + ...
                [0, rotSpanR(2)], AHatRot(3)/2 + [0, rotSpanR(3)],...
                'r-','linewidth',3)
            plot3(AHatRot(1)/2 + [0, rotSpanL(1)], AHatRot(2)/2 + ...
                [0, rotSpanL(2)], AHatRot(3)/2 + [0, rotSpanL(3)],....
                'b-','linewidth',3)
            % chords
            plot3(AHatRot(1)/2 + rotSpanR(1)/2 + [0, rotChordR(1)]/2, ...
                AHatRot(2)/2 + rotSpanR(2)/2 + [0, rotChordR(2)]/2,...
                AHatRot(3)/2 + rotSpanR(3)/2 + [0, rotChordR(3)]/2,...
                'b--','linewidth',3)
            plot3(AHatRot(1)/2 + rotSpanL(1)/2 + [0, rotChordL(1)]/2, ...
                AHatRot(2)/2 + rotSpanL(2)/2 + [0, rotChordL(2)]/2,...
                AHatRot(3)/2 + rotSpanL(3)/2 + [0, rotChordL(3)]/2,...
                'r--','linewidth',3)
            % alt chords
            plot3(AHatRot(1)/2 + rotSpanR(1)/2 + [0, rotChordAltR(1)]/2, ...
                AHatRot(2)/2 + rotSpanR(2)/2 + [0, rotChordAltR(2)]/2,...
                AHatRot(3)/2 + rotSpanR(3)/2 + [0, rotChordAltR(3)]/2,...
                'b:','linewidth',3)
            plot3(AHatRot(1)/2 + rotSpanL(1)/2 + [0, rotChordAltL(1)]/2, ...
                AHatRot(2)/2 + rotSpanL(2)/2 + [0, rotChordAltL(2)]/2,...
                AHatRot(3)/2 + rotSpanL(3)/2 + [0, rotChordAltL(3)]/2,...
                'r:','linewidth',3)
            
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            grid on
            box on
            axis equal
            title(['Frame ' num2str(i,'%04d')])
            % wing centers of mass
            %             figure ; hold on
            %             plot3([0,rotCR(1)],[0,rotCR(2)],[0,rotCR(3)],'r-','linewidth',3)
            %             plot3([0,rotCL(1)],[0,rotCL(2)],[0,rotCL(3)],'b-','linewidth',3)
            %             axis equal
            %             grid on
            %             box on
            %             xlabel('X')
            %             ylabel('Y')
            %             zlabel('Z')
        end
        
        %----------------------------------------------------------------------
        %% perform heuristic checks on rotated fly
        %---------------------------------------
        % first check if we need to swap wings
        rightCheckVec = cross(AHatRot, rotCR) ;
        rightCheck = sign(rightCheckVec(3)) ;
        
        leftCheckVec = cross(AHatRot, rotCL) ;
        leftCheck = sign(leftCheckVec(3)) ;
        
        % if true, swap 'em
        if (rightCheck > 0) && (leftCheck < 0)
            swapFlag(i) = true ;
            [rotCR, rotSpanR, rotChordR, rotChordAltR, rotCL, rotSpanL,...
                rotChordL, rotChordAltL] = swapRotVecs(rotCR, rotSpanR, rotChordR, ...
                rotChordAltR, rotCL, rotSpanL, rotChordL, rotChordAltL) ;
        end
        
        %---------------------------------------
        % next, calculate wing velocity, if possible. use to see if chords
        % need swapping and/or if cm is wayyyy off
        
        % first try to get an estimate of the wing's previous position,
        % either by reading it from array or interpolating
        if (i >=2)
            %{
            % right wing center of mass from previous frame
            if ~any(isnan(rightCMRot(i-1,:)))
                rightCMRot_prev = rightCMRot(i-1,:) ;
            elseif any(isnan(rightCMRot(i-1,:))) && (i > 4)
                rightCMRot_prev = interp_3D_vec(rightCMRot, i-1) ;
            else
                rightCMRot_prev = nan(1,3) ;
            end
            % left wing center of mass from previous frame
            if any(isnan(leftCMRot(i-1,:)))
                leftCMRot_prev = leftCMRot(i-1,:) ;
            elseif any(isnan(leftCMRot(i-1,:))) && (i > 4)
                leftCMRot_prev = interp_3D_vec(leftCMRot, i-1) ;
            else
                leftCMRot_prev = nan(1,3) ;
            end
            %}
            rightCMRot_prev = rightCMRot(i-1,:) ;
            leftCMRot_prev = leftCMRot(i-1,:) ;
        end
        % the reason for the dumb-looking double if statement is that i
        % can't think of a clever way to deal with the fact that sometimes
        % we won't be able to interp, and i don't want to just throw a
        % 'continue'
        if (i >= 2) && ~any(isnan(rightCMRot_prev)) && ...
                ~any(isnan(leftCMRot_prev))
            % then perform simple velocity estimation
            vRotR = rotCR - rightCMRot_prev' ;
            vRotL = rotCL - leftCMRot_prev' ;
            
            vRotR_nrm = norm(vRotR) ;
            vRotL_nrm = norm(vRotL) ;
            %-------------------------
            % first check magnitude
            
            % if both velocities are too large, it's likely they
            % switched--not really that accurate actually
            %{
            if (vRotR_nrm > velThresh) && (vRotL_nrm > velThresh) && ~swapFlag(i)
                swapFlag(i) = true ;
                [rotCR, rotSpanR, rotChordR, rotChordAltR, rotCL, rotSpanL,...
                    rotChordL, rotChordAltL] = swapRotVecs(rotCR, rotSpanR, rotChordR, ...
                    rotChordAltR, rotCL, rotSpanL, rotChordL, rotChordAltL) ;
                vRotR = rotCR - rightCMRot(i-1,:)' ;
                vRotL = rotCL - leftCMRot(i-1,:)' ;
                vRotR_nrm = norm(vRotR) ;
                vRotL_nrm = norm(vRotL) ;
            end
            %}
            % but if they're both large and we already swapped, probably
            % just a bad frame
            rightVelErrFlag(i) = (vRotR_nrm > velThresh) ;
            leftVelErrFlag(i) = (vRotL_nrm > velThresh) ;
            
            %-------------------------
            % then compare direction of chord and velocity
            if vRotR_nrm >= velThreshSwap
                [rotChordR, rotChordAltR, cSwapFlagR(i), cInvFlagR(i)] = ...
                    checkVelAndChord(rotChordR, rotChordAltR, vRotR,...
                    vRotR_nrm, velThreshSwap, rotSpanR) ;
            end
            if vRotL_nrm >= velThreshSwap
                [rotChordL, rotChordAltL, cSwapFlagL(i), cInvFlagL(i)] = ...
                    checkVelAndChord(rotChordL, rotChordAltL, vRotL,...
                    vRotL_nrm, velThreshSwap, rotSpanL) ;
            end
            
            %-------------------------
            % finally, use vel to check for rhoTimes. we want to make sure
            % that the spans are pointed out to the left/right, and that
            % the wings are moving backward
            %{
            if (vRotR_nrm ~= 0) && (vRotL_nrm ~= 0)
                spanRollCheck = (dot(rotSpanR, rotSpanL) < -0.90) ;
                %velRollCheck = (abs(dot(rotSpanR, vRotR)/vRotR_nrm) < 0.1)...
                %    & (abs(dot(rotSpanL, vRotL)/vRotL_nrm) < 0.1) ;
                if spanRollCheck %&& velRollCheck
                    rhoTimes = sort([rhoTimes , i]) ;
                    rollVecRot = rotCL - rotCR ;
                    rollVecRot = rollVecRot./norm(rollVecRot) ;
                    rollVec = M'*rollVecRot ;
                    rollVec = rollVec' - AHat(i,:)*dot(rollVec, AHat(i,:)) ;
                    rollVec = rollVec./norm(rollVec) ;
                    
                    rollVectors(i,:) = rollVec ;
                    [bodyRoll_new, ~, ~, rhoSamp] = calcBodyRoll(rhoTimes, ...
                        rollVectors, t, rotMats_YP, data.params) ;
                    rollTemp = deg2rad*rhoSamp(end) ;
                    if abs(rollTemp-rollCurr) < rollDiffThresh
                        rollCurr = rollTemp ;
                        newRhoFlag(i) = true ;
                    else
                        rollVectors(i,:) = 0 ;
                        rhoTimes = rhoTimes(1:end-1) ;
                        disp('new roll vector is wildly different--not using')
                    end
                    if (0)
                        figure ; hold on
                        plot(bodyRoll_new)
                        plot(rhoTimes, bodyRoll_new(rhoTimes),'ko')
                    end
                end
            end
            %}
            % make sure to store velocity
            rightVelRot(i,:) = vRotR ;
            leftVelRot(i,:) = vRotL ;
        end
        
        %% check for body roll angle (this time w/o velocity measurement)
        spanRollCheck = (dot(rotSpanR, rotSpanL) < -0.85) ;
        %velRollCheck = (abs(dot(rotSpanR, vRotR)/vRotR_nrm) < 0.1)...
        %    & (abs(dot(rotSpanL, vRotL)/vRotL_nrm) < 0.1) ;
        if spanRollCheck && ~ismember(i,rhoTimes) 
            rhoTimes = sort([rhoTimes , i]) ;
            rollVecRot = rotCL - rotCR ;
            rollVecRot = rollVecRot./norm(rollVecRot) ;
            % this function should find the change in roll angle from
            % previous estimation
            [deltaRoll, rollVecOut] = ...
                calcRollRotatedFrame(rollVecRot, AHatRot, M2) ; 
            % store the lab frame roll vector for future reference
            rollVec = M1'*rollVecOut ;
            rollVec = rollVec' - AHat(i,:)*dot(rollVec, AHat(i,:)) ;
            rollVec = rollVec./norm(rollVec) ;
            rollVectors(i,:) = rollVec ;
            
%             [bodyRoll_new, ~, ~, rhoSamp] = calcBodyRoll(rhoTimes, ...
%                 rollVectors, t, rotMats_YP, data.params) ;
%             rollTemp = deg2rad*rhoSamp(end) ;
            if abs(deltaRoll) < rollDiffThresh
                rollCurr = rollCurr + deltaRoll ;
                rhoDiff = [rhoDiff, deltaRoll] ; 
                newRhoFlag(i) = true ;
            else
                rollVectors(i,:) = 0 ;
                rhoTimes = rhoTimes(1:end-1) ;
                disp('new roll estimate is wildly different--not using')
            end
            if (0)
                bodyRoll_new = rad2deg*cumsum(rhoDiff) + bodyRoll(1) ; 
                frames = 1:i ; 
                bodyRollInterp = interp1(rhoTimes, bodyRoll_new, frames,...
                    'spline') ; 
                figure ; hold on
                plot(frames, bodyRollInterp,'r-')
                plot(rhoTimes, bodyRoll_new,'ko')
                
            end
        end
        
        %% quick sign check on chord z component (should be positive)
        if rotChordR(3) < 0
            rotChordR = - rotChordR ;
            cInvFlagR(i) = true ;
        end
        if rotChordL(3) < 0
            rotChordL = - rotChordL ;
            cInvFlagL(i) = true ;
        end
        %% store new wing vectors
        
        % right wing
        rightCMRot(i,:) = rotCR ;
        rightSpanRot(i,:) = rotSpanR ;
        rightChordRot(i,:) = rotChordR ;
        rightChordAltRot(i,:) = rotChordAltR ;
        
        %left wing
        leftCMRot(i,:) = rotCL ;
        leftSpanRot(i,:) = rotSpanL ;
        leftChordRot(i,:) = rotChordL ;
        leftChordAltRot(i,:) = rotChordAltL ;
        
        
    end
end

%% calculate wing angles based on the vectors from the transformed frame


%% update data structure and take care of wing L<->R swapping
disp('Updating data structure...')
data_new = data ;
%bad_frames = R_nan_idx | L_nan_idx | rightVelErrFlag | leftVelErrFlag ;
bad_frames = R_nan_idx | L_nan_idx ;
for k = 1:N_frames
    % perform the necessary swaps (the rotated data has been swapped
    % already, but this does some of the other stuff e.g. wing tips
    if swapFlag(k)
        data_new = swapWingLeftRight(data_new, k) ;
    end
    
    % update wing vectors by inverting rotation and translating
    if ~bad_frames(k)
        M_curr = squeeze(rotMats(:,:,k)) ;
        Minv = M_curr' ;
        % right wing
        [data_new.rightWingCM(k,:), data_new.rightSpanHats(k,:), ...
            data_new.rightChordHats(k,:), data_new.chord1AltHats(k,:)] = ...
            unrotateWingVecs(rightCMRot(k,:), rightSpanRot(k,:), ...
            rightChordRot(k,:), rightChordAltRot(k,:), bodyCM(k,:), Minv) ;
        % left wing
        [data_new.leftWingCM(k,:), data_new.leftSpanHats(k,:), ...
            data_new.leftChordHats(k,:), data_new.chord2AltHats(k,:)] = ...
            unrotateWingVecs(leftCMRot(k,:), leftSpanRot(k,:), ...
            leftChordRot(k,:), leftChordAltRot(k,:), bodyCM(k,:), Minv) ;
    end
end

% add in other fields that we've modified
% data_new.rhoTimes = rhoTimes ;
% data_new.rollVectors = rollVectors ;
data_new.rhoTimes = find(newRhoFlag)' ;
rollVectors_new = rollVectors ;
rollVectors_new(~newRhoFlag,:) = 0 ;
data_new.rollVectors = rollVectors_new ;
data_new.newRhoSamp = rad2deg*cumsum(rhoDiff) + bodyRoll(1) ; 

% make sure to ignore bad frames later
ignoreFramesNew = find(bad_frames) ;
if  isfield(data_new,'ignoreFrames')
    data_new.ignoreFrames = sort([data_new.ignoreFrames, ignoreFramesNew']) ;
else
    data_new.ignoreFrames = ignoreFramesNew' ;
end

% save the various flags as diagnostics for later

flagStruct.R_nan_idx = R_nan_idx ;
flagStruct.L_nan_idx = L_nan_idx ;
flagStruct.swapFlag = swapFlag ;
flagStruct.clusterRightFlag = clusterRightFlag ;
flagStruct.clusterLeftFlag = clusterLeftFlag ;
flagStruct.clusterError1Flag = clusterError1Flag ;
flagStruct.clusterError2Flag = clusterError2Flag ;
flagStruct.rightVelErrFlag = rightVelErrFlag ;
flagStruct.leftVelErrFlag = leftVelErrFlag ;
flagStruct.cSwapFlagR = cSwapFlagR ;
flagStruct.cSwapFlagL = cSwapFlagL ;
flagStruct.cInvFlagR = cInvFlagR ;
flagStruct.cInvFlagL = cInvFlagL ;
flagStruct.newRhoFlag = newRhoFlag ;

disp('Done updating data structure')

%% calculate wing/body angles since we have the coordinate transforms
% assigne empty arrays for angle data
[anglesLabFrame, anglesBodyFrame, ~, ~, ~, ~, ~, ~, ~] = ...
    calcAnglesRaw_Sam(data_new, false ,true) ; 


% these arrays in data structure 
data_new.anglesLabFrame = anglesLabFrame ; 
data_new.anglesBodyFrame = anglesBodyFrame ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------
% perfrom rotation
function [rotCM, rotSpan, rotChord, rotChordAlt] = ...
    rotateWingVecs(wingCM, wingSpan, wingChord, wingChordAlt, bodyCM, rotM)

% wing CM
rotCM = rotM * (wingCM - bodyCM)' ;

% wing chord
rotTempChord   = rotM * (wingCM + wingChord - bodyCM)' ;
rotChord = rotTempChord - rotCM ;

% wing chord alt
rotTempChordAlt   = rotM * (wingCM + wingChordAlt - bodyCM)' ;
rotChordAlt = rotTempChordAlt - rotCM ;

% wing span
rotTempSpan   = rotM * ( wingCM + wingSpan - bodyCM)' ;
rotSpan = rotTempSpan - rotCM ;

end
%---------------------------
% undo rotation
function [wingCM, wingSpan, wingChord, wingChordAlt] = ...
    unrotateWingVecs(rotCM, rotSpan, rotChord, rotChordAlt, bodyCM, rotMinv)

%delta = 1000*eps ;
% make sure everybody is a column vector
if size(rotCM,2) > size(rotCM,1)
    rotCM = rotCM' ;
end
if size(rotSpan,2) > size(rotSpan,1)
    rotSpan = rotSpan' ;
end
if size(rotChord,2) > size(rotChord,1)
    rotChord = rotChord' ;
end
if size(rotChordAlt,2) > size(rotChordAlt,1)
    rotChordAlt = rotChordAlt' ;
end
if size(bodyCM,2) > size(bodyCM,1)
    bodyCM = bodyCM' ;
end

% wing CM
wingCM = rotMinv * rotCM + bodyCM ;

% wing chord
wingChord = (rotMinv * rotChord)  ;

% wing chord alt
wingChordAlt = (rotMinv * rotChordAlt) ;

% wing span
wingSpan = (rotMinv * rotSpan) ;

% make sure wingChord is perpendicular to wingSpan
wingChord = wingChord - dot(wingChord, wingSpan).*wingSpan ;
wingChord = wingChord ./ norm(wingChord) ;

end

%-------------------------------------
% swap left and right rotated vectors
function [rotCR, rotSpanR, rotChordR, rotChordAltR, rotCL, rotSpanL,...
    rotChordL, rotChordAltL]  =  swapRotVecs(rotCR, rotSpanR, rotChordR,...
    rotChordAltR, rotCL, rotSpanL, rotChordL, rotChordAltL)

tempCM = rotCR ;
tempSpan = rotSpanR ;
tempChord = rotChordR ;
tempChordAlt = rotChordAltR ;

rotCR = rotCL ;
rotSpanR = rotSpanL ;
rotChordR = rotChordL ;
rotChordAltR = rotChordAltL ;

rotCL = tempCM ;
rotSpanL = tempSpan ;
rotChordL = tempChord ;
rotChordAltL = tempChordAlt ;
end

%-------------------------------------
% check chord against velocity vec
function [chord, chordAlt, chordSwapFlag, chordInvertFlag] = ...
    checkVelAndChord(chord, chordAlt, vWing, vWing_nrm, velThresh, span)

% normalize velocity vector
vWingHat = vWing./vWing_nrm ;

% project out component that is parallel to span
vWingHat = vWingHat - span * dot(span, vWingHat) ;
vWingHat = vWingHat./norm(vWingHat) ;

% now compare with chord and alternative chord
dot1 = dot(chord, vWingHat) ;
dot2 = dot(chordAlt, vWingHat) ;

% determine whether wing velocity suggests swapping
if (dot2>dot1) % swap
    velSwapFlag = true ;
else
    velSwapFlag = false ;
end
% if neither chord aligns well with velocity, then just invert current
% chord
if (dot1<0) && (dot2<0) && (vWing_nrm >= velThresh) && ~velSwapFlag
    chordInvertFlag = true ;
    chord = -1*chord ;
else
    chordInvertFlag = false ;
end

% if the other chord aligns better with other velocity AND the velocity is
% sufficiently large, swap chord and chordAlt
chordSwapFlag = (velSwapFlag && vWing_nrm >= velThresh) ;

if chordSwapFlag
    tmp = chord ;
    chord = chordAlt ;
    chordAlt = tmp ;
end

end

%----------------------------------------------------------------
% simple function to interpolate through nans for e.g. wing cm
% **NB: only returns a single vector currently
function out_vec = interp_3D_vec(in_vecs, ind)

% check that "in_vecs" has 3D vectors as rows
if size(in_vecs,2) > size(in_vecs,1)
    in_vecs = in_vecs' ;
end

N_frames = size(in_vecs,1) ;
N_dim = size(in_vecs,2) ; % should be 3 for all of our uses

% assign interp values. we'll also do a hampel filter to try to ignore
% outliers
frames = 1:N_frames ;
out_vec = zeros(1,N_dim) ;
for i = 1:N_dim
    nan_ind = isnan(in_vecs(:,i)) ;
    [~, hampel_ind] = hampel(in_vecs(:,i)) ;
    ignore_ind = nan_ind | hampel_ind ;
    out_vec(i) = interp1(frames(~ignore_ind), in_vecs(~ignore_ind,i), ind,...
        'spline','extrap') ;
end

if (0)
    figure
    for j = 1:N_dim
        subplot(N_dim,1,j)
        hold on
        plot(frames, in_vecs(:,i), 'b-')
        plot(ind, out_vec(i),'rx')
        ylabel(num2str(j))
    end
    disp('plot check')
end

end

%----------------------------------------------------------------
% function to calculate roll in the rotated frame
function [rho_curr, rollVecOut] = ...
    calcRollRotatedFrame(rollVecRot, AHatRot, rotM_B0)

% first undo the rotation that pitches the fly up from the x axis
rollVecRot2 = rotM_B0' * rollVecRot ; 
AHatRot2 = rotM_B0' * AHatRot ; 

%now make sure to project out any component of the roll vector in the body
%axis direction
rollVecRot2 = rollVecRot2 - dot(rollVecRot2,AHatRot2).*AHatRot2 ; 
rollVecRot2 = rollVecRot2 ./ norm(rollVecRot2) ; 

% now find roll angle 
rho_curr = acos( rollVecRot2(2) ) * sign(rollVecRot2(3)) ; % RADIANS!
rollVecOut = rollVecRot2 ; 
end