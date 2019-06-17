%--------------------------------------------------------------------------
% In cases when the fly undergoes large pitch perturbations, use a
% quaternion conversion to avoid gimbal lock.
%
% Version _mk2 tries to deal with the fact that, when solving for the
% quaternions that rotate frame to frame, you get some rotation about the
% roll axis, which we want to avoid. we want to try to constrain that
%--------------------------------------------------------------------------

function [bodyPitch, bodyYaw, rotM_array] = ...
    calcPitchLargePert(data, debugFlag)
%--------------------------------------------------------------------------
% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ;
end

rad2deg = 180 / pi ;
deg2rad = pi / 180 ;
tol = 1e-12 ; % tolerance for dot product check
N_frames = data.Nimages ;
flipThresh = 65 ; % limit on raw angle values that, if exceeded, indicate need for large pert calculation
%smoothingParams = setSmoothingParams() ;
%--------------------------------------------------------------------------
% read out relevant data
AHat = data.AHat ;
rawPsi  = atan2(AHat(:,2),AHat(:,1));
rawBeta = asin(AHat(:,3));
rawPsiDeg = rad2deg*rawPsi ;
rawBetaDeg = rad2deg*rawBeta ;

% initialize data containers for rotated body axes
rotM_array = nan(3,3,N_frames) ;
AHat_rot = nan(N_frames, 3) ;
bodyPitch_rot = zeros(N_frames,1) ;
bodyYaw_rot = zeros(N_frames,1) ;

% initial rotation matrix to get AHat in first frame to align with x axis
rotM_init = eulerRotationMatrix(rawPsi(1), rawBeta(1),0) ;
AHat_rot_init = rotM_init*AHat(1,:)' ;
bodyPitch_rot(1) = rawBeta(1) ;
bodyYaw_rot(1) = rawPsi(1) ;

if abs(AHat_rot_init(1) - 1) > tol
    disp('error finding initial rotation matrix')
    keyboard
end

rotM_array(:,:,1) = rotM_init ;
AHat_rot(1,:) = AHat_rot_init' ; % we'll store this so we can check it works

% check to see if this is actually a large pert--if not we can avoid
% wind-up error by just using the raw values
if (max(abs(rawBetaDeg)) < flipThresh)
    bodyPitch = rawBetaDeg ;
    bodyYaw = rawPsiDeg ;
    for i = 2:N_frames
        dPitch = deg2rad*bodyPitch(i) ;
        dYaw = deg2rad*bodyYaw(i) ;
        rotM_curr = eulerRotationMatrix(dYaw, dPitch,0) ;
        rotM_array(:,:,i) = rotM_curr ;
    end
else
    % this means we need to do the frame-to-frame rotation estimate
    % loop through frames and try to rotate everybody onto x axis just
    % using yaw and pitch. then, hopefully, we can back out rest of the angles
    for i = 2:N_frames
        AHat_curr = AHat(i,:)' ;
        rotM_prev = squeeze(rotM_array(:,:,i-1)) ;
        AHat_rot1 = rotM_prev*AHat_curr ;
        
        dYaw = atan2(AHat_rot1(2),AHat_rot1(1));
        dPitch = asin(AHat_rot1(3)) ;
        
        rotM_curr = eulerRotationMatrix(dYaw, dPitch,0) ;
        AHat_rot(i,:) = (rotM_curr * AHat_rot1)' ;
        bodyYaw_rot(i) = dYaw ;
        bodyPitch_rot(i) = dPitch ;
        
        rotM_array(:,:,i) = rotM_curr*rotM_prev ;
    end
    
    % make sure to convert to degrees
    
    bodyPitch = cumsum(rad2deg*bodyPitch_rot) ;
    bodyYaw = cumsum(rad2deg*bodyYaw_rot) ;
end


%--------------------------------------------------------------------------
if debugFlag
    figure ;
    subplot(1,2,1)
    hold on
    plot(rawBetaDeg)
    plot(bodyPitch)
    %plot(bodyPitch_refined)
    xlabel('Frame')
    ylabel('Body Pitch Angle (deg)')
    legend({'Raw','From Quat.','refined'})
    title('Pitch')
    
    subplot(1,2,2)
    hold on
    plot(rawPsiDeg)
    plot(bodyYaw)
    %plot(bodyYaw_refined)
    xlabel('Frame')
    ylabel('Body Yaw Angle (deg)')
    legend({'Raw','From Quat.','refined'})
    title('Yaw')
    %plot(pitch_new_test)
    
    %     figure ;
    %     hold on
    %     plot(AHat_rot(:,1),'-')
    %     plot(AHat_rot(:,2),':')
    %     plot(AHat_rot(:,3),'--')
    %     title('Rotated AHat Check')
end

%--------------------------------------------------------------------------
%{
% check the sign of the pitch cumulative sum--for some reason it seems to
% switch oddly. just check which one matches best up until extreme rotation
extremePitchIdx = find((rawBeta > 70) | (rawBeta <= 0),1,'first') ;
if isempty(extremePitchIdx)
    extremePitchIdx = length(rawBeta) ;
end

extremeYawIdx = find((rawPsi > 170) | (rawPsi <10),1,'first') ;
if isempty(extremeYawIdx)
    extremeYawIdx = length(rawPsi) ;
end

%--------------------------------------------------------------------------
% check pitch
bodyPitch_pos = [rawBeta(1); pitch_cumsum_deg + rawBeta(1)] ;
bodyPitch_neg = [rawBeta(1); -1*pitch_cumsum_deg + rawBeta(1)] ;

pos_diff = norm(bodyPitch_pos(1:extremePitchIdx) - rawBeta(1:extremePitchIdx)) ;
neg_diff = norm(bodyPitch_neg(1:extremePitchIdx) - rawBeta(1:extremePitchIdx)) ;

if pos_diff < neg_diff
    bodyPitch = bodyPitch_pos ;
else
    bodyPitch = bodyPitch_neg ;
end

%--------------------------------------------------------------------------
% check yaw
bodyYaw_pos = [rawPsi(1); yaw_cumsum_deg + rawPsi(1)] ;
bodyYaw_neg = [rawPsi(1); -1*yaw_cumsum_deg + rawPsi(1)] ;

pos_diff = norm(bodyYaw_pos(1:extremeYawIdx) - rawPsi(1:extremeYawIdx)) ;
neg_diff = norm(bodyYaw_neg(1:extremeYawIdx) - rawPsi(1:extremeYawIdx)) ;

if pos_diff < neg_diff
    bodyYaw = bodyYaw_pos ;
else
    bodyYaw = bodyYaw_neg ;
end

%--------------------------------------------------------------------------
% refine estimates of yaw and pitch?
if refineFlag
    bodyPitch_refined = bodyPitch ;
    bodyYaw_refined = bodyYaw ;
    deltaPitch = zeros(N_frames, 1) ;
    deltaYaw = zeros(N_frames, 1) ;
    for i= 1:data.Nimages
        AHat_curr = data.AHat(i,:)' ;
        yaw_curr = deg2rad*bodyYaw(i) ;
        pitch_curr = deg2rad*bodyPitch(i) ;
        
        % unyaw matrix
        unyaw = eulerRotationMatrix(yaw_curr, 0, 0) ;
        
        % unpitch matrix
        unpitch = eulerRotationMatrix(0, pitch_curr, 0) ;
        
        % first undo the yaw angle rotation. if there is any residual yaw,
        % add it to the "refined" array
        AHat_unyawed = unyaw*AHat_curr ;
        delta_yaw  = rad2deg*atan2(AHat_unyawed(2),AHat_unyawed(1));
        bodyYaw_refined(i) = bodyYaw_refined(i) + delta_yaw ;

        % then rotated by the corrected angle so we can check pitch
        unyaw_new = eulerRotationMatrix(deg2rad*bodyYaw_refined(i), 0, 0) ;
        AHat_unyawed_new = unyaw_new * AHat_curr ;

        AHat_unpitched = unpitch * AHat_unyawed_new ;
        delta_pitch = rad2deg*asin(AHat_unpitched(3)) ;
%         Ahat_unrot = unpitch*unyaw*AHat_curr ;
%         ax_curr = cross([1; 0; 0],Ahat_unrot) ;
%
%         ang_curr = acos(dot([1; 0; 0],Ahat_unrot)) ;
%         axang_curr = [ax_curr', ang_curr] ;
%         quat_curr = axang2quat(axang_curr) ;
%
%         [delta_yaw, delta_pitch, ~] = quat2angle(quat_curr,'ZYX') ;
        bodyYaw_refined(i) = bodyYaw_refined(i) + delta_yaw ;
        bodyPitch_refined(i) = bodyPitch_refined(i) + delta_pitch ;
        deltaPitch(i) = delta_pitch ;
        deltaYaw(i) = delta_yaw ;
    end
else
    bodyPitch_refined = [] ;
    bodyYaw_refined = [] ;
    deltaPitch = [] ;
    deltaYaw = [] ;
end
%--------------------------------------------------------------------------
%}

end
