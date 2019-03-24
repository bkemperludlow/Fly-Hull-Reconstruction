%--------------------------------------------------------------------------
% In cases when the fly undergoes large pitch perturbations, use a
% quaternion conversion to avoid gimbal lock
%--------------------------------------------------------------------------

function [bodyPitch, bodyYaw, deltaPitch, deltaYaw] = ...
    calcPitchLargePert(data, debugFlag, refineFlag)
%--------------------------------------------------------------------------
% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ;
end
if ~exist('refineFlag','var')
    % if true, try to account for (what i assume are) wind up errors in the
    %  angle calculation
    refineFlag = true ;
end

rad2deg = 180 / pi ;
deg2rad = pi / 180 ;
N_frames = data.Nimages ; 
%smoothingParams = setSmoothingParams() ;
%--------------------------------------------------------------------------
% read out relevant data
AHat = data.AHat ;
rawPsi  = rad2deg*atan2(AHat(:,2),AHat(:,1));
rawBeta = rad2deg*asin(AHat(:,3));

% estimate pitch angle via conversion to quaternions
angleMat = acos(dot(AHat(2:end,:),AHat(1:end-1,:),2)) ;
axisMat = cross(AHat(2:end,:),AHat(1:end-1,:),2) ;
axAngMat = [axisMat, angleMat] ;
quatMat = axang2quat(axAngMat) ;

quatMatCum= zeros(size(quatMat)) ;
quatMatCum(1,:) = quatMat(1,:) ;
for q = 2:size(quatMatCum,1)
    quatMatCum(q,:) = quatmultiply(quatMatCum(q-1,:),quatMat(q,:)) ;
end

% get cumulative sum of rotations, rather than individual frame readout
[yaw_diff_temp, pitch_diff_temp, ~] = quat2angle(quatMat,'ZYX') ; %'ZYX'
pitch_cumsum = cumsum(pitch_diff_temp) ;
pitch_cumsum_deg = rad2deg*pitch_cumsum ;

yaw_cumsum = cumsum(yaw_diff_temp) ;
yaw_cumsum_deg = rad2deg*yaw_cumsum ;

%--------------------------------------------------------------------------
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
if debugFlag
    figure ;
    subplot(1,2,1)
    hold on
    plot(rawBeta)
    plot(bodyPitch)
    plot(bodyPitch_refined)
    xlabel('Frame')
    ylabel('Body Pitch Angle (deg)')
    legend({'Raw','From Quat.','refined'})
    title('Pitch')
    
    subplot(1,2,2)
    hold on
    plot(rawPsi)
    plot(bodyYaw)
    plot(bodyYaw_refined)
    xlabel('Frame')
    ylabel('Body Yaw Angle (deg)')
    legend({'Raw','From Quat.','refined'})
    title('Yaw')
    %plot(pitch_new_test)
end


end
