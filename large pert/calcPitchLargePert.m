%--------------------------------------------------------------------------
% In cases when the fly undergoes large pitch perturbations, use wind-up
% calculations to avoid weirdness at Euler angle limits
%
% NB: this is the cleaned up version--more commented code (that may have
% useful snippets) is in the "old" folder
%
% NB2: originally thought we should use full data struct, but allowing just
% AHat as input also
%--------------------------------------------------------------------------
function [bodyPitch, bodyYaw, rotM_cumarray, largePertFlag] = ...
    calcPitchLargePert(data, debugFlag)
%--------------------------------------------------------------------------
%% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ;
end

% deal with "data" input
if isstruct(data)
    AHat = data.AHat ;
elseif isnumeric(data)
    AHat = data ;
else
    fprintf('Cannot interpret input \n')
    keyboard
end
N_frames = size(AHat,1) ;
% general params
RAD2DEG = 180 / pi ;
DEG2RAD = pi / 180 ;
tol = 1e-12 ; % tolerance for dot product check


% params for excluding bad frames
filtFlag = true ;
diffThresh = 10 ; % in deg. max allowable frame-to-frame change in body angle
minChunkLen = 50 ; % in frames. minimum length for good chunk

% assume that we have a large perturbation, but change if checks fail
largePertFlag = true ;
%--------------------------------------------------------------------------
%% read out relevant data
rawBeta = asin(AHat(:,3));
rawPhi  = atan2(AHat(:,2),AHat(:,1));
rawPhiDeg = RAD2DEG*rawPhi ;
rawBetaDeg = RAD2DEG*rawBeta ;

% initialize data containers for rotated body axes
rotM_cumarray = nan(3,3,N_frames) ;
rotM_array = nan(3,3,N_frames) ;
AHat_rot = nan(N_frames, 3) ;
AHat_pre_rot = nan(N_frames, 3) ;
bodyPitch_rot = zeros(N_frames,1) ;
bodyYaw_rot = zeros(N_frames,1) ;

% eulEstMat = zeros(N_frames,6) ; % [phi1, phi2, theta1, theta2] ( only getting yaw and pitch right now)
% set up a boolean value that should tell us if we need to bail on the
% wind-up (because it's failing us)
failFlag = false ;
% -------------------------------------------------------------------------
%% initial rotation matrix to get AHat in first frame to align with x axis

rotM_init = eulerRotationMatrix(rawPhi(1), rawBeta(1),0) ;
bodyPitch_rot(1) = rawBeta(1) ;
bodyYaw_rot(1) = rawPhi(1) ;

AHat_rot_init = rotM_init*AHat(1,:)' ;

if abs(AHat_rot_init(1) - 1) > tol
    disp('error finding initial rotation matrix')
    keyboard
end

rotM_cumarray(:,:,1) = rotM_init ;
rotM_array(:,:,1) = rotM_init ;
AHat_rot(1,:) = AHat_rot_init' ; % we'll store this so we can check it works

% eulEstMat(1,:) = [ rawPhi(1),  rawPhi(1), rawBeta(1), rawBeta(1), nan, nan ] ;
% ------------------------------------------------------------------------
%% check to see if this is actually a large pert
%  ...if not we can avoid wind-up error by just using the raw values
largePertFlag  = largePertFlag & checkAHat(AHat) ;
if ~largePertFlag
    bodyPitch = rawBetaDeg ;
    bodyYaw = rawPhiDeg ;
    for i = 2:N_frames
        dPitch = DEG2RAD*bodyPitch(i) ;
        dYaw = DEG2RAD*bodyYaw(i) ;
        rotM_curr = eulerRotationMatrix(dYaw, dPitch,0) ;
        rotM_cumarray(:,:,i) = rotM_curr ;
    end
else
    % this means we need to do the frame-to-frame rotation estimate
    % loop through frames and try to rotate everybody onto x axis just
    % using yaw and pitch. then, hopefully, we can back out rest of the angles
    loop_int = 1 ;
    for i = (1+loop_int):loop_int:N_frames
        AHat_curr = AHat(i,:)' ;
        rotM_prev = squeeze(rotM_cumarray(:,:,i-loop_int)) ;
        AHat_rot1 = rotM_prev*AHat_curr ;
        
        if (sign(AHat_rot1(1)) < 0)
            failFlag = true ;
        end
        %dYaw = asin(AHat_rot1(2));
        dYaw = atan2(AHat_rot1(2),AHat_rot1(1));
        dPitch = asin(AHat_rot1(3)) ;
        
        rotM_curr = eulerRotationMatrix(dYaw, dPitch,0) ;
        
        bodyYaw_rot(i) = dYaw ;
        bodyPitch_rot(i) = dPitch ;
        
        AHat_rot(i,:) = (rotM_curr * AHat_rot1)' ;
        AHat_pre_rot(i,:) = AHat_rot1 ;
        rotM_array(:,:,i) = rotM_curr ;
        rotM_cumarray(:,:,i) = rotM_curr*rotM_prev ;
        
        %         [phi1, phi2, theta1, theta2, psi1, psi2] = myRot2Eul(rotM_curr*rotM_prev) ;
        %         eulEstMat(i,:) = [phi1, phi2, -1*theta1, -1*theta2, psi1, psi2] ;
    end
    
    % body pitch (needs to be in degrees)
    bodyPitchVel = bodyPitch_rot(2:end) ; %detrend(bodyPitch_rot(2:end)) ;
    bodyPitchVel = [0 ; bodyPitchVel] ;
    bodyPitch = RAD2DEG*(bodyPitch_rot(1) + cumtrapz(bodyPitchVel)) ;
    
    % body yaw (needs to be in degrees)
    bodyYawVel = bodyYaw_rot(2:end) ; %detrend(bodyPitch_rot(2:end)) ;
    bodyYawVel = [0 ; bodyYawVel] ;
    bodyYaw = RAD2DEG*(bodyYaw_rot(1) + cumtrapz(bodyYawVel)) ;
end

% ----------------------------------------------------------------
%% check for large wind-up failures
if failFlag
    % if we did catch some sort of failure, just use non-large-pert angles
    bodyPitch = rawBetaDeg ;
    bodyYaw = rawPhiDeg ;
    for i = 2:N_frames
        dPitch = DEG2RAD*bodyPitch(i) ;
        dYaw = DEG2RAD*bodyYaw(i) ;
        rotM_curr = eulerRotationMatrix(dYaw, dPitch,0) ;
        rotM_cumarray(:,:,i) = rotM_curr ;
    end
    largePertFlag = false ;
    
end
% -------------------------------------------------------------------------
%% filter out any large jumps?
% NB: this will not change rotation matrices, so wing angles should be
% unaffected
if filtFlag
    % get frame-to-frame change in body angles
    bodyYawDiff = [0; diff(bodyYaw)] ;
    bodyPitchDiff = [0; diff(bodyPitch)] ;
    
    % get indices where diff is larger than diffThresh
    bad_idx_yaw = (abs(bodyYawDiff) > diffThresh) ;
    bad_idx_pitch = (abs(bodyPitchDiff) > diffThresh) ;
    
    good_idx_yaw = idx_by_thresh(~bad_idx_yaw) ;
    good_idx_pitch = idx_by_thresh(~bad_idx_pitch) ;
    
    % get lengths of chunks where diff is find
    good_idx_yaw_len = cellfun(@(y) length(y), good_idx_yaw) ;
    good_idx_pitch_len = cellfun(@(y) length(y), good_idx_pitch) ;
    
    good_idx_yaw = good_idx_yaw(good_idx_yaw_len > minChunkLen) ;
    good_idx_pitch = good_idx_pitch(good_idx_pitch_len > minChunkLen) ;
    
    % initialize filtered arrays to then fill with good values
    yawFilt = nan(size(bodyYaw)) ;
    pitchFilt = nan(size(bodyPitch)) ;
    
    for ii = 1:length(good_idx_yaw)
        idx = good_idx_yaw{ii} ;
        yawFilt(idx) = bodyYaw(idx) ;
    end
    for jj = 1:length(good_idx_pitch)
        idx = good_idx_pitch{jj} ;
        pitchFilt(idx) = bodyPitch(idx) ;
    end
    
    % now interpolate over remaining nan values (frames identified as having too
    % large a jump)
    frames = 1:N_frames ;
    nan_idx_yaw = isnan(yawFilt) ;
    yawFilt_interp = interp1(frames(~nan_idx_yaw), ...
        yawFilt(~nan_idx_yaw), frames, 'pchip') ;
    
    nan_idx_pitch = isnan(pitchFilt) ;
    pitchFilt_interp = interp1(frames(~nan_idx_pitch), ...
        pitchFilt(~nan_idx_pitch), frames, 'pchip') ;
    
    bodyYaw = yawFilt_interp ;
    bodyPitch = pitchFilt_interp ;
end

% -------------------------------------------------------------------------
%% plot normal vs large pert calculation?
if debugFlag
    figure ;
    subplot(1,2,1)
    hold on
    plot(rawBetaDeg)
    plot(bodyPitch)
    %plot(bodyPitch_refined)
    xlabel('Frame')
    ylabel('Body Pitch Angle (deg)')
    legend({'Raw','From Wind Up'})
    title('Pitch')
    
    subplot(1,2,2)
    hold on
    plot(rawPhiDeg)
    plot(bodyYaw)
    %plot(bodyYaw_refined)
    xlabel('Frame')
    ylabel('Body Yaw Angle (deg)')
    legend({'Raw','From Wind Up'})
    title('Yaw')
    %plot(pitch_new_test)
    %{
    figure ;
    hold on
    plot(AHat_rot(:,1),'-')
    plot(AHat_rot(:,2),':')
    plot(AHat_rot(:,3),'--')
    title('Rotated AHat Check')
    %}
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------
%% function to check if AHat approaches singularity
function singularityFlag = checkAHat(AHat)
% get z component of AHat and see if it approaches +/- 1
z_tol = 0.05 ; % distance from z component of AHat approaching +/- 1
[minDist, ~] = min(abs(abs(AHat(:,3)) - 1)) ;

% also check for jumps in yaw
phi_tol = 0.25 ; % maximum allowed frame-to-frame change in yaw
rawPhi  = atan2(AHat(:,2),AHat(:,1));
maxPhiDiff = max(abs(diff(rawPhi))) ; 

% check for both conditions
singularityFlag = (minDist < z_tol) || (maxPhiDiff > phi_tol) ;
end

%% failed attempts:
%{
% -------------------------------------------------------------------------
%% experimenting with integrating rotational velocity about x, y, z axes
if (0)
%     % rotate AHat so first frame is along x axis
%     AHat_body = (rotM_init*AHat')' ;
    
    % take derivatives of AHat in x, y, and z coords
    movslope_order = 2 ;
    movslope_len = 100 ;
    
    v_x = movingslope(AHat(:,1),movslope_len, movslope_order) ;
    v_y = movingslope(AHat(:,2),movslope_len, movslope_order) ;
    v_z = movingslope(AHat(:,3),movslope_len, movslope_order) ;
    
    v_all = [v_x, v_y, v_z] ;
    
    omega = cross(AHat, v_all) ;
    omega_rot = zeros(size(omega)) ;
    % try to rotate omega into body frame?
    for k = 1:N_frames
        omega_rot(k,:) = squeeze(rotM_cumarray(:,:,k))*omega(k,:)' ;
        %omega_rot(k,:) = omega(k,:)' ;
    end
    % it's not always going to be one component of omega--need to think
    % about body frame of reference
    pitchEst = -1*RAD2DEG*cumtrapz(omega_rot(:,2)) + rawBetaDeg(1) ;
    yawEst = RAD2DEG*cumtrapz(omega_rot(:,3)) + rawPsiDeg(1) ;
    
    figure ;
    hold on
    plot(rawBetaDeg)
    plot(bodyPitch)
    plot(pitchEst)

    figure ;
    hold on
    plot(rawPsiDeg)
    plot(bodyYaw)
    plot(yawEst)
    
end

% -------------------------------------------------------------------------
%% check integrated angle vs standard
if (0)
    % initialize array of body axis vectors aligned to x axis
    AHat_test = [ones(N_frames,1), zeros(N_frames,2)] ;
    AHat_test_rot_lp = zeros(size(AHat_test)) ;
    AHat_test_rot_reg = zeros(size(AHat_test)) ;
    for j = 1:N_frames
        R_lp = eulerRotationMatrix(DEG2RAD*bodyYaw(j), ...
            DEG2RAD*bodyPitch(j),0) ;
        R_reg = eulerRotationMatrix(rawPsi(j), ...
            rawBeta(j),0) ;

%         R_lp = eulerRotationMatrix(DEG2RAD*bodyYaw(j),0,0) ;
%         R_reg = eulerRotationMatrix(rawPsi(j), 0,0) ;
        %AHat_test_rot(j,:) = squeeze(rotM_cumarray(:,:,j))'*AHat_test(j,:)' ;
        AHat_test_rot_lp(j,:) = R_lp'*AHat_test(j,:)' ;
        AHat_test_rot_reg(j,:) = R_reg'*AHat_test(j,:)' ;
    end
    
    test_dot_lp = dot(AHat_test_rot_lp,AHat,2) ;
    test_dot_reg = dot(AHat_test_rot_reg,AHat,2) ;
    test_dot = dot(AHat_test_rot_reg,AHat_test_rot_lp,2) ;
    figure
    hold on
    plot(test_dot_reg)
    plot(test_dot_lp)
    legend({'Regular','LargePert'})
    
    figure
    plot(test_dot)
    
    unyaw_init = eulerRotationMatrix(rawPsi(1),0,0) ;
    AHat_body = (unyaw_init*AHat')' ;
    AHat_XZ = [AHat_body(:,1), zeros(N_frames,1), AHat_body(:,3)]  ;
    AHat_XZ = AHat_XZ ./ myNorm(AHat_XZ) ;
    checkPitch = unwrap(atan2(AHat_XZ(:,3), AHat_XZ(:,1))) ;
    figure ;
    hold on
    plot(checkPitch) ;
    plot(rawBeta) ;
%     toXZPlane = zeros(3,3,N_frames) ;
%     toXZPlane(:,:,1) = rotM_init ;
    for k = 1:N_frames
        angle_to_x = acos(AHat_body(k,1)) ;
    end
end
%}