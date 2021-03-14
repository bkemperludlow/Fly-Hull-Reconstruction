%--------------------------------------------------------------------------
% function to smooth euler angles describing the flie's body orientation.
% essentially does the same thing as the individual controller code bits,
% but combined here
%--------------------------------------------------------------------------
function [pitchSmooth, yawSmooth, rollSmooth] = ...
    smoothBodyAngles(data,largePertFlag, pitch_lvl, roll_lvl, yaw_lvl)
%--------------------------------------------------------------------------
%% params and inputs 
if ~exist('largePertFlag','var')
   largePertFlag = false ; % has wind-up error for small perts
end
% low pass filter levels for each body angle
if ~exist('pitch_lvl','var')
    pitch_lvl = [] ; 
end
if ~exist('roll_lvl','var')
    roll_lvl = [] ; 
end
if ~exist('yaw_lvl','var')
    yaw_lvl = [] ; 
end
defineConstantsScript
smoothingParams = setSmoothingParams() ; 
init_window = 20 ; 

RAD2DEG = (180/pi) ; 
DEG2RAD = (pi/180) ; 
% ------------------------------------------------------------
%% filter params
% determine whether we'll use default filter levels or inputs
if isempty(pitch_lvl)
    pitch_filt_lvl = smoothingParams.pitch_filt_lvl ;
else
    pitch_filt_lvl = pitch_lvl ; 
end
if isempty(roll_lvl)
    roll_filt_lvl = smoothingParams.roll_filt_lvl ;
else
    roll_filt_lvl = roll_lvl ; 
end
if isempty(yaw_lvl)
    yaw_filt_lvl = smoothingParams.yaw_filt_lvl ;
else
    yaw_filt_lvl = yaw_lvl ; 
end
%------------------------------------------------
%% load angle data
if largePertFlag 
    t = (1/data.params.fps) * ...
        (data.params.startTrackingTime : data.params.endTrackingTime) ;
   [bodyPitch, bodyYaw, rotM_YP] = calcPitchLargePert(data) ;
   [bodyRoll, ~, ~, ~, ~] = ...
        calcBodyRoll(data.rhoTimes, data.rollVectors, t, rotM_YP, data.params,...
         largePertFlag) ;
else
    bodyPitch = data.anglesLabFrame(:,BETA) ; 
    bodyYaw = data.anglesLabFrame(:,PHIB) ;
    bodyRoll = data.anglesLabFrame(:,RHO) ; 
end

% -----------------------------------------------
%% unwrap body yaw if we're not using wind-up
% (otherwise try to identify jumps, but it's less clear cut)
if ~largePertFlag
    % first subtract off early value of yaw, in case we start around a flip
    % point
    i1 = 1 ; 
    i2 = min([init_window, length(bodyYaw)]) ; 
    bodyYawInit = nanmedian(bodyYaw(i1:i2)) ; 
    bodyYaw = bodyYaw - bodyYawInit ;
    
    % then restrict yaw to be in the range [-180, 180]
    for i = 1:length(bodyYaw)
        while bodyYaw(i) < -180
            bodyYaw(i) = bodyYaw(i) + 360 ;
        end
        while bodyYaw(i) > 180
            bodyYaw(i) = bodyYaw(i) - 360 ;
        end
    end
    
    % add back initial value
    bodyYaw = bodyYaw + bodyYawInit ; 
else
    bodyYaw = RAD2DEG*unwrap(DEG2RAD*bodyYaw) ; 
end

%------------------------------------------------
%% perform angle smoothing
pitchSmooth = filterEulerAngle(bodyPitch, pitch_filt_lvl ) ;
yawSmooth = filterEulerAngle(bodyYaw, yaw_filt_lvl ) ;
rollSmooth = filterEulerAngle(bodyRoll, roll_filt_lvl ) ;

%------------------------------------------------
%% make sure arrays are column vectors
if size(pitchSmooth,1) < size(pitchSmooth,2)
    pitchSmooth = pitchSmooth' ; 
end
if size(yawSmooth,1) < size(yawSmooth,2)
    yawSmooth = yawSmooth' ; 
end
if size(rollSmooth,1) < size(rollSmooth,2)
    rollSmooth = rollSmooth' ; 
end

end