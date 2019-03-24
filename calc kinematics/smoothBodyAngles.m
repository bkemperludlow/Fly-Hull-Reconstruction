%--------------------------------------------------------------------------
% function to smooth euler angles describing the flie's body orientation.
% essentially does the same thing as the individual controller code bits,
% but combined here
%--------------------------------------------------------------------------
function [pitchSmooth, yawSmooth, rollSmooth] = ...
    smoothBodyAngles(data,largePertFlag)
%--------------------------------------------------------------------------
%% params and inputs 
if ~exist('largePertFlag','var')
   largePertFlag = true ; % gives equivalent results for small perts
end

defineConstantsScript
smoothingParams = setSmoothingParams() ; 

%------------------------------------------------
%% load angle data
if largePertFlag 
    t = (1/data.params.fps) * ...
        (data.params.startTrackingTime : data.params.endTrackingTime) ;
   [bodyPitch, bodyYaw, rotM_YP] = calcPitchLargePert(data) ;
   [bodyRoll, ~, ~, ~] = ...
    calcBodyRoll(data.rhoTimes, data.rollVectors, t, rotM_YP, data.params) ;
else
    bodyPitch = data.anglesLabFrame(:,BETA) ; 
    bodyYaw = data.anglesLabFrame(:,PHIB) ;
    bodyRoll = data.anglesLabFrame(:,RHO) ; 
end

%------------------------------------------------
%% perform angle smoothing
pitchSmooth = filterEulerAngle(bodyPitch, smoothingParams.pitch_filt_lvl) ;
yawSmooth = filterEulerAngle(bodyYaw, smoothingParams.yaw_filt_lvl) ;
rollSmooth = filterEulerAngle(bodyRoll, smoothingParams.roll_filt_lvl) ;

end