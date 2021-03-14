% -------------------------------------------------------------------------
% function to calculate the velocity of a wing tip in body frame 
% coordinates. used in the calucation of aerodynamic forces. contributing 
% to this are:
%   - the wing's rotation relative to the body 
%   - the body's rotation
%   - the body's center of mass velocity
%   - the velocity of the wing hinge location relative to the body cm
%
% -----------
% INPUTS:
%   -wingAngleMat: Nx3 matrix containing wing Euler angles in the body frame,
%   in the order [stroke, deviation, (wing) pitch]. *SHOULD BE IN RADIANS
%   ***NB: assuming here that sign is already flipped on stroke angle!
%
%   -t: Nx1 vector with time of each frame (in seconds)
%
%   -wingSide: char indicating whether we're doing the right or left wing
%   (either 'R' or 'L')
%
%   -params: stucture containing wing morphological parameters. This
%   structure can be set to default values via XX.m, but can be adjusted as
%   necessary to accomodate different wing types
%
%   -bodyCM: Nx3 matrix contaning the body position in lab frame. Should be
%   in meters / second
%
%   -bodyYPR: Nx3 matrix containing body Euler angles in the order [yaw,
%   pitch, roll].
%
% -------------
% OUTPUTS:
%   - U_t: wing tip velocity expressed in body frame coordinates
% -------------------------------------------------------------------------
function [U_t, omegaWingBody, bodyVel_body, hingeVel_body, omegaWingWing] = ...
    calcWingTipVelocity(wingAngleMat, t, params, wingSide, bodyCM, bodyYPR)
% -------------------
%% inputs and params
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams() ;
end
if ~exist('bodyCM','var') || isempty(bodyCM)
    bodyCM = [] ;
end
if ~exist('bodyYPR','var') || isempty(bodyYPR)
    bodyYPR = [] ;
end

% ---------------------------------------------------------------
% read in params from struct
span = params.span ; %in meters (was .002)
hinge_vec = params.hinge_vec ; 
thorax_width = params.thorax_width ; 
beta_0 = params.beta_0 ; % resting pitch angle in radians
% fitType = params.fitType ; % 'cubicspline' | 'smoothingspline' (for alpha)
smoothWin = params.bodyVelSmoothWin ; % window size for smoothing body frame velocity
diffType = params.wingKinDiffType ; % method to differentiate wing kinematics

dt = nanmean(diff(t)) ; % timing difference
N_frames = size(wingAngleMat,1) ; % also get number of points

% ---------------------------------------------------------------------
%% need to swap the sign on some stuff depending on wing side
if contains(wingSide,'R','IgnoreCase',true)
    wingSign = -1 ;
elseif contains(wingSide,'L','IgnoreCase',true)
    wingSign = +1 ;
else
    fprintf('Invalid wingSide entry: %s \n', wingSide)
    keyboard
end

% matrix to pitch down by thetaB0 (~45deg) w.r.t body axis
% (transforms from body to stroke plane coords)
body2sp = eulerRotationMatrix(0, -1*beta_0, 0) ;

% -----------------------------------------------------------------------
%% read out wing Euler angles and get rotational velocities
% read out euler angles from matrix, fill in nan values, get derivatives
wingAngleMat = interpAngleNans(wingAngleMat) ; 

phi = wingAngleMat(:,1) ;
theta = wingAngleMat(:,2) ;
psi = wingAngleMat(:,3) ;

[angleVel, ~] = diffWingEulerAngles(wingAngleMat, dt, diffType) ; 

phi_dot = angleVel(:,1) ;
theta_dot = angleVel(:,2) ;
psi_dot = angleVel(:,3) ;

% need to account for sign differences in left vs right pitch angle
if contains(wingSide,'L','IgnoreCase',true)
   psi_dot = -1.*psi_dot ; 
   psi = pi - psi ; 
end

%psi_dot = wingSign*psi_dot ; % account for right-handed vs left-handed
[omegaWingBody, omegaWingWing] = calcWingRotVel(phi, theta, psi, ...
    phi_dot, theta_dot, psi_dot) ;

% ---------------------------------------------------
%% get span and chord hat vectors (body coordinates)
[spanHat, ~] = getWingVecsFromAng(wingAngleMat) ;

%--------------------------------------------------------------------------
%% get body cm velocity and body rotational velocity in body frame
% (if data is provided)
if ~isempty(bodyCM) && ~isempty(bodyYPR)
    % read out body Euler angles
    phiB = bodyYPR(:,1) ; % yaw NB: all should be in radians
    thetaB = bodyYPR(:,2) ; % pitch
    rhoB = bodyYPR(:,3) ;  % roll
    
    % -------------------------
    % BODY CM VELOCITY
    % -------------------------
    % get velocity of body CM in lab frame
    [~, bodyVel_lab, ~] = smoothBodyCM(bodyCM) ;
    
    % to transform lab frame CM velocity to stroke plane coords, rotate
    bodyVel_body = zeros(size(bodyVel_lab)) ;
    lab2body_all = nan(3,3,N_frames) ;
    for i = 1:N_frames
        M1 = eulerRotationMatrix(phiB(i),thetaB(i),rhoB(i) ) ; % to strict body axis
        lab2body = body2sp * M1 ; % full rotation (M2 rotates to thetaB0rad)
        
        % rotate body velocity vector
        bodyVel_body(i,:) = lab2body*bodyVel_lab(i,:)' ;
        
        % store rotation matrix
        lab2body_all(:,:,i) = lab2body ;
    end
    
    % rotating by body angles introduces extra noise--smooth body frame
    % velocity
    for dim = 1:3
       bodyVel_body(:,dim) = smooth(bodyVel_body(:,dim),smoothWin,'rloess') ;  
    end
    % -------------------------
    % BODY ANGULAR VELOCITY
    % -------------------------
    % (for body frame wingtip velocity, need to add in contribution of body
    % frame rotation in the form cross(omega, r_tip) )
    % get derivatives of euler angles
    angleMat = [phiB, thetaB, rhoB] ; 
    [angleVel, ~] = diffBodyEulerAngles(angleMat, dt) ; 
    
    phiB_dot   = angleVel(:,1) ;
    thetaB_dot = angleVel(:,2) ;
    rhoB_dot   = angleVel(:,3) ;
    
    % use tsevi's/attila's/goldstein's method for calulating omega in the
    % body-fixed frame (omega = rotation of body frame w.r.t. lab frame)
    cphB=cos(phiB)    ; sphB=sin(phiB)   ; %#ok<NASGU>
    cthB=cos(-thetaB) ; sthB=sin(-thetaB);
    crhB=cos(rhoB)    ; srhB=sin(rhoB)   ;
    
    modThetaB_dot    = -thetaB_dot ;
    
    % body angular velocity as measured in the body frame of reference.
    omegaBodyBody = [ rhoB_dot - phiB_dot.*sthB   , ... % note minus sign is definition of sth above, so minus here is ok
        modThetaB_dot.*crhB  + phiB_dot.*cthB.*srhB  , ...
        -modThetaB_dot.*srhB + phiB_dot.*cthB.*crhB ] ;
    
    % -----------------------------------------------------
    % add body angular velocity to wing angular velocity
    omegaWingBody = omegaWingBody + omegaBodyBody ;
    
    % -------------------------
    % HINGE VELOCITY
    % -------------------------
    % hinge location changes depending on whether or not we're using right
    % vs left wing. here expressed in body coords
    hinge_vec = hinge_vec + wingSign.*[0 ; thorax_width ; 0] ; 
    hinge_vec_sp = body2sp*hinge_vec ; 
    % hinge velocity due to body rotation + it being offset from CM
    hingeVel_body = cross(omegaBodyBody, repmat(hinge_vec_sp',N_frames,1));
    
else
    bodyVel_body = 0 ;
    hingeVel_body = 0 ;
end

% ---------------------------------------------------------------
%% get wing tip lab-frame velocity (but expressed in body coords)
% get wing tip position (in body-tied frame)
r_tip = span.*spanHat ;

% get tip velocity via omega cross r and add in both body and hinge
% velocity (if we have them)
U_t = cross(omegaWingBody, r_tip) + bodyVel_body + hingeVel_body ;


end