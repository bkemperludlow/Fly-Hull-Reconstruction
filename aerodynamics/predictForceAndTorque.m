% -------------------------------------------------------------------------
% function to plug in terms to the wingbeat-averaged equation of motion
% (which should match with calculated quasisteady forces/torques). Assumes
% we're in the BODY FIXED frame of reference
%
% NB: we have to define roll slightly differently here, since we
% need orthogonal axes for torque. So here, we'll refer to the forward
% direction (x) as the roll axis, as opposed to the long body axis
%
% For reference, see:
%   -Sun, 2014 review of insect flight
%       (http://dx.doi.org/10.1103/RevModPhys.86.615)
%   - Tsevi's wing pitch paper
%       (http://dx.doi.org/10.1103/PhysRevE.92.022712)
% -------------------------------------------------------------------------
function [forcePred, torquePred, wb] = predictForceAndTorque(bodyCM, ...
    bodyYPR, wbTimes, t, rotM_YP, rotM_roll, params, largePertFlag, ...
    smoothAnglesFlag, debugFlag)
% -------------------------------------------
%% inputs and params
if ~exist('rotM_YP','var') || isempty(rotM_YP)
    rotM_YP = [] ;
end
if ~exist('rotM_roll','var') || isempty(rotM_roll)
    rotM_roll = [] ;
end
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams() ;
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ;
end
if ~exist('smoothAnglesFlag','var') || isempty(smoothAnglesFlag)
    smoothAnglesFlag = false ;
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ;
end
% general params for force calc
g           = params.g ; % m/s^2 ; gravitational accel
body_mass   = params.body_mass ; % kg ; body mass
Iyy         = params.Iyy ;  % N*m*s^2 ; moment of inertia (pitch axis)
Ixx         = params.Ixx ;  % N*m*s^2 ; moment of inertia (roll axis)
Izz         = params.Izz ;  % N*m*s^2 ; moment of inertia (yaw axis)
Ixz         = 0 ; %params.Ixz ;  % N*m*s^2 ; moment of inertia (yaw/roll mixed)
beta_0      = params.beta_0 ;  % angle defining stroke plane horizontal
span        = params.span ;    % m ; wing span length

M2          = eulerRotationMatrix(0, -beta_0, 0) ;
I           = [Ixx, 0, Ixz; 0, Iyy, 0 ; Ixz, 0, Izz] ;

% indices for body angles
defineConstantsScript ;

% low pass filter freq for body angles
angle_filt_lvl = 100 ; % Hz % was 50 on 1/13/20

% % low pass filter for body angular velocity (if doing things numerically)
% d1 = designfilt('lowpassiir','FilterOrder',3,'SampleRate',8000, ...
% 'HalfPowerFrequency',400,'DesignMethod','butter');
% moving slope parameters
movslope_order = 2 ;
movslope_len = 100 ;

% timing info
N_frames = length(t) ;
dt = nanmean(diff(t)) ;

% make sure body angles are in radians
maxPitch = nanmax(bodyYPR(:,2)) ;
if (maxPitch > 10)
    bodyYPR = (pi/180).*bodyYPR ;
end
% -----------------------------------------------------------------
%% get body translational and rotational motion (velocity, accel)
% get body center of mass velocity and acceleration
[~, bodyVel, bodyAccel] = smoothBodyCM(bodyCM,'lowpass') ;

% smooth body euler angles -- note that we have to mess with yaw a bit so
% we don't smooth over a jump. this is the same as in "smoothBodyAngles.m"
if smoothAnglesFlag
    bodyYPR = smoothBodyEulerAng(bodyYPR, angle_filt_lvl, largePertFlag) ;
end

% get derivatives of euler angles w.r.t. time (vel, accel)
[angleVel, angleAccel] = diffBodyEulerAngles(bodyYPR, dt) ;

% -----------------------------------------------------------------
%% get rotation matrices for lab -> body-tied frame
% note that, for large perturbations, our wind-up calculations are
% incompatible with the standard form of euler rotation matrix, so we need
% to import the rotation matrices from large pert calculations

% first check if we have rot matrices as inputs or not
if isempty(rotM_YP) || isempty(rotM_roll)
    % if it's a large perturbation, we can't rotation matrices due to wind
    % up issues, so throw an exception
    if largePertFlag
        fprintf('Cannot calculate rotation matrices for large pert. \n')
        keyboard
    else
        % otherwise calculate rotation matrices from euler angles
        rotM_YP = zeros(N_frames, 3, 3) ;
        rotM_roll = zeros(N_frames, 3, 3) ;
        for k = 1:N_frames
            % current rotation matrix
            phiB   =  bodyYPR(k, 1) ;   % yaw, in radians
            thetaB =  bodyYPR(k, 2) ;   % pitch, in radians
            rhoB   =  bodyYPR(k, 3) ;   % roll, in radians
            
            % get yaw-pitch and roll rotation matrices (keeping them
            % separate for angular velocity calculations
            rotM_YP(k,:,:) = eulerRotationMatrix(phiB, thetaB, 0) ;
            rotM_roll(k,:,:) = eulerRotationMatrix(0,0,rhoB) ;
        end
    end
end

% -------------------------------------------------------------------------
%% apply lab -> body transformation to velocity and acceleration
% also do gravitational acceleration while we're at it
gravAccel = -g.*[zeros(N_frames, 2), ones(N_frames,1)] ;
% loop over frames
for m = 1:N_frames
    % read out rotation matrix
    %     M = M2*squeeze(rotM_roll(m,:,:))*squeeze(rotM_YP(m,:,:)) ;
    M = squeeze(rotM_roll(m,:,:))*squeeze(rotM_YP(m,:,:)) ;
    % rotate cm vel/accel
    bodyVel(m,:) = M*bodyVel(m,:)' ;
    bodyAccel(m,:) = M*bodyAccel(m,:)' ;
    
    % rotate gravitational accel
    gravAccel(m,:) = M*gravAccel(m,:)' ;
end

% % ------------------------------------------------------------------------
% % commented out for now: analytic form for gravitational acceleration in
% % body frame. doesn't work with wind-up matrices :(
% a_g_x = cos(beta_0).*sin(bodyYPR(:,2)) - ...
%     sin(beta_0).*cos(bodyYPR(:,2)).*cos(bodyYPR(:,3)) ;
% a_g_y = cos(bodyYPR(:,2)).*sin(bodyYPR(:,3)) ;
% a_g_z = sin(beta_0).*sin(bodyYPR(:,2)) + ...
%     cos(beta_0).*cos(bodyYPR(:,2)).*cos(bodyYPR(:,3)) ;
% gravAccel = -1.*g.*[a_g_x, a_g_y, a_g_z] ;

% -------------------------------------------------------------------
%% get body angular velocity and accerlation in body frame
% we should be able to calculate body angular velocity in body frame
% (omegaBodyBody) using just rotation matrices and time derivative of euler
% angles. if we're not doing a wind-up calculation (not large pert) we have
% an analytic form for this. otherwise use numerical approach
% use tsevi's/attila's/goldstein's method for calulating omega in the
% body-fixed frame (omega = rotation of body frame w.r.t. lab frame)

% regardless of method, need euler angle derivatives w.r.t. time
phiB_dot        = angleVel(:,1) ; phiB_ddot      = angleAccel(:,1) ;
thetaB_dot      = angleVel(:,2) ; thetaB_ddot    = angleAccel(:,2) ;
rhoB_dot        = angleVel(:,3) ; rhoB_ddot      = angleAccel(:,3) ;
modThetaB_dot   = -thetaB_dot   ; modThetaB_ddot = -thetaB_ddot ;

% now switch methods depending on pert size
if largePertFlag
    % get omegaBodyBody (body ang. vel. in body frame) by applying rotation
    % matrix as in goldstein
    omegaBodyBody = zeros(N_frames,3) ;
    for n = 1:N_frames
        % read out rotation matrices
        R_yp = squeeze(rotM_YP(n,:,:)) ; % yaw-pitch rotation
        R_r  = squeeze(rotM_roll(n,:,:)) ; % roll rotation
        
        % get phi, theta, and rho components of omega in body coordinates
        omega_phi = R_r*R_yp*[0; 0; phiB_dot(n)] ;
        omega_theta = R_r*[0 ; modThetaB_dot(n) ; 0] ;
        omega_rho = [rhoB_dot(n) ; 0 ; 0] ;
        
        % sum each component and add to array
        omega_sum = (omega_phi + omega_theta + omega_rho) ;
        omegaBodyBody(n,:) = omega_sum ;
    end
    
    % because we're doing this numerically, need to smooth (rotation
    % introduces noise)
    for dim = 1:3
        omegaBodyBody(:,dim) = smooth(omegaBodyBody(:,dim),7) ;
    end
    
    % need to also get angular acceleration. probably just need to do this
    % numerically :(
    omegaBodyBodyDot = zeros(N_frames,3) ;
    movslope_len = min([movslope_len, N_frames-1]) ;
    for dim = 1:3
        omegaBodyBodyDot(:,dim) = (1/dt)*movingslope(omegaBodyBody(:,dim),...
            movslope_len, movslope_order) ;
    end
    
    if (0)
        figure ;
        for dim = 1:3
            subplot(3,1,dim)
            hold on
            plot(t(2:end), (1/dt).*diff(omegaBodyBody(:,dim)))
            plot(t, omegaBodyBodyDot(:,dim))
            
            %plot(t, omegaBodyBodyDotOld(:,dim))
        end
    end
else
    % calculate sines and cosines of body euler angles for ease of writing
    % things out
    cphB=cos(bodyYPR(:,1))    ; sphB=sin(bodyYPR(:,1))   ; %#ok<NASGU>
    cthB=cos(-bodyYPR(:,2))   ; sthB=sin(-bodyYPR(:,2));
    crhB=cos(bodyYPR(:,3))    ; srhB=sin(bodyYPR(:,3))   ;
    
    % body angular velocity as measured in the body frame of reference.
    % note minus sign indefinition of sth above, so minus here is ok
    omegaBodyBody = [ rhoB_dot - phiB_dot.*sthB   , ...
        modThetaB_dot.*crhB  + phiB_dot.*cthB.*srhB  , ...
        -modThetaB_dot.*srhB + phiB_dot.*cthB.*crhB ] ;
    
    % derivative of body angular velocity w.r.t. time in body frame
    omegaBodyBodyDot = [rhoB_ddot - phiB_ddot.*sthB - phiB_dot.*modThetaB_dot.*cthB , ... % x component
        modThetaB_ddot.*crhB - modThetaB_dot.*rhoB_dot.*srhB + ...                        % y component
        phiB_ddot.*cthB.*srhB + phiB_dot.*modThetaB_dot.*sthB.*srhB + ...
        phiB_dot.*rhoB_dot.*cthB.*crhB, ...
        -1.*modThetaB_ddot.*srhB - modThetaB_dot.*rhoB_dot.*crhB + ...                   % z component
        phiB_ddot.*cthB.*crhB + phiB_dot.*modThetaB_dot.*sthB.*crhB - ...
        phiB_dot.*rhoB_dot.*cthB.*srhB] ;
end

% % need to rotate up into stroke plane horizontal frame (?)
% for q = 1:N_frames
%     omegaBodyBody(q,:) = M2 * omegaBodyBody(q,:)' ;
%     omegaBodyBodyDot(q,:) = M2 * omegaBodyBodyDot(q,:)' ;
% end
% ----------------------------------------------------------------------
%% average over wingbeats
% get wingbeat number info
N_wb = length(wbTimes) - 1 ;
wb0_ind = find(wbTimes < 0, 1, 'last') ;
if isempty(wb0_ind)
    wb0_ind = 1 ;
end
wb = (1:N_wb) - wb0_ind ;


% initialize arrays for wb-averaged quantities
bodyAccel_wba       = zeros(N_wb, 3) ;
bodyVel_wba         = zeros(N_wb, 3) ;
gravAccel_wba       = zeros(N_wb, 3) ;
omegaBody_wba       = zeros(N_wb, 3) ;
omegaBodyDot_wba    = zeros(N_wb, 3) ;

% loop over wingbeats (probably a better way to do averaging, but don't
% overthink it)
for j = 1:N_wb
    % get indices for cutOffTimes in data matrices
    [~, ind1] = min(abs(t - wbTimes(j))) ;
    [~, ind2] = min(abs(t - wbTimes(j+1))) ;
    
    ind = ind1:ind2 ;
    
    % average each variable
    bodyVel_wba(j,:) = nanmean(bodyVel(ind,:)) ;
    bodyAccel_wba(j,:) = nanmean(bodyAccel(ind,:)) ;
    
    gravAccel_wba(j,:) = nanmean(gravAccel(ind,:)) ;
    
    omegaBody_wba(j,:) = nanmean(omegaBodyBody(ind,:)) ;
    omegaBodyDot_wba(j,:) = nanmean(omegaBodyBodyDot(ind,:)) ;
    
end
% -----------------------------------------------------------------------
%% combine values via Newton-Euler equations
% perform calculations for linear acceleration
% forcePred = body_mass.*(bodyAccel_wba + ...
%     cross(omegaBody_wba, bodyVel_wba,2) - gravAccel_wba) ;
forcePred = body_mass.*(bodyAccel_wba - gravAccel_wba) ;
% perform calculations for rotational acceleration
torquePred = zeros(N_wb, 3) ;
for m = 1:N_wb
    torquePred(m,:) = cross(omegaBody_wba(m,:), (I*omegaBody_wba(m,:)')') + ...
        (I*omegaBodyDot_wba(m,:)')' ;
%         torquePred(m,:) =   (I*omegaBodyDot_wba(m,:)')' ;
end

% ---------------------------------------------------------
%% plot results?
if debugFlag
    labelList = {'x', 'y', 'z'} ;
    % plot forces
    figure ;
    hold on
    plot([wb(1), wb(end)], [0,0],'k--','HandleVisibility','off')
    for dim = 1:3
        %subplot(3,1,dim)
        plot(wb, forcePred(:,dim)./(body_mass*g))
    end
    legend(labelList)
    xlabel('Wingbeat')
    ylabel('F/mg')
    axis tight
    title('Forces')
    
    % plot forces
    figure ;
    hold on
    plot([wb(1), wb(end)], [0,0],'k--','HandleVisibility','off')
    for dim = 1:3
        %subplot(3,1,dim)
        plot(wb, torquePred(:,dim)./(body_mass*g*span))
    end
    legend(labelList)
    xlabel('Wingbeat')
    ylabel('T/mgl')
    axis tight
    title('Torques')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------
%% smooth body angles
% ("smoothBodyAngles.m" takes a data structure as input, so we'll just copy
% the methods for the data here)
function bodyYPR_smooth = smoothBodyEulerAng(bodyYPR, filt_lvl, ...
    largePertFlag)
% ------------------------------
% read out body angles
bodyYaw = bodyYPR(:,1) ;
bodyPitch = bodyYPR(:,2) ;
bodyRoll = bodyYPR(:,3) ;
% ------------------------------------------------------------------------
% The main thing to worry about is unwrapping body yaw (if we're not using
% wind-up). Otherwise try to identify jumps, but it's less clear cut
if ~largePertFlag
    % first subtract off early value of yaw, in case we start around a flip
    % point
    i1 = 1 ;
    i2 = min([10, length(bodyYaw)]) ;
    bodyYawInit = nanmedian(bodyYaw(i1:i2)) ;
    bodyYaw = bodyYaw - bodyYawInit ;
    
    % then restrict yaw to be in the range [-180, 180]
    for i = 1:length(bodyYaw)
        while (bodyYaw(i) < -pi)
            bodyYaw(i) = bodyYaw(i) + 2*pi ;
        end
        while (bodyYaw(i) > pi)
            bodyYaw(i) = bodyYaw(i) -  2*pi ;
        end
    end
    
    % add back initial value
    bodyYaw = bodyYaw + bodyYawInit ;
else
    bodyYaw = unwrap(bodyYaw) ;
end

%------------------------------------------------
% perform angle smoothing
pitchSmooth = filterEulerAngle(bodyPitch, filt_lvl ) ;
yawSmooth = filterEulerAngle(bodyYaw, filt_lvl ) ;
rollSmooth = filterEulerAngle(bodyRoll, filt_lvl ) ;

% -----------------------------------------------
% combine into output matrix
bodyYPR_smooth = [yawSmooth, pitchSmooth, rollSmooth] ;

end