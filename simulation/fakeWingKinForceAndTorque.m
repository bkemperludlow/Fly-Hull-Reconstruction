% -------------------------------------------------------------------------
% function to calculate body coordinate forces and torques from
% parameterized (fake) wing kinematics.
%
% ** try for full 3D or just longitudinal flight?
%
% INPUTS:
%   - wing_kin_params: 12x2 matrix with parameter values for fake wing
%   kinenmatics of the right and left wing. Columns 1 and 2 correspond to
%   right and left wing (respectively). Each column is of the form:
%     wing_kin_params(:,i) = [omega; phi_f; phi_b; K; theta_0; theta_m;...
%                             del_theta; psi_0; psi_m; del_psi; C]
%
%   - body_params: vector containing values for fly morphological
%   measurements in the form:
%        body_params = [span, hinge_x, thorax_width]
%   - t: Nx1 vector of time values over which to evaluate the force/torque
%   - velBodyBody: Tx3 matrix of body center of mass velocity, in body
%     coordinates
%   - bodyYPR: Tx3 matrix of body Euler angles in the form:
%       [phiB, thetaB, rhoB]
%   - bodyYPR_dot:  Tx3 matrix of body Euler angle derivatives w.r.t. time
%       in the form: [phiB_dot, thetaB_dot, rhoB_dot]
%
% OUTPUTS:
%   - F_tot: total quasi-steady force vector, expressed in body coords
%   - T_tot: total quasi-steady torque vector, expressed in body coords
% -------------------------------------------------------------------------
function [F_tot, T_tot] = fakeWingKinForceAndTorque(wing_kin_params, ...
    body_params, t, velBodyBody, bodyYPR, bodyYPR_dot, params, ...
    addedMassFlag, plotFlag)
% ----------------------------
%% inputs and params
if ~exist('velBodyBody','var') || isempty(velBodyBody)
    velBodyBody = zeros(length(t), 3) ;
end
if ~exist('bodyYPR','var') || isempty(bodyYPR)
    bodyYPR = zeros(length(t), 3) ;
end
if ~exist('bodyYPR_dot','var') || isempty(bodyYPR_dot)
    bodyYPR_dot = zeros(length(t), 3) ;
end
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams('span',body_params(1), ...
        'r_hinge', body_params(2)) ; 
end
if ~exist('addedMassFlag','var') || isempty(addedMassFlag)
    addedMassFlag = false ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = false ;
end

% plot wing kin?
debugFlag = false ;

% smooth wing tip velocity? probably unnecessary for fake data
smoothFlag = false ; 

% make list of right and left side to help looping over them
wingSideList = {'R', 'L'} ;
wingSignList = [-1, 1] ;   % for adding y component to hinge vec

% initialize outputs
N_pts = length(t) ;
F_tot = zeros(N_pts, 3) ;
T_tot = zeros(N_pts, 3) ;

% structure of indices for params vecs
ind_struct = defineParamIndices() ;

% read out body params
span = body_params(ind_struct.span) ;
hinge_x = body_params(ind_struct.hinge_x) ;
thorax_width = body_params(ind_struct.thorax_width) ;

% quasi-steady parameter struct
hinge_x_vec = params.hinge_vec ; % read out updated hinge vec
thetaB0 = params.beta_0 ; % stroke plane to body coords angle. also "rest" angle
body2sp = eulerRotationMatrix(0, -thetaB0, 0) ; 
sp2body = eulerRotationMatrix(0, thetaB0, 0) ; % stroke plane to body rotation matrix

% if no body Euler angles provided, set pitch to equilibrium value
if all(bodyYPR == 0, 'all') && all(bodyYPR_dot == 0, 'all')
   bodyYPR(:,2) = thetaB0.*ones(length(t),1) ; 
end

% if wing_kin_params has only a single column, copy it so that we have same
% kinematics for left and right wing
if (size(wing_kin_params,2) < 2)
    wing_kin_params = repmat(wing_kin_params,1,2) ; 
end
% ----------------------------------------------------------------------
%% get body rotational velocity from Euler angles and their derivatives
% read out body Euler angles and their derivatives w.r.t. time
phiB = bodyYPR(:,1) ; 
thetaB = bodyYPR(:,2) ; 
rhoB = bodyYPR(:,3) ; 

phiB_dot = bodyYPR_dot(:,1) ; 
thetaB_dot = bodyYPR_dot(:,2) ; 
rhoB_dot = bodyYPR_dot(:,3) ; 

% write out sines and cosines for convenience and to do sign swap on theta
cphB=cos(phiB)    ; sphB=sin(phiB)   ; %#ok<NASGU>
cthB=cos(-thetaB) ; sthB=sin(-thetaB);
crhB=cos(rhoB)    ; srhB=sin(rhoB)   ;

modThetaB_dot    = -thetaB_dot ;

% body angular velocity as measured in the body frame of reference.
omegaBodyBody = [ rhoB_dot - phiB_dot.*sthB   , ... % note minus sign is definition of sth above, so minus here is ok
    modThetaB_dot.*crhB  + phiB_dot.*cthB.*srhB  , ...
    -modThetaB_dot.*srhB + phiB_dot.*cthB.*crhB ] ;

% ----------------------------------------------------------
%% loop over right and left wing side to get force/torque
for k = 1:length(wingSideList)
    % current wing side (right vs left)
    wingSide = wingSideList{k} ;
    wingSign = wingSignList(k) ;
    
    % ----------------------------
    % get current wing kin
    [wingAngleMat, wingVelMat] = getFakeWingKin(wing_kin_params(:,k), ...
        wingSide, t, debugFlag) ;
    
    % get span and chord unit vectors based on these angles
    [spanHat, chordHat] = getWingVecsFromAng(wingAngleMat) ;
    
%     % convert span and chord to body coordinates (will be output in stroke
%     % plane coordinates)
%     spanHatBody = (M_thetaB0'*spanHat')' ; 
%     chordHatBody = (M_thetaB0'*chordHat')' ; 
    % ------------------------------------------------------------------
    %% calculate wing tip velocity (based on wing AND body kinematics)
    % wing rotational velocity
    [omegaWingBody, omegaWingWing] = calcWingRotVel(wingAngleMat(:,1), ...
        wingAngleMat(:,2), wingAngleMat(:,3), wingVelMat(:,1), ...
        wingVelMat(:,2), wingVelMat(:,3)) ;
    
    % add body rotational velocity to wing rotational velocity
    omegaWingBody = omegaWingBody + omegaBodyBody ;
    
    % also get rotational velocity of hinge
    hinge_vec = (hinge_x_vec + wingSign.*[0 ; thorax_width ; 0])' ;
    hinge_vec = repmat(hinge_vec, size(omegaBodyBody,1), 1) ;
    hinge_vec_sp = (body2sp*hinge_vec')' ; 
    velHingeBody = cross(omegaBodyBody, hinge_vec_sp) ;
    
    % rotate body frame cm velocity to stroke plane 
    velBodySp = (body2sp* velBodyBody')' ;
    
    % wing tip positions
    % wing_tip = span.*spanHat ;
    wing_tip = span.*spanHat ;
    
    % sum terms to get total wing tip velocity IN STROKE PLANE COORDINATES
%     U_t = cross(omegaWingBody, wing_tip) + velHingeBody + velBodyBody ;
    U_t = cross(omegaWingBody, wing_tip) + velHingeBody + velBodySp ;
    
%     omegaWingSp = (body2sp*omegaWingBody')' ; 
%     velBodySp = (body2sp* velBodyBody')' ;
%     velHingeSp = (body2sp* velHingeBody')' ;
%     U_t = cross(omegaWingSp, wing_tip) + velHingeSp + velBodySp ;
    % ------------------------------------------------------------------
    %% get angle of attack and lift/drag coefficients
    % unit vectors giving wing tip velocity direction
    Uhat = U_t./myNorm(U_t) ;
    
    % angle of attack
    alpha = acos(dot(Uhat, chordHat,2)) ;
%     alpha = acos(dot(Uhat, chordHatBody,2)) ;
    % ------------------------------------------------------------------
    %% calculate forces
    [F_sum, ~, ~, ~, ~, ~, ~] = calcQuasiSteadyTerms(U_t, spanHat,...
        chordHat, omegaWingWing, omegaWingBody, alpha, params, ...
        addedMassFlag, smoothFlag) ;
%     [F_sum, ~, ~, ~, ~, ~, ~] = calcQuasiSteadyTerms(U_t, spanHatBody,...
%         chordHatBody, omegaWingWing, omegaWingBody, alpha, params, ...
%         addedMassFlag, smoothFlag) ;
    % ------------------------------------------------------------------
    %% get torques from forces and wing/body vecs
    % center of pressure in stroke plane frame
    CoP = 0.7*span.*spanHat ; %centers of pressure are 70% along span
    
    % transform to body frame coordinates
    CoP = (sp2body*CoP')' ;
    
    % add (center of mass -> hinge) vector
    CoP = hinge_vec + CoP ;
    
    % calculate torque
    T_sum = cross(CoP, F_sum) ;
    
    % ------------------------------------------
    %% add current wing side to total
    F_tot = F_tot + F_sum ;
    T_tot = T_tot + T_sum ;
end

% ---------------------------------
%% plot results?
if plotFlag
   fprintf('Plotting under construction! \n')
   return 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------
%% define indices for param vecs
function ind_struct = defineParamIndices()
% write out indices in wing params vec
omega_ind = 1 ;

phi_f_ind = 2 ;
phi_b_ind = 3 ;
K_ind = 4 ;

theta_0_ind = 5 ;
theta_m_ind = 6 ;
del_theta_ind = 7 ;

psi_0_ind = 8 ;
psi_m_ind = 9 ;
del_psi_ind = 10 ;
C_ind = 11 ;

% write out indices in body params vec
span_ind = 1 ;
hinge_x_ind = 2 ;
thorax_width_ind = 3 ;

% add indices to struct
ind_struct = struct() ;
ind_struct.omega = omega_ind ;

ind_struct.phi_f = phi_f_ind ;
ind_struct.phi_b = phi_b_ind ;
ind_struct.K = K_ind ;

ind_struct.theta_0 = theta_0_ind ;
ind_struct.theta_m = theta_m_ind ;
ind_struct.del_theta = del_theta_ind ;

ind_struct.psi_0 = psi_0_ind ;
ind_struct.psi_m = psi_m_ind ;
ind_struct.del_psi = del_psi_ind ;
ind_struct.C = C_ind ;

ind_struct.span = span_ind ;
ind_struct.hinge_x = hinge_x_ind ;
ind_struct.thorax_width = thorax_width_ind ;


end