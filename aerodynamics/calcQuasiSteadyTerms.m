% -------------------------------------------------------------------------
% calculate quasi-steady force terms (translational, rotational, and added
% mass)
%
% INPUTS:
%
%
% OUTPUTS:
%   - F_t: translational force
%   - F_rot: rotational force
%   - F_a: added-mass force
% -------------------------------------------------------------------------
function [F_tot, F_t, F_rot, F_a, Lhat, Dhat, wingNormalHat] = ...
    calcQuasiSteadyTerms(U_t, spanHat, chordHat, omegaWingWing, ...
        omegaWingBody, alpha, params, addedMassFlag, smoothFlag)
% ----------------------------
%% inputs and params
if ~exist('alpha','var') || isempty(alpha)
    % if angle of attack (alpha) not provided as input, calculate it
    Uhat = U_t./myNorm(U_t) ;
    alpha = acos(Uhat, chordHat) ;
end
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams ;
end
if ~exist('addedMassFlag','var') || isempty(addedMassFlag)
    addedMassFlag = false ; % calculate added mass force?
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = false ; % smooth wing tip velocity?
end

% read out params from struct
airDensity = params.airDensity ; %kg/m^3
span = params.span ; %in meters (was .002)
mean_chord = params.mean_chord ; % in meters
S = params.wing_area ; % wing area in m^2
r22_S = params.r22_S ; % second moment of wing area
beta_0 = params.beta_0 ; % resting pitch angle in radians

C_rot_theo = params.C_rot_theo ; % coefficient of rotational force (Sane & Dickinson, 2002)
c2_r_int = params.c2_r_int ; % area moment of some kind 
c2_int = params.c2_int ; % area moment of some kind

% matrix to pitch down by thetaB0 (~45deg) w.r.t body axis
%M2 = eulerRotationMatrix(0, -1*beta_0, 0) ;
sp2body = eulerRotationMatrix(0, beta_0, 0) ; % stroke plane to body rotation matrix
% ------------------------------------------------
%% translational force
% lift and drag coeffs
[C_L, C_D] = getForceCoeffs(alpha, params) ;

% squared magnitude of wing tip velocity
Usq = sum(U_t.*U_t,2) ; 

% smooth Usq?
if smoothFlag
   Usq = smooth(Usq, 3) ;  
end

% ---------------------
% force magnitudes
fLift = 0.5 * airDensity * S * r22_S * Usq .* C_L ;
fDrag = 0.5 * airDensity * S * r22_S * Usq .* C_D ;

% direction vectors
Uhat = U_t./myNorm(U_t) ; % unit vector of wing tip velocity
Lhat = cross(Uhat, spanHat) ; % perpendicular to Uhat
Lhat = Lhat ./ myNorm(Lhat) ; % but still need to determine sign
% in body fixed frame, should always have positive z component
Lsign = sign(Lhat(:,3)) ;
Lhat((Lsign < 0),:) = -1*Lhat((Lsign < 0),:) ;

Dhat = -Uhat ;

% combine to get full force vectors
F_L = fLift.*Lhat ;
F_D = fDrag.*Dhat ;

F_t = F_L + F_D ;

% ------------------------------------------------
%% rotational force
% wing rotation about span from (Sane & Dickinson, 2002)
omega_rot = omegaWingWing(:,1) ;

% magnitude of rotational force
fRot = C_rot_theo * airDensity * sqrt(Usq) .* abs(omega_rot) .* ...
    mean_chord^2 .* span .* c2_r_int ;

% should act perpendicular to wing surface
wingNormalHat = cross(spanHat, chordHat) ; % need to check sign
wingNormalHat = wingNormalHat./myNorm(wingNormalHat) ;

% one way of ensuring we have the right sign for wing normal is to compare
% it to the total translational force vector--this should also be normal to
% wing, in general
F_t_hat = F_t./myNorm(F_t) ;
wingNormalHat = sign(dot(wingNormalHat, F_t_hat,2)).*wingNormalHat ;

% get full rotational force vector
F_rot = fRot .* wingNormalHat ;

% ------------------------------------------------
%% added mass force
if addedMassFlag
    % added mass force has 2 terms. From (Sane & Dickinson, 2002). Note 
    % that the "phi" in their notation refers to "angular position" of 
    % wing, not just stroke angle (I think). So first thing is to get 
    % (d/dt)omegaBodyWing :
    c_omega_wb_x = fit(t, omegaWingBody(:,1), fitType) ;
    c_omega_wb_y = fit(t, omegaWingBody(:,2), fitType) ;
    c_omega_wb_z = fit(t, omegaWingBody(:,3), fitType) ;
    
    omega_wb_x_dot = differentiate(c_omega_wb_x,t) ;
    omega_wb_y_dot = differentiate(c_omega_wb_y,t) ;
    omega_wb_z_dot = differentiate(c_omega_wb_z,t) ;
    omega_wb_dot = [omega_wb_x_dot,omega_wb_y_dot,omega_wb_z_dot] ;
    
    % get norm
    omega_wb_norm = myNorm(omegaWingBody) ;
    omega_wb_dot_norm = myNorm(omega_wb_dot) ;
    
    % also need derivatives of angle of attack w.r.t. time
    c_alpha = fit(t, alpha, fitType) ;
    [alpha_dot, alpha_ddot] = differentiate(c_alpha, t) ;
    
    % now calculate the two terms in the added mass force (magnitudes)
    fAddedMass1 = (airDensity * (pi/4) * span^2 * mean_chord^2) .*...
        c2_r_int .*(omega_wb_dot_norm .* sin(alpha) + ...
        omega_wb_norm.*alpha_dot.*cos(alpha))  ;
    fAddedMass2 = -1 * alpha_ddot .* airDensity .* (pi/16) .*...
        mean_chord^3 .* span .* c2_int ;
    
    % sum to get total added-mass force
    fAddedMass = fAddedMass1 + fAddedMass2 ;
    
    % added-mass force also acts normal to wing
    F_a = fAddedMass .* wingNormalHat ;
else
    F_a = zeros(size(F_t)) ;
end

% ----------------------------------------
%% sum components to get total force
F_tot = F_t + F_rot + F_a ; 

% --------------------------------------------------------------------
% %% rotate force terms from stroke plane frame to body frame
F_tot = (sp2body*F_tot')' ; 
F_t = (sp2body*F_t')' ; 
F_rot = (sp2body*F_rot')' ; 
F_a = (sp2body*F_a')' ; 

% also rotate lift and drag unit vectors
Lhat = (sp2body*Lhat')' ; 
Dhat = (sp2body*Dhat')' ; 

end